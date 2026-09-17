!> Coupling of finite-rate chemistry (generated fcmech module) to the AMR compressible solver amrcomp:
!> explicit species source terms evaluated at the solver's current primitive state, transport properties from
!> the simplified fctransport model with Vreman SGS contributions, chemical time-scale monitoring for the time
!> step, and heat-release diagnostics. The dependency is one-way: amrchem uses amrcomp, never the reverse.
!> The kinetics module orders species its own way while the solver carries ns-1 species in the material's order
!> with the last material species as closure; both maps are built by species name at initialization.
module amrchem_class
   use precision,         only: WP
   use string,            only: str_medium
   use messager,          only: die
   use amrgrid_class,     only: amrgrid
   use amrdata_class,     only: amrdata
   use amrcomp_class,     only: amrcomp
   use fctransport_class, only: fctransport
   use fcmech,            only: nS
   implicit none
   private

   public :: amrchem

   real(WP), parameter :: tau_inactive=1.0e10_WP   !< Chemical time scale reported where chemistry is inactive [s]

   type :: amrchem
      character(len=str_medium) :: name='UNNAMED_AMRCHEM'
      class(amrgrid), pointer :: amr=>null()
      type(fctransport) :: tr
      ! Species maps
      integer, dimension(:), allocatable :: kmat   !< kmat(n): fcmech index of material species n
      integer, dimension(:), allocatable :: iQ     !< iQ(k): Q component of fcmech species k (0 for the closure species)
      integer :: kclose=0                          !< fcmech index of the solver's closure species
      ! Settings
      real(WP) :: T_min=0.0_WP           !< No chemistry below this temperature (0: react everywhere)
      real(WP) :: Y_tau_min=1.0e-10_WP   !< Species mass-fraction floor for the chemical time-scale estimate
      real(WP) :: CFLchem_max=1.0_WP     !< Allowed chemical CFL dt*max(D_k); the RK2 limit is 2, accuracy degrades above ~1.5
      ! Diagnostic fields (piecewise-constant interpolation so they survive regrids)
      type(amrdata) :: hrr               !< Heat release rate -sum_k h_k rho dY_k/dt [W/m3]
      type(amrdata) :: tauchem           !< Chemical time scale 1/max_k(D_k), D_k the destruction rate coefficient [s]
      ! Statistics of the last time step (global after get_stats)
      real(WP) :: CFLchem=0.0_WP         !< Chemical CFL dt/tau_chem
      real(WP) :: tau_min=tau_inactive   !< Smallest chemical time scale [s]
      real(WP) :: hrr_int=0.0_WP         !< Integrated heat release rate [W]
      real(WP) :: hrr_max=0.0_WP         !< Maximum heat release rate [W/m3]
      real(WP) :: wtime=0.0_WP           !< Wall time spent in add_source (max over ranks) [s]
      integer  :: ncell_react=0          !< Cell evaluations with active chemistry
      ! Rank-local accumulators
      real(WP), private :: CFLchem_loc=0.0_WP,wtime_loc=0.0_WP
      integer,  private :: ncell_loc=0
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: get_Y
      procedure :: add_source
      procedure :: update_properties
      procedure :: get_stats
      procedure :: print
   end type amrchem

contains

   !> Build the species maps by name and register the diagnostic fields (call before amr%init_from_scratch)
   subroutine initialize(this,amr,fs,name)
      use amrdata_class, only: interp_const
      use fcmech,        only: fcmech_get_speciesnames
      class(amrchem), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      class(amrcomp), intent(in) :: fs
      character(len=*), intent(in), optional :: name
      character(len=str_medium), dimension(nS) :: names
      integer :: n,k
      if (present(name)) this%name=trim(name)
      this%amr=>amr
      if (fs%mat%ns.ne.nS) call die('[amrchem initialize] the material must carry the nS species of the kinetics module')
      ! Species maps by name
      call fcmech_get_speciesnames(names)
      allocate(this%kmat(nS),this%iQ(nS)); this%kmat=0; this%iQ=0
      do n=1,nS
         do k=1,nS
            if (trim(adjustl(fs%mat%species_names(n))).eq.trim(adjustl(names(k)))) then; this%kmat(n)=k; exit; end if
         end do
         if (this%kmat(n).eq.0) call die('[amrchem initialize] species '//trim(fs%mat%species_names(n))// &
         &                              ' is not in the kinetics module')
      end do
      do n=1,nS-1
         if (this%iQ(this%kmat(n)).ne.0) call die('[amrchem initialize] duplicate species name in the material')
         this%iQ(this%kmat(n))=fs%Y_lo+n-1
      end do
      this%kclose=this%kmat(nS)
      if (this%iQ(this%kclose).ne.0) call die('[amrchem initialize] duplicate species name in the material')
      ! Transport model
      call this%tr%initialize(name=this%name)
      ! Diagnostic fields
      call this%hrr%initialize(amr,name='hrr',ncomp=1,ng=0,interp=interp_const); call this%hrr%register()
      call this%tauchem%initialize(amr,name='tauchem',ncomp=1,ng=0,interp=interp_const); call this%tauchem%register()
   end subroutine initialize

   subroutine finalize(this)
      class(amrchem), intent(inout) :: this
      call this%hrr%finalize()
      call this%tauchem%finalize()
      call this%tr%finalize()
      if (allocated(this%kmat)) deallocate(this%kmat)
      if (allocated(this%iQ)) deallocate(this%iQ)
      nullify(this%amr)
   end subroutine finalize

   !> Full composition in fcmech order from the solver's cached (clipped) mass fractions of one cell
   subroutine get_Y(this,Ysolver,Y)
      class(amrchem), intent(in) :: this
      real(WP), dimension(:), intent(in) :: Ysolver     !< ns-1 mass fractions in material order
      real(WP), dimension(nS), intent(out) :: Y
      integer :: n
      do n=1,nS-1
         Y(this%kmat(n))=Ysolver(n)
      end do
      Y(this%kclose)=max(0.0_WP,1.0_WP-sum(Ysolver(1:nS-1)))
   end subroutine get_Y

   !> Add the explicit chemical source rho*dY_k/dt to the species components of dQdt at the solver's current
   !> primitive state (call right after fs%get_dQdt). The internal energy carries the formation enthalpies, so
   !> no energy source is needed. Also fills the heat-release and chemical time-scale diagnostics and
   !> accumulates the chemical CFL dt/tau_chem for get_stats. The time scale is the inverse of the largest
   !> species destruction-rate coefficient D_k = (destruction rate of k)/Y_k, i.e. the magnitude of the chemical
   !> Jacobian diagonal that bounds the stable explicit time step (net rates would miss partial-equilibrium
   !> stiffness and spike when a species crosses the floor).
   subroutine add_source(this,fs,dQdt,dt)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use fcmech,           only: fcmech_get_ydot_ddot
      use parallel,         only: parallel_time
      class(amrchem), intent(inout) :: this
      class(amrcomp), intent(in) :: fs
      type(amrdata), intent(inout) :: dQdt
      real(WP), intent(in) :: dt
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pT,pP,pY,pdQ,pH,pTau
      real(WP), dimension(nS) :: y,ydot,ddot,hk
      real(WP) :: T,rho,tau,Dmax,t0
      integer :: lvl,i,j,k,n
      t0=parallel_time()
      do lvl=0,fs%amr%clvl()
         call fs%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ  =>fs%Q%mf(lvl)%dataptr(mfi)
            pT  =>fs%T%mf(lvl)%dataptr(mfi)
            pP  =>fs%P%mf(lvl)%dataptr(mfi)
            pY  =>fs%Y%mf(lvl)%dataptr(mfi)
            pdQ =>dQdt%mf(lvl)%dataptr(mfi)
            pH  =>this%hrr%mf(lvl)%dataptr(mfi)
            pTau=>this%tauchem%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pH(i,j,k,1)=0.0_WP; pTau(i,j,k,1)=tau_inactive
               T=pT(i,j,k,1)
               if (T.lt.this%T_min) cycle
               rho=pQ(i,j,k,1)
               if (.not.(T.gt.0.0_WP.and.T.lt.huge(1.0_WP).and.pP(i,j,k,1).gt.0.0_WP.and.rho.gt.0.0_WP)) call bad_state()
               call this%get_Y(pY(i,j,k,:),y)
               call fcmech_get_ydot_ddot(pP(i,j,k,1),T,y,ydot,ddot)
               ! Species sources (the closure species receives -sum implicitly)
               do n=1,nS
                  if (this%iQ(n).gt.0) pdQ(i,j,k,this%iQ(n))=pdQ(i,j,k,this%iQ(n))+rho*ydot(n)
               end do
               ! Heat release rate
               call this%tr%get_hk(T,hk)
               pH(i,j,k,1)=-rho*sum(hk(:)*ydot(:))
               ! Chemical time scale: inverse of the largest destruction-rate coefficient among species above the floor
               Dmax=0.0_WP
               do n=1,nS
                  if (y(n).gt.this%Y_tau_min) Dmax=max(Dmax,ddot(n)/y(n))
               end do
               tau=tau_inactive; if (Dmax.gt.0.0_WP) tau=min(tau,1.0_WP/Dmax)
               pTau(i,j,k,1)=tau
               this%CFLchem_loc=max(this%CFLchem_loc,dt/tau)
               this%ncell_loc=this%ncell_loc+1
            end do; end do; end do
         end do
         call fs%amr%mfiter_destroy(mfi)
      end do
      this%wtime_loc=this%wtime_loc+parallel_time()-t0
   contains
      subroutine bad_state()
         use string, only: str_long
         character(len=str_long) :: msg
         write(msg,'(a,i0,a,3(i0,1x),a,es12.5,a,es12.5,a,es12.5,a,es12.5)') '[amrchem add_source] invalid state at level ',lvl, &
         &  ' cell ',i,j,k,': rho=',pQ(i,j,k,1),' rhoe=',pQ(i,j,k,5),' T=',pT(i,j,k,1),' P=',pP(i,j,k,1)
         call die(trim(msg))
      end subroutine bad_state
   end subroutine add_source

   !> Fill the solver's transport properties on grown tiles (ghosts are needed by the face fluxes): molecular
   !> mu, lambda, rho*D from the transport model at the cached primitives, zero bulk viscosity, then optional
   !> Vreman SGS contributions rho*nu_t, rho*nu_t*cp/Pr_t, rho*nu_t/Sc_t and optional artificial bulk viscosity
   subroutine update_properties(this,fs,dt,Cs,Cartif)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use amrsgs,           only: get_vreman
      class(amrchem), intent(inout) :: this
      class(amrcomp), intent(inout) :: fs
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs,Cartif
      type(amrdata) :: cpf
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pT,pY,pVisc,pBeta,pDiff,pDiffY,pCp
      real(WP), dimension(nS) :: y
      real(WP) :: cp,mu,lambda,rhoD
      integer :: lvl,i,j,k
      ! Molecular properties
      call cpf%initialize(fs%amr,name='cpf',ncomp=1,ng=fs%nover); call cpf%reset()
      do lvl=0,fs%amr%clvl()
         call fs%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pT    =>fs%T%mf(lvl)%dataptr(mfi)
            pY    =>fs%Y%mf(lvl)%dataptr(mfi)
            pVisc =>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta =>fs%beta%mf(lvl)%dataptr(mfi)
            pDiff =>fs%diff%mf(lvl)%dataptr(mfi)
            pDiffY=>fs%diffY%mf(lvl)%dataptr(mfi)
            pCp   =>cpf%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               call this%get_Y(pY(i,j,k,:),y)
               call this%tr%get_props(T=pT(i,j,k,1),Y=y,cp=cp,mu=mu,lambda=lambda,rhoD=rhoD)
               pVisc(i,j,k,1)=mu
               pBeta(i,j,k,1)=0.0_WP
               pDiff(i,j,k,1)=lambda
               pDiffY(i,j,k,1)=rhoD
               pCp(i,j,k,1)=cp
            end do; end do; end do
         end do
         call fs%amr%mfiter_destroy(mfi)
      end do
      ! Vreman SGS contributions
      if (present(Cs)) then
         add_sgs: block
            type(amrdata) :: nut
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pNut
            real(WP) :: rho
            call nut%initialize(fs%amr,name='nut',ncomp=1,ng=fs%nover); call nut%reset()
            call get_vreman(dt=dt,visc=nut,U=fs%UVW,V=fs%UVW,W=fs%UVW,Ucomp=1,Vcomp=2,Wcomp=3,Cs=Cs)
            call nut%sync()
            do lvl=0,fs%amr%clvl()
               call fs%amr%mfiter_build(lvl,mfi)
               do while (mfi%next())
                  pQ    =>fs%Q%mf(lvl)%dataptr(mfi)
                  pNut  =>nut%mf(lvl)%dataptr(mfi)
                  pVisc =>fs%visc%mf(lvl)%dataptr(mfi)
                  pDiff =>fs%diff%mf(lvl)%dataptr(mfi)
                  pDiffY=>fs%diffY%mf(lvl)%dataptr(mfi)
                  pCp   =>cpf%mf(lvl)%dataptr(mfi)
                  bx=mfi%growntilebox(fs%nover)
                  do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     rho=max(pQ(i,j,k,1),fs%rho_floor)
                     pVisc(i,j,k,1) =pVisc(i,j,k,1) +rho*pNut(i,j,k,1)
                     pDiff(i,j,k,1) =pDiff(i,j,k,1) +rho*pNut(i,j,k,1)*pCp(i,j,k,1)/this%tr%Pr_t
                     pDiffY(i,j,k,1)=pDiffY(i,j,k,1)+rho*pNut(i,j,k,1)/this%tr%Sc_t
                  end do; end do; end do
               end do
               call fs%amr%mfiter_destroy(mfi)
            end do
            call nut%finalize()
         end block add_sgs
      end if
      call cpf%finalize()
      ! Artificial bulk viscosity (shock capturing)
      if (present(Cartif)) call fs%add_viscartif(dt=dt,Cartif=Cartif)
   end subroutine update_properties

   !> Reduce the chemistry statistics of the last time step across ranks and evaluate the heat-release
   !> integral (composite, fine-masked) and extrema; then reset the rank-local accumulators
   subroutine get_stats(this,fs)
      use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy,amrex_mfiter,amrex_box
      use amrex_interface,  only: amrmask_make_fine
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX,MPI_SUM,MPI_INTEGER
      use parallel,         only: MPI_REAL_WP
      class(amrchem), intent(inout) :: this
      class(amrcomp), intent(in) :: fs
      type(amrex_imultifab) :: mask
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pH
      integer, dimension(:,:,:,:), contiguous, pointer :: pMask
      integer :: lvl,i,j,k,ierr
      ! Reduce accumulators
      this%CFLchem=this%CFLchem_loc; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLchem,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      this%wtime=this%wtime_loc;     call MPI_ALLREDUCE(MPI_IN_PLACE,this%wtime,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      this%ncell_react=this%ncell_loc; call MPI_ALLREDUCE(MPI_IN_PLACE,this%ncell_react,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
      this%CFLchem_loc=0.0_WP; this%wtime_loc=0.0_WP; this%ncell_loc=0
      ! Extrema (parallel reductions inside get_max/get_min)
      this%hrr_max=-huge(1.0_WP); this%tau_min=tau_inactive
      do lvl=0,this%amr%clvl()
         this%hrr_max=max(this%hrr_max,this%hrr%get_max(lvl=lvl))
         this%tau_min=min(this%tau_min,this%tauchem%get_min(lvl=lvl))
      end do
      ! Composite integral of the heat release rate
      this%hrr_int=0.0_WP
      do lvl=0,this%amr%clvl()
         if (lvl.lt.this%amr%clvl()) then
            call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
            call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
         end if
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pH=>this%hrr%mf(lvl)%dataptr(mfi)
            if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (lvl.lt.this%amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
               this%hrr_int=this%hrr_int+pH(i,j,k,1)*this%amr%cell_vol(lvl)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%hrr_int,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
   end subroutine get_stats

   subroutine print(this)
      use messager, only: log
      use string,   only: str_long
      use fcmech,   only: fcmech_get_speciesnames
      class(amrchem), intent(in) :: this
      character(len=str_long) :: msg
      character(len=str_medium), dimension(nS) :: names
      call fcmech_get_speciesnames(names)
      write(msg,'(a,a)') '[amrchem] ',trim(this%name); call log(msg)
      write(msg,'(2x,a,i0,a,a)') 'species = ',nS,'  closure species = ',trim(names(this%kclose)); call log(msg)
      write(msg,'(2x,a,es12.5,a,es12.5,a,es12.5)') 'T_min=',this%T_min,'  Y_tau_min=',this%Y_tau_min, &
      &  '  CFLchem_max=',this%CFLchem_max
      call log(msg)
      call this%tr%print()
   end subroutine print

end module amrchem_class
