!> ============================================================================
!> relax1D: 1D relaxation-propagation study
!> (The default 'input' file is the original smoke test and doubles as a
!> fast regression check: ~1 s serial run, p-relaxation only.)
!> ============================================================================
!> Purpose: verify, with minimal new code, that the three ingredients needed
!> for the planned 1D Sembian-relaxation study work together inside NGA2's
!> amrmpcomp framework:
!>   (1) quasi-1D domain on a single AMR level (no refinement),
!>   (2) NASG liquid + multi-species ideal-gas mixture (vapor + air), which
!>       adds a transported vapor partial density Q(8) to the state vector,
!>   (3) the relax_igmix_nasg relaxation model coupled to the solver
!>       (mechanical 'p' relaxation for the smoke test; 'pT' and 'pTg' are
!>       selectable from the input file for later stages).
!>
!> Physical setup: closed(ish) 1D tube, liquid water on the left, moist air
!> (vapor mass fraction Yv0) on the right, at different pressures and equal
!> or different temperatures. On the first steps the interfacial cell relaxes
!> toward mechanical equilibrium and the resulting pressure waves propagate
!> into both phases -- exactly the mechanism the full study will quantify.
!>
!> This file is a modified copy of examples/amr_shocktube/src/simulation.f90.
!> Every modification relative to that reference is flagged with "SMOKE:".
!> ============================================================================
module simulation
   use precision,              only: WP
   use string,                 only: str_medium
   use amrgrid_class,          only: amrgrid
   use amrmpcomp_class,        only: amrmpcomp
   use amrviz_class,           only: amrviz
   use amrdata_class,          only: amrdata
   use timetracker_class,      only: timetracker
   use event_class,            only: event
   use monitor_class,          only: monitor
   ! SMOKE: EOS classes change relative to amr_shocktube:
   !   liquid: stiffened_gas -> nasg  (Noble-Abel SG, adds co-volume b)
   !   gas   : ideal_gas -> igmix     (mixture of calorically-perfect species;
   !                                   here ns=2: species 1 = H2O vapor,
   !                                              species 2 = air (closure))
   use nasg_class,             only: nasg
   use igmix_class,            only: igmix
   ! SMOKE: relaxation model changes: relax_ig_sg -> relax_igmix_nasg.
   ! The model enum constants (Prelax/PTrelax/PTgrelax) live in the parent
   ! module relax_igmix_sg_class; the NASG extension re-uses them.
   use relax_igmix_nasg_class, only: relax_igmix_nasg
   use relax_igmix_sg_class,   only: Prelax,PTrelax,PTgrelax
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid (used on a single level here -- Max level: 0 in the input)
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,Umag,Mach

   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile

   !> Materials
   !> SMOKE: water is NASG (not SG); the gas is a 2-species ideal-gas mixture
   type(nasg),  target :: water
   type(igmix), target :: gasmix

   !> Relaxation model
   !> SMOKE: igmix+NASG-aware relaxation; supports p, pT and pTg models
   type(relax_igmix_nasg), target :: relax_model

   !> Species indexing convention for the gas mixture
   !> The solver transports the first (ns-1) gas species in trailing slots of
   !> Q; the LAST species is recovered by closure (Y_last = 1 - sum(others)).
   !> We therefore put VAPOR FIRST (transported, Q(8)) and AIR LAST (closure).
   integer, parameter :: indV=1   !< Vapor index in gasmix
   integer, parameter :: indA=2   !< Air (carrier) index in gasmix

   !> Flow parameters
   real(WP) :: pG_init,TG_init,Yv0     !< Gas state (right side): pressure, temperature, vapor mass fraction
   real(WP) :: pL_init,TL_init         !< Liquid state (left side): pressure, temperature
   real(WP) :: rhoG,rhoL               !< Densities derived from (p,T) through each EOS
   real(WP) :: x_int                   !< Interface position
   real(WP) :: relax_tstart            !< Time after which relaxation is applied (0 = from the start)

   !> Molecular transport parameters (all read from input; zero disables the mechanism)
   real(WP) :: muL,muG                 !< Dynamic viscosities [Pa.s]
   real(WP) :: PrL,PrG                 !< Phasic Prandtl numbers Pr = mu*cp/lambda
   real(WP) :: ScV                     !< Vapor Schmidt number Sc = mu/(rho*D)
   real(WP) :: Cartif                  !< Artificial bulk-viscosity coefficient (0 = off)

   !> Optional second gas state for single-mechanism verification tests:
   !> for x >= x_step the gas is initialized at (pG, TG2, Yv2) instead of
   !> (pG, TG, Yv0). Same pressure on both sides keeps acoustics weak, so a
   !> temperature step (conduction test) or vapor step (diffusion test)
   !> evolves primarily by diffusion and can be compared to the erf solution.
   !> Defaults (TG2=TG, Yv2=Yv0, x_step huge) make this completely inert.
   real(WP) :: TG2,Yv2,x_step
   real(WP) :: rhoG2                   !< Density of gas state 2 from EOS

   !> CSV profile dump counter (shared with viz Output period schedule)
   integer :: nprof=0

   !> Volume fraction of the interfacial mixture cell (VFlo < VF < VFhi).
   !> For a single-interface 1D run there should be exactly one such cell;
   !> if several exist we keep the most mixed (max VF*(1-VF)); if none,
   !> VFmix is set to -1.
   real(WP) :: VFmix=-1.0_WP

contains

   !> Locate the mixed cell(s) and update the module-level VFmix monitor value.
   subroutine update_VFmix()
      use amrmpcomp_class, only: VFlo,VFhi
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: lvl,i,j,k,ierr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      real(WP) :: vf,score
      real(WP), dimension(2) :: inpair,outpair
      ! Local best: (score, VF); score = VF*(1-VF) in [0, 1/4], or -1 if none
      inpair=[-1.0_WP,-1.0_WP]
      lvl=amr%maxlvl
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pVF=>fs%VF%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            vf=pVF(i,j,k,1)
            if (vf.ge.VFlo.and.vf.le.VFhi) then
               score=vf*(1.0_WP-vf)
               if (score.gt.inpair(1)) then
                  inpair(1)=score
                  inpair(2)=vf
               end if
            end if
         end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      ! Rank with the highest score contributes its VF (MPI_MAXLOC on first entry)
      call MPI_ALLREDUCE(inpair,outpair,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,amr%comm,ierr)
      VFmix=outpair(2)
   end subroutine update_VFmix

   !> Gather the 1D cell-centered profiles to rank 0 and write
   !> data/profiles_NNNNNN.csv (plus a time comment in the header).
   !> Called at t=0 and whenever the Output-period event fires, so the CSV
   !> cadence matches the AMReX plotfiles.
   subroutine write_profiles(t)
      use mpi_f08
      use parallel,         only: MPI_REAL_WP
      use string,           only: str_medium
      use filesys,          only: makedir,isdir
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      real(WP), intent(in) :: t
      integer :: i,j,k,lvl,iunit,ierr,nlocal,ntotal
      integer, dimension(:), allocatable :: recvcount,displs
      real(WP), dimension(:,:), allocatable :: local_d,global_d
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOL,pRHOG,pPL,pPG,pYG,pUVW,pTL,pTG
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP) :: x_cc,dx
      real(WP), dimension(9) :: tmp
      character(len=str_medium) :: filename
      integer, parameter :: nvar=9   !< x,VF,PL,PG,TL,TG,Yv,rho_mix,U
      lvl=0 ! Single level, no AMR
      dx=amr%dx(lvl)
      ! First pass: count local cells
      nlocal=0
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         nlocal=nlocal+(bx%hi(1)-bx%lo(1)+1)
      end do
      call amr%mfiter_destroy(mfi)
      allocate(local_d(nvar,nlocal))
      ! Second pass: fill local array
      nlocal=0
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
         pRHOL=>fs%RHOL%mf(lvl)%dataptr(mfi)
         pRHOG=>fs%RHOG%mf(lvl)%dataptr(mfi)
         pPL  =>fs%PL%mf(lvl)%dataptr(mfi)
         pPG  =>fs%PG%mf(lvl)%dataptr(mfi)
         pTL  =>fs%TL%mf(lvl)%dataptr(mfi)
         pTG  =>fs%TG%mf(lvl)%dataptr(mfi)
         pYG  =>fs%Yg%mf(lvl)%dataptr(mfi)
         pUVW =>fs%UVW%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         j=bx%lo(2); k=bx%lo(3)
         do i=bx%lo(1),bx%hi(1)
            x_cc=amr%xlo+(real(i,WP)+0.5_WP)*dx
            nlocal=nlocal+1
            local_d(1,nlocal)=x_cc
            local_d(2,nlocal)=pVF(i,j,k,1)
            local_d(3,nlocal)=pPL(i,j,k,1)
            local_d(4,nlocal)=pPG(i,j,k,1)
            local_d(5,nlocal)=pTL(i,j,k,1)
            local_d(6,nlocal)=pTG(i,j,k,1)
            local_d(7,nlocal)=pYG(i,j,k,1)
            local_d(8,nlocal)=pVF(i,j,k,1)*pRHOL(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*pRHOG(i,j,k,1)
            local_d(9,nlocal)=pUVW(i,j,k,1)
         end do
      end do
      call amr%mfiter_destroy(mfi)
      ! Gather to root
      allocate(recvcount(amr%nproc),displs(amr%nproc))
      call MPI_GATHER(nlocal*nvar,1,MPI_INTEGER,recvcount,1,MPI_INTEGER,0,amr%comm,ierr)
      ntotal=0
      if (amr%amRoot) then
         displs(1)=0
         do i=2,amr%nproc
            displs(i)=displs(i-1)+recvcount(i-1)
         end do
         ntotal=(displs(amr%nproc)+recvcount(amr%nproc))/nvar
      end if
      allocate(global_d(nvar,ntotal))
      call MPI_GATHERV(local_d,nlocal*nvar,MPI_REAL_WP,global_d,recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
      ! Sort by x and write on root
      if (amr%amRoot) then
         do i=2,ntotal
            if (global_d(1,i).lt.global_d(1,i-1)) then
               j=i-1
               do while (j.ge.1)
                  if (global_d(1,j).le.global_d(1,i)) exit
                  j=j-1
               end do
               j=j+1
               tmp=global_d(:,i); global_d(:,j+1:i)=global_d(:,j:i-1); global_d(:,j)=tmp
            end if
         end do
         if (.not.isdir('data')) call makedir('data')
         nprof=nprof+1
         write(filename,'(a,i6.6,a)') 'data/profiles_',nprof,'.csv'
         open(newunit=iunit,file=trim(filename),form='formatted',status='replace',iostat=ierr)
         write(iunit,'(a,es16.8)') '# time = ',t
         write(iunit,'(a)') 'x,VF,PL,PG,TL,TG,Yv,RHOmix,U'
         do i=1,ntotal
            write(iunit,'(8(g0.17,","),g0.17)') global_d(:,i)
         end do
         close(iunit)
      end if
      deallocate(local_d,global_d,recvcount,displs)
   end subroutine write_profiles

   !> Levelset function for planar interface at x = x_int
   !> Returns positive for liquid (x < x_int), negative for gas (x >= x_int)
   function planar_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=x_int-xyz(1)
   end function planar_levelset

   !> Set molecular transport coefficients from the input-file parameters.
   !>
   !> The solver consumes four cell-centered property fields inside get_dQdt:
   !>   fs%visc   mixture dynamic (shear) viscosity  mu     [Pa.s]
   !>   fs%beta   mixture bulk viscosity             beta   [Pa.s]
   !>   fs%diffL  LIQUID thermal conductivity        lam_L  [W/(m.K)]
   !>   fs%diffG  GAS    thermal conductivity        lam_G  [W/(m.K)]
   !>   fs%diffYl LIQUID species diffusivity         rho*D  [kg/(m.s)]
   !>   fs%diffYg GAS    species diffusivity         rho*D  [kg/(m.s)]
   !>
   !> They enter the fluxes as (see amrmpcomp_class::get_dQdt):
   !>   momentum:  tau = mu*(grad u + grad u^T) + (beta - 2/3 mu)*(div u) I
   !>   phasic heat (own-phase aperture-weighted):  q_K = -VF_K * lam_K * grad T_K
   !>   phasic species (Fickian):                   j_K = -VF_K * (rho D)_K * grad Y_K
   !> plus the interdiffusion enthalpy flux sum_k (h_k - h_ns) * j_k in the
   !> phasic energy equation.
   !>
   !> The coefficients are built from the classic dimensionless groups:
   !>   lam_K   = mu_K * cp_K / Pr_K     (Prandtl number,  per phase)
   !>   (rho D) = mu_G / Sc_V            (Schmidt number, vapor in gas)
   !> so setting Pr=0 or Sc=0 in the input cleanly disables that mechanism
   !> (matches the NGA2_SH cavitation-case convention).
   !>
   !> The gas cp is composition-dependent: cp_G = Yv*cpV + (1-Yv)*cpA, taken
   !> from the local transported vapor mass fraction.
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pYg,pVisc,pBeta,pDiffL,pDiffG,pDiffYl,pDiffYg
      real(WP) :: cp_g
      real(WP), parameter :: myeps=1.0e-15_WP
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pVF    =>fs%VF%mf(lvl)%dataptr(mfi)
            pYg    =>fs%Yg%mf(lvl)%dataptr(mfi)
            pVisc  =>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta  =>fs%beta%mf(lvl)%dataptr(mfi)
            pDiffL =>fs%diffL%mf(lvl)%dataptr(mfi)
            pDiffG =>fs%diffG%mf(lvl)%dataptr(mfi)
            pDiffYl=>fs%diffYl%mf(lvl)%dataptr(mfi)
            pDiffYg=>fs%diffYg%mf(lvl)%dataptr(mfi)
            ! Loop over grown tilebox (properties needed in ghost cells for face averaging)
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Mixture shear viscosity: harmonic VF-blend of the phasic
               ! viscosities (harmonic = series coupling, the appropriate
               ! average for a stress transmitted across an interface)
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(muL,myeps)+(1.0_WP-pVF(i,j,k,1))/max(muG,myeps))
               ! No physical bulk viscosity (artificial part added separately)
               pBeta(i,j,k,1)=0.0_WP
               ! Liquid thermal conductivity: lam_L = muL*cpL/PrL
               if (PrL.gt.0.0_WP) then
                  pDiffL(i,j,k,1)=muL*water%cp/PrL
               else
                  pDiffL(i,j,k,1)=0.0_WP
               end if
               ! Gas thermal conductivity: lam_G = muG*cpG(Yv)/PrG with the
               ! local mixture cp from the transported vapor mass fraction
               if (PrG.gt.0.0_WP) then
                  cp_g=pYg(i,j,k,1)*gasmix%cp(indV)+(1.0_WP-pYg(i,j,k,1))*gasmix%cp(indA)
                  pDiffG(i,j,k,1)=muG*cp_g/PrG
               else
                  pDiffG(i,j,k,1)=0.0_WP
               end if
               ! Liquid species diffusivity: single-species liquid, no diffusion
               pDiffYl(i,j,k,1)=0.0_WP
               ! Vapor-in-air mass diffusivity: rho*D = muG/ScV
               if (ScV.gt.0.0_WP) then
                  pDiffYg(i,j,k,1)=muG/ScV
               else
                  pDiffYg(i,j,k,1)=0.0_WP
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback - set Q and VF for the 1D two-phase tube
   !> Liquid on left (x < x_int), gas on right (x >= x_int)
   !>
   !> SMOKE: state vector layout with liq%ns=1 and gas%ns=2 (nQ=8):
   !>   Q(1) = VF*rhoL            liquid partial density
   !>   Q(2) = (1-VF)*rhoG        gas partial density
   !>   Q(3) = VF*rhoL*eL         liquid internal energy density
   !>   Q(4) = (1-VF)*rhoG*eG     gas internal energy density
   !>   Q(5:7) = mixture momentum (zero initially)
   !>   Q(8) = (1-VF)*rhoG*Yv     VAPOR partial density (new vs amr_shocktube)
   !> Air is NOT stored: Ya = 1 - Yv by closure inside the solver/EOS calls.
   subroutine tube_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpcomp_class, only: VFlo
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,IEG,IEG2,x_cc,myrhoG,myIEG,myYv
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Get phasic specific internal energies from the initial (p,rho) pairs.
      ! rhoL/rhoG were derived from (p,T) in simulation_init, so this closes
      ! a thermodynamically consistent loop: e = e(p,rho) is the exact
      ! inverse of the pressure law used by the solver.
      ! SMOKE: the gas EOS call needs the composition vector y=[Yv, Ya].
      IEL =water%get_e_from_p_rho (p=pL_init,rho=rhoL ,y=[1.0_WP])
      IEG =gasmix%get_e_from_p_rho(p=pG_init,rho=rhoG ,y=[Yv0,1.0_WP-Yv0])
      IEG2=gasmix%get_e_from_p_rho(p=pG_init,rho=rhoG2,y=[Yv2,1.0_WP-Yv2])
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Compute VF from planar levelset
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=planar_levelset,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
            ! Store volume fraction
            pVF(i,j,k,1)=myVF
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
            ! Select gas state: state 1 (TG,Yv0) or state 2 (TG2,Yv2) past
            ! x_step (verification tests only; identical states by default)
            x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            if (x_cc.ge.x_step) then
               myrhoG=rhoG2; myIEG=IEG2; myYv=Yv2
            else
               myrhoG=rhoG;  myIEG=IEG;  myYv=Yv0
            end if
            ! Set conserved variables (see layout in the header comment)
            pQ(i,j,k,1)=(       myVF)*rhoL
            pQ(i,j,k,2)=(1.0_WP-myVF)*myrhoG
            pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
            pQ(i,j,k,4)=pQ(i,j,k,2)*myIEG
            pQ(i,j,k,5)=0.0_WP
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
            ! SMOKE: vapor partial density wherever there is gas
            pQ(i,j,k,8)=pQ(i,j,k,2)*myYv
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine tube_init

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log,die
         use string,   only: str_long
         character(len=str_long)   :: message
         character(len=str_medium) :: relaxstr
         real(WP) :: GammaL,PinfL,CvL,bL,qL,qpL     !< NASG liquid parameters
         real(WP) :: GammaV,CvV,qV,qpV              !< Vapor (H2O) ideal-gas parameters
         real(WP) :: GammaA,CvA,qA,qpA              !< Air ideal-gas parameters
         ! ------------------------------------------------------------------
         ! Liquid EoS: Noble-Abel stiffened gas.
         ! q (energy shift) and qp (entropy shift) matter as soon as phase
         ! change (pTg) is on, because they set the saturation curve; we read
         ! them from the input from day one so the same case file drives all
         ! three relaxation models.
         ! ------------------------------------------------------------------
         call param_read('GammaL',GammaL)
         call param_read('PinfL',PinfL)
         call param_read('CvL',CvL)
         call param_read('bL',bL)
         call param_read('qL',qL,default=0.0_WP)
         call param_read('qpL',qpL,default=0.0_WP)
         call water%initialize(gamma=GammaL,pinf=PinfL,b=bL,cv=CvL,q=qL,qp=qpL,name='water')
         ! ------------------------------------------------------------------
         ! Gas EoS: 2-species ideal-gas mixture.
         ! Species 1 = H2O vapor (transported in Q(8)), species 2 = air
         ! (closure species, never stored). Order must match indV/indA above.
         ! ------------------------------------------------------------------
         call param_read('GammaV',GammaV)
         call param_read('CvV',CvV)
         call param_read('qV',qV,default=0.0_WP)
         call param_read('qpV',qpV,default=0.0_WP)
         call param_read('GammaA',GammaA)
         call param_read('CvA',CvA)
         call param_read('qA',qA,default=0.0_WP)
         call param_read('qpA',qpA,default=0.0_WP)
         ! Carrier species is labeled N2 (input GammaA/CvA); historically this
         ! slot was "air". Vapor remains H2O. Order must match indV/indA.
         call gasmix%initialize(gamma=[GammaV,GammaA],cv=[CvV,CvA],q=[qV,qA],qp=[qpV,qpA], &
         &                      species_names=[character(len=str_medium)::'H2O','N2'],name='H2O_N2')
         ! ------------------------------------------------------------------
         ! Initial states: we read (p,T[,Yv]) per side and derive densities
         ! through each phase's EOS, so the initial condition is exactly on
         ! the EOS surface (no spurious startup transient from inconsistent
         ! rho/p/T triplets).
         ! ------------------------------------------------------------------
         call param_read('Gas pressure',pG_init)
         call param_read('Gas temperature',TG_init)
         call param_read('Gas vapor mass fraction',Yv0)
         call param_read('Liquid pressure',pL_init)
         call param_read('Liquid temperature',TL_init)
         rhoL=water%get_rho_from_p_T (p=pL_init,T=TL_init,y=[1.0_WP])
         rhoG=gasmix%get_rho_from_p_T(p=pG_init,T=TG_init,y=[Yv0,1.0_WP-Yv0])
         ! Optional second gas state past x_step (verification tests; inert by default)
         call param_read('Gas step location',x_step,default=1.0e10_WP)
         call param_read('Gas temperature 2',TG2,default=TG_init)
         call param_read('Gas vapor mass fraction 2',Yv2,default=Yv0)
         rhoG2=gasmix%get_rho_from_p_T(p=pG_init,T=TG2,y=[Yv2,1.0_WP-Yv2])
         ! Interface location
         call param_read('Interface location',x_int)
         ! ------------------------------------------------------------------
         ! Relaxation model selection: 'p' (mechanical), 'pT' (+thermal),
         ! 'pTg' (+chemical/phase change). The smoke test uses 'p'.
         ! ------------------------------------------------------------------
         call param_read('Relaxation type',relaxstr,default='p')
         call param_read('Relax start time',relax_tstart,default=0.0_WP)
         ! ------------------------------------------------------------------
         ! Molecular transport: dynamic viscosities plus Prandtl and Schmidt
         ! numbers per phase. Zero disables the mechanism (see the
         ! get_viscosities documentation for the coefficient definitions).
         ! ------------------------------------------------------------------
         call param_read('Liquid viscosity',muL,default=0.0_WP)
         call param_read('Gas viscosity',muG,default=0.0_WP)
         call param_read('Liquid Prandtl number',PrL,default=0.0_WP)
         call param_read('Gas Prandtl number',PrG,default=0.0_WP)
         call param_read('Vapor Schmidt number',ScV,default=0.0_WP)
         call param_read('Artificial viscosity Cartif',Cartif,default=5.0_WP)
         ! Log
         call water%print(); call gasmix%print()
         write(message,'("[Gas]    rhoG=",es12.5," pG=",es12.5," TG=",es12.5," Yv=",es12.5)') rhoG,pG_init,TG_init,Yv0; call log(message)
         if (x_step.lt.1.0e9_WP) then
            write(message,'("[Gas2]   rhoG2=",es12.5," TG2=",es12.5," Yv2=",es12.5," x_step=",es12.5)') &
            &    rhoG2,TG2,Yv2,x_step; call log(message)
         end if
         write(message,'("[Liquid] rhoL=",es12.5," pL=",es12.5," TL=",es12.5)') rhoL,pL_init,TL_init; call log(message)
         write(message,'("[Interface] x_int=",es12.5)') x_int; call log(message)
         write(message,'("[Transport] muL=",es12.5," muG=",es12.5," PrL=",es12.5," PrG=",es12.5," ScV=",es12.5," Cartif=",es12.5)') &
         &    muL,muG,PrL,PrG,ScV,Cartif; call log(message)
         ! Build the relaxation model: liquid=NASG water, gas=igmix, with the
         ! vapor/air species indices matching the gasmix ordering above.
         call relax_model%initialize(liq=water,gas=gasmix,indV=indV,indA=indA)
         select case (trim(adjustl(relaxstr)))
         case ('p');   relax_model%model=Prelax
         case ('pT');  relax_model%model=PTrelax
         case ('pTg'); relax_model%model=PTgrelax
         case default; call die('[relax1D] Unknown Relaxation type: '//trim(relaxstr))
         end select
         write(message,'("[Relax] model=",a," start time=",es12.5)') trim(relaxstr),relax_tstart; call log(message)
      end block init_eos_and_flow

      ! Initialize AMR grid (x extents from input; default preserves the
      ! original [-1.5, 1.5] m shocktube domain)
      create_amrgrid: block
         amr%name='relax1D'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         call param_read('xlo',amr%xlo,default=-1.5_WP)
         call param_read('xhi',amr%xhi,default=+1.5_WP)
         amr%ylo=-1.0_WP; amr%yhi=+1.0_WP
         amr%zlo=-1.0_WP; amr%zhi=+1.0_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         ! AMReX requires each base dimension to be divisible by blocking_factor
         ! (default 8). Shrink nbloc to the largest power-of-two divisor of nx
         ! so cases like nx=100 (blocking 4) work without changing the grid.
         do while (mod(amr%nx,amr%nbloc).ne.0.and.amr%nbloc.gt.1)
            amr%nbloc=amr%nbloc/2
         end do
         ! Enable quasi-1D: set y extent to match one cell at finest level
         if (amr%ny.eq.1) then
            amr%ylo=-0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
            amr%yhi=+0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         end if
         ! Enable quasi-2D/1D: set z extent to match one cell at finest level
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_foextrap,amrex_bc_reflect_odd
         use amrmpcomp_class,  only: BC_REFLECT
         use amrdata_class,    only: interp_face_lin
         use messager,         only: die,log
         use string,           only: str_long
         character(len=str_medium) :: xbc
         character(len=str_long)   :: message
         ! Assign materials and create flow solver.
         ! SMOKE: with gas%ns=2 the solver automatically sizes nQ=8 and
         ! allocates the Yg primitive field + species advection/diffusion
         ! machinery -- no solver-side change needed.
         fs%liq=>water; fs%gas=>gasmix; call fs%initialize(amr=amr,name='relax1D')
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         if (amr%ny.eq.1) fs%interp_vel=interp_face_lin
         ! SMOKE: hand the igmix/NASG relaxation model to the solver; the
         ! solver calls relax_model%apply(dt,VF,Q,Pjump) cell-by-cell inside
         ! fs%apply_relax, on mixture cells at the finest level only.
         fs%relax=>relax_model
         ! Set initial conditions
         fs%user_init=>tube_init
         ! ----------------------------------------------------------------
         ! x-boundary conditions. Default is a CLOSED domain: rigid walls
         ! (U=0) via reflect_odd on normal momentum/velocity, matching
         ! amrcomp_impact. Set "X boundary type: extrapolate" only for
         ! open/Neumann tests (the old foextrap behaviour).
         ! ----------------------------------------------------------------
         if (.not.amr%xper) then
            call param_read('X boundary type',xbc,default='wall')
            select case (trim(adjustl(xbc)))
            case ('wall','closed','reflect')
               ! VOF symmetry at both walls
               fs%lo_bc(1)=BC_REFLECT; fs%hi_bc(1)=BC_REFLECT
               ! Scalars / energies / species: Neumann; normal momentum: odd
               fs%Q%lo_bc(1,:)=amrex_bc_foextrap; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
               fs%Q%lo_bc(1,5)=amrex_bc_reflect_odd; fs%Q%hi_bc(1,5)=amrex_bc_reflect_odd
               ! Normal face velocity U = 0 (anti-symmetric); slip on V,W
               fs%U%lo_bc(1,:)=amrex_bc_reflect_odd; fs%U%hi_bc(1,:)=amrex_bc_reflect_odd
               fs%V%lo_bc(1,:)=amrex_bc_foextrap;    fs%V%hi_bc(1,:)=amrex_bc_foextrap
               fs%W%lo_bc(1,:)=amrex_bc_foextrap;    fs%W%hi_bc(1,:)=amrex_bc_foextrap
               write(message,'("[BC] X boundary type = wall (closed domain, U=0)")'); call log(message)
            case ('extrapolate','open','neumann')
               fs%lo_bc(1)=BC_REFLECT; fs%hi_bc(1)=BC_REFLECT
               fs%Q%lo_bc(1,:)=amrex_bc_foextrap; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
               fs%U%lo_bc(1,:)=amrex_bc_foextrap; fs%U%hi_bc(1,:)=amrex_bc_foextrap
               fs%V%lo_bc(1,:)=amrex_bc_foextrap; fs%V%hi_bc(1,:)=amrex_bc_foextrap
               fs%W%lo_bc(1,:)=amrex_bc_foextrap; fs%W%hi_bc(1,:)=amrex_bc_foextrap
               write(message,'("[BC] X boundary type = extrapolate (open/Neumann)")'); call log(message)
            case default
               call die('[relax1D] Unknown X boundary type (use wall or extrapolate): '//trim(xbc))
            end select
         end if
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         ! SMOKE: dQdt must have fs%nQ (=8) components, not the hardcoded 7
         ! of the original shocktube case, or the vapor slot would be lost.
         call dQdt%initialize(amr,name='dQdt',ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize grid (no regridding for uniform mesh)
      init_grid: block
         ! Fresh start
         call amr%init_from_scratch(time=time%t)
         ! Build PLIC
         call fs%build_plic(time%t)
         call fs%build_subVF()
         ! Initialize primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! Initialize face velocities
         call fs%get_face_velocity()
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         ! Molecular transport coefficients from input (Pr/Sc = 0 disables)
         call get_viscosities()
         ! Optional artificial bulk viscosity for shock capturing
         if (Cartif.gt.0.0_WP) call fs%add_viscartif(dt=time%dt,Cartif=Cartif)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_grid

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'relax1D',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%TL,1,'TL')
         call viz%add_scalar(fs%TG,1,'TG')
         ! SMOKE: vapor mass fraction (component 1 of the Yg primitive field;
         ! there are gas%ns-1 = 1 transported gas species)
         call viz%add_scalar(fs%Yg,1,'Yv')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%C,1,'C')
         call viz%add_surfmesh(fs%smesh,'plic')
         ! Create visualization output event. The same period drives both
         ! AMReX plotfiles and the CSV profile dumps (see write_profiles).
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Always dump the initial state (plotfile + CSV at t = 0).
         ! Do not gate on viz_evt%occurs(): that test can miss t=0 depending
         ! on the current dt vs Output period.
         call viz%write(time=time%t)
         call write_profiles(t=time%t)
      end block create_viz

      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%RHOLmin,'rhoLmin')
         call mfile%add_column(fs%RHOLmax,'rhoLmax')
         call mfile%add_column(fs%PLmin,'PLmin')
         call mfile%add_column(fs%PLmax,'PLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         ! SMOKE: track the vapor mass fraction bounds; for pure 'p'
         ! relaxation Yv must stay at Yv0 in gas cells (no mass transfer),
         ! so any drift here flags an advection/relaxation bug immediately.
         call mfile%add_column(fs%Ygmin(1),'Yvmin')
         call mfile%add_column(fs%Ygmax(1),'Yvmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         ! Single-interface mixture-cell VF (see update_VFmix)
         call update_VFmix()
         call mfile%add_column(VFmix,'VFmix')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid IntEnergy')
         call consfile%add_column(fs%Qint(4),'Gas IntEnergy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         ! SMOKE: integral of the vapor partial density Q(8). With 'p'
         ! relaxation and closed-ish BCs this must be conserved exactly;
         ! under 'pTg' it will change only through evaporation/condensation.
         call consfile%add_column(fs%Qint(8),'Vapor Mass')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%write()
         ! Create timing monitor
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         call tfile%add_column(fs%wtmax_prim,'prim_max')
         call tfile%write()
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   !> The RK2 time loop is copied verbatim from amr_shocktube; the only
   !> difference is the relaxation gate, which uses the input-controlled
   !> relax_tstart instead of a hardcoded 5e-6 s.
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old state
         call fs%store_old()

         ! ======================= RK2 Stage 1: Q*=Q[n]+dt/2*dQdt(t,Q[n]) =======================
         ! Increment Q without pressure gradient
         call fs%get_dQdt(dQdt=dQdt,dt=0.5_WP*time%dt,time=time%tmid)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=0.5_WP*time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         ! Apply relaxation (mechanical/thermal/chemical depending on model)
         if (time%t.ge.relax_tstart) call fs%apply_relax(dt=0.5_WP*time%dt,time=time%tmid)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities and ensure C/F consistency
         call fs%get_face_velocity(); call fs%average_down_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=0.5_WP*time%dt)
         ! Add surface tension term
         call fs%add_surface_tension(scale=0.5_WP*time%dt)
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%tmid)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! ======================= RK2 Stage 2: Q[n+1]=Q[n]+dt*dQdt(t,Q*) =======================
         ! Increment Q without pressure gradient
         call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         ! Apply relaxation (mechanical/thermal/chemical depending on model)
         if (time%t.ge.relax_tstart) call fs%apply_relax(dt=time%dt,time=time%t)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities and ensure C/F consistency
         call fs%get_face_velocity(); call fs%average_down_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=time%dt)
         ! Add surface tension term
         call fs%add_surface_tension(scale=time%dt)
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! ======================================================================================

         ! Compute viscosities (zero molecular transport for the smoke test)
         call get_viscosities()

         ! Add artificial shock-capturing viscosity
         if (Cartif.gt.0.0_WP) call fs%add_viscartif(dt=time%dt,Cartif=Cartif)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization + CSV profiles on the Output-period schedule
         if (viz_evt%occurs()) then
            call viz%write(time%t)
            call write_profiles(t=time%t)
         end if

         ! Perform and output monitoring
         call fs%get_info()
         call update_VFmix()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()

      end do

      ! Always dump the final state (covers the case where Max time is not
      ! an exact multiple of Output period)
      call write_profiles(t=time%t)

   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize materials
      call water%finalize()
      call gasmix%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
   end subroutine simulation_final

end module simulation
