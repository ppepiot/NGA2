module amrscalar_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid,dispatch_fillbc
   use amrdata_class,    only: amrdata
   use amrflux_class,    only: amrflux
   use amrsolver_class,  only: amrsolver
   use amrex_amr_module, only: amrex_multifab,amrex_boxarray,amrex_distromap
   implicit none
   private

   ! Expose type/constructor/methods
   public :: amrscalar

   !> Constant density scalar solver object definition
   type, extends(amrsolver) :: amrscalar

      ! This is our amrgrid
      class(amrgrid), pointer :: amr=>null()

      ! Scalar variable definition
      integer :: nscalar
      character(len=str_medium), dimension(:), allocatable :: SCname

      ! Scalar data containers
      type(amrdata) :: SC
      type(amrdata) :: SCold

      ! Flux register for interlevel conservation
      type(amrflux) :: flux

      ! Overlap size
      integer :: nover=2

      ! Monitoring quantities
      real(WP), dimension(:), allocatable :: SCmax,SCmin,SCint

   contains
      ! Basic procedures
      procedure :: initialize
      procedure :: finalize
      procedure :: get_info
      ! Physics procedures
      procedure :: get_dSCdt
      procedure :: copy2old
      procedure :: reflux_avg_lvl
      procedure :: reflux_avg
      ! Implement amrsolver abstract methods
      procedure :: fillbc => amrscalar_fillbc
      procedure :: on_regrid => amrscalar_on_regrid
   end type amrscalar

contains


   !> Initialization for amrscalar solver
   subroutine initialize(this,amr,nscalar,name)
      use iso_c_binding,    only: c_loc
      use messager,         only: die
      use amrex_amr_module, only: amrex_bc_int_dir,amrex_bc_reflect_even,amrex_interp_cell_cons
      implicit none
      class(amrscalar), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      integer, intent(in) :: nscalar
      character(len=*), optional :: name

      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))

      ! Point to amrgrid object
      this%amr=>amr

      ! Set the number of scalars
      this%nscalar=nscalar
      if (this%nscalar.le.0) call die('[amrscalar initialize] At least 1 scalar is required')

      ! Initialize scalar names
      allocate(this%SCname(1:this%nscalar))
      this%SCname=''

      ! Initialize info storage
      allocate(this%SCmin(1:this%nscalar)); this%SCmin=+huge(1.0_WP)
      allocate(this%SCmax(1:this%nscalar)); this%SCmax=-huge(1.0_WP)
      allocate(this%SCint(1:this%nscalar)); this%SCint= 0.0_WP

      ! Configure SC amrdata
      this%SC%name='SC'
      this%SC%ncomp=this%nscalar
      this%SC%ng=0
      this%SC%interp=amrex_interp_cell_cons
      allocate(this%SC%mf(0:this%amr%nlvl))
      allocate(this%SC%lo_bc(3,this%nscalar),this%SC%hi_bc(3,this%nscalar))
      this%SC%lo_bc=amrex_bc_reflect_even
      if (this%amr%xper) this%SC%lo_bc(1,:)=amrex_bc_int_dir
      if (this%amr%yper) this%SC%lo_bc(2,:)=amrex_bc_int_dir
      if (this%amr%zper) this%SC%lo_bc(3,:)=amrex_bc_int_dir
      this%SC%hi_bc=this%SC%lo_bc

      ! Configure SCold amrdata
      this%SCold%name='SCold'
      this%SCold%ncomp=this%nscalar
      this%SCold%ng=0
      this%SCold%interp=amrex_interp_cell_cons
      allocate(this%SCold%mf(0:this%amr%nlvl))
      allocate(this%SCold%lo_bc(3,this%nscalar),this%SCold%hi_bc(3,this%nscalar))
      this%SCold%lo_bc=this%SC%lo_bc
      this%SCold%hi_bc=this%SC%hi_bc

      ! Configure flux register
      this%flux%name='SC_flux'
      this%flux%ncomp=this%nscalar
      allocate(this%flux%fr(1:this%amr%nlvl))

      ! Register callbacks with amrgrid (pass self as context)
      ! Note: Need select type since c_loc requires non-polymorphic type
      select type (this)
       type is (amrscalar)
         call this%amr%add_on_init(on_init_level,c_loc(this))
         call this%amr%add_on_coarse(on_coarse_level,c_loc(this))
         call this%amr%add_on_remake(on_remake_level,c_loc(this))
         call this%amr%add_on_clear(on_clear_level,c_loc(this))
      end select

   end subroutine initialize


   !> Callback: create level from scratch
   subroutine on_init_level(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx,this)
      ! Allocate SC and SCold
      call this%SC%define(lvl,ba,dm)
      call this%SCold%define(lvl,ba,dm)
      ! Build flux register for fine levels
      if (lvl.ge.1) call this%flux%build_level(lvl,ba,dm,this%amr%rref(lvl-1))
   end subroutine on_init_level


   !> Callback: create level from coarse
   subroutine on_coarse_level(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx,this)
      ! Allocate SC and SCold
      call this%SC%define(lvl,ba,dm)
      call this%SCold%define(lvl,ba,dm)
      ! Fill SC from coarse level
      call this%amr%fill_from_coarse(this%SC,lvl,time)
      ! Build flux register for fine levels
      if (lvl.ge.1) call this%flux%build_level(lvl,ba,dm,this%amr%rref(lvl-1))
   end subroutine on_coarse_level


   !> Callback: remake level (after regrid)
   !> NOTE: This is only called for levels that already exist; new levels go through on_coarse_level
   subroutine on_remake_level(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_ptr,c_f_pointer
      use amrex_amr_module, only: amrex_multifab_build,amrex_multifab_destroy
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrscalar), pointer :: this
      type(amrex_multifab) :: SCnew
      call c_f_pointer(ctx,this)

      ! Step 1: Build temp MultiFab with new geometry
      call amrex_multifab_build(SCnew,ba,dm,this%nscalar,0)

      ! Step 2: Fill temp from old data (old MultiFabs still valid)
      call this%amr%fill_mfab(SCnew,this%SC,lvl,time)

      ! Step 3: Destroy old data
      call amrex_multifab_destroy(this%SC%mf(lvl))
      call amrex_multifab_destroy(this%SCold%mf(lvl))
      if (lvl.ge.1) call this%flux%destroy_level(lvl)

      ! Step 4: Build new data on new geometry
      call amrex_multifab_build(this%SC%mf(lvl),ba,dm,this%nscalar,this%SC%ng)
      call amrex_multifab_build(this%SCold%mf(lvl),ba,dm,this%nscalar,this%SCold%ng)
      if (lvl.ge.1) call this%flux%build_level(lvl,ba,dm,this%amr%rref(lvl-1))

      ! Step 5: Copy from temp to new SC
      call this%SC%mf(lvl)%copy(SCnew,1,1,this%nscalar,0)

      ! Step 6: Destroy temp
      call amrex_multifab_destroy(SCnew)
   end subroutine on_remake_level



   !> Callback: clear level
   subroutine on_clear_level(ctx,lvl)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx,this)
      call this%SC%clear_level(lvl)
      call this%SCold%clear_level(lvl)
      if (lvl.ge.1) call this%flux%destroy_level(lvl)
   end subroutine on_clear_level


   !> Finalization for amrscalar solver
   impure elemental subroutine finalize(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      ! Destroy containers
      call this%SC%destroy()
      call this%SCold%destroy()
      call this%flux%destroy()
      ! Deallocate arrays
      if (allocated(this%SCmax)) deallocate(this%SCmax)
      if (allocated(this%SCmin)) deallocate(this%SCmin)
      if (allocated(this%SCint)) deallocate(this%SCint)
      if (allocated(this%SCname)) deallocate(this%SCname)
      nullify(this%amr)
   end subroutine finalize


   !> Calculate various information on our amrscalar object
   subroutine get_info(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: lvl,nsc
      this%SCmin=+huge(1.0_WP)
      this%SCmax=-huge(1.0_WP)
      this%SCint= 0.0_WP
      do nsc=1,this%nscalar
         do lvl=0,this%amr%clvl()
            this%SCmin(nsc)=min(this%SCmin(nsc),this%SC%mf(lvl)%min(comp=nsc))
            this%SCmax(nsc)=max(this%SCmax(nsc),this%SC%mf(lvl)%max(comp=nsc))
         end do
         this%SCint(nsc)=this%SC%mf(0)%sum(comp=nsc)*(this%amr%dx(0)*this%amr%dy(0)*this%amr%dz(0))/this%amr%vol
      end do
   end subroutine get_info


   !> Calculate dSC/dt at level (lvl)
   subroutine get_dSCdt(this,lvl,dSCdt,SC,U,V,W)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: dSCdt
      type(amrex_multifab), intent(in) :: SC,U,V,W
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab), dimension(3) :: flx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: rhs,pSC,FX,FY,FZ,pU,pV,pW
      integer :: i,j,k,nsc

      ! Prepare flux multifabs
      call this%amr%mfab_build(lvl=lvl,mfab=flx(1),ncomp=this%nscalar,nover=0,atface=[.true. ,.false.,.false.])
      call this%amr%mfab_build(lvl=lvl,mfab=flx(2),ncomp=this%nscalar,nover=0,atface=[.false.,.true. ,.false.])
      call this%amr%mfab_build(lvl=lvl,mfab=flx(3),ncomp=this%nscalar,nover=0,atface=[.false.,.false.,.true. ])

      ! Loop over boxes
      call this%amr%mfiter_build(lvl,mfi)
      do while(mfi%next())
         bx=mfi%tilebox()
         pSC=>SC%dataptr(mfi)
         rhs=>dSCdt%dataptr(mfi)
         pU=>U%dataptr(mfi)
         pV=>V%dataptr(mfi)
         pW=>W%dataptr(mfi)
         FX=>flx(1)%dataptr(mfi)
         FY=>flx(2)%dataptr(mfi)
         FZ=>flx(3)%dataptr(mfi)

         do nsc=1,this%nscalar
            ! X fluxes
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)+1
                     FX(i,j,k,nsc)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i-2,j,k,nsc)+5.0_WP/6.0_WP*pSC(i-1,j,k,nsc)+2.0_WP/6.0_WP*pSC(i  ,j,k,nsc)) &
                     &             -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i-1,j,k,nsc)+5.0_WP/6.0_WP*pSC(i  ,j,k,nsc)-1.0_WP/6.0_WP*pSC(i+1,j,k,nsc))
                  end do
               end do
            end do
            ! Y fluxes
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)+1
                  do i=bx%lo(1),bx%hi(1)
                     FY(i,j,k,nsc)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i,j-2,k,nsc)+5.0_WP/6.0_WP*pSC(i,j-1,k,nsc)+2.0_WP/6.0_WP*pSC(i,j  ,k,nsc)) &
                     &             -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i,j-1,k,nsc)+5.0_WP/6.0_WP*pSC(i,j  ,k,nsc)-1.0_WP/6.0_WP*pSC(i,j+1,k,nsc))
                  end do
               end do
            end do
            ! Z fluxes
            do k=bx%lo(3),bx%hi(3)+1
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)
                     FZ(i,j,k,nsc)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i,j,k-2,nsc)+5.0_WP/6.0_WP*pSC(i,j,k-1,nsc)+2.0_WP/6.0_WP*pSC(i,j,k  ,nsc)) &
                     &             -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i,j,k-1,nsc)+5.0_WP/6.0_WP*pSC(i,j,k  ,nsc)-1.0_WP/6.0_WP*pSC(i,j,k+1,nsc))
                  end do
               end do
            end do
            ! RHS
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)
                     rhs(i,j,k,nsc)=(FX(i+1,j,k,nsc)-FX(i,j,k,nsc))/this%amr%dx(lvl)+&
                     &              (FY(i,j+1,k,nsc)-FY(i,j,k,nsc))/this%amr%dy(lvl)+&
                     &              (FZ(i,j,k+1,nsc)-FZ(i,j,k,nsc))/this%amr%dz(lvl)
                  end do
               end do
            end do
         end do
         ! Scale fluxes for refluxing
         FX=FX*this%amr%dy(lvl)*this%amr%dz(lvl)
         FY=FY*this%amr%dz(lvl)*this%amr%dx(lvl)
         FZ=FZ*this%amr%dx(lvl)*this%amr%dy(lvl)
      end do
      call this%amr%mfiter_destroy(mfi)

      ! Handle refluxing
      if (lvl.gt.0)               call this%flux%fr(lvl  )%fineadd (flx,-1.0_WP)
      if (lvl.lt.this%amr%clvl()) call this%flux%fr(lvl+1)%crseinit(flx,+1.0_WP)

      ! Cleanup
      call this%amr%mfab_destroy(flx(1))
      call this%amr%mfab_destroy(flx(2))
      call this%amr%mfab_destroy(flx(3))

   end subroutine get_dSCdt


   !> Copy SC in SCold
   subroutine copy2old(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: lvl
      do lvl=0,this%amr%clvl()
         call this%SCold%mf(lvl)%copy(srcmf=this%SC%mf(lvl),srccomp=1,dstcomp=1,nc=this%nscalar,ng=0)
      end do
   end subroutine copy2old


   !> Perform refluxing and averaging at a given level
   subroutine reflux_avg_lvl(this,lvl,dt)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer, intent(in) :: lvl
      interface
         subroutine amrex_fi_fluxregister_reflux(fr,mf,scale,geom) bind(c)
            import :: c_ptr,WP
            type(c_ptr), value :: fr,mf,geom
            real(WP), value :: scale
         end subroutine amrex_fi_fluxregister_reflux
      end interface
      call amrex_fi_fluxregister_reflux(this%flux%fr(lvl+1)%p,this%SC%mf(lvl)%p,dt,this%amr%geom(lvl)%p)
      call this%amr%average_downto(this%SC%mf,lvl)
   end subroutine reflux_avg_lvl


   !> Perform successive refluxing and averaging at all levels
   subroutine reflux_avg(this,dt)
      implicit none
      class(amrscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: lvl
      do lvl=this%amr%clvl()-1,0,-1
         call this%reflux_avg_lvl(lvl,dt)
      end do
   end subroutine reflux_avg


   !> Implement amrsolver%fillbc - apply physical BCs to scalar data
   subroutine amrscalar_fillbc(this,mf_ptr,geom_ptr,time,scomp,ncomp)
      use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrscalar), intent(inout) :: this
      type(c_ptr), value :: mf_ptr
      type(c_ptr), value :: geom_ptr
      real(WP), value :: time
      integer, value :: scomp,ncomp
      type(amrex_geometry) :: geom
      type(amrex_multifab) :: mf
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer, dimension(4) :: plo,phi
      ! Skip if fully periodic
      if (all([this%amr%xper,this%amr%yper,this%amr%zper])) return
      ! Convert pointers
      geom=geom_ptr; mf=mf_ptr
      ! Loop over boxes and apply filcc
      call amrex_mfiter_build(mfi,mf)
      do while(mfi%next())
         p=>mf%dataptr(mfi)
         if (.not.geom%domain%contains(p)) then
            plo=lbound(p); phi=ubound(p)
            call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,&
            &   geom%get_physical_location(plo),this%SC%lo_bc,this%SC%hi_bc)
         end if
      end do
   end subroutine amrscalar_fillbc


   !> Implement amrsolver%on_regrid - called after regrid for each level
   subroutine amrscalar_on_regrid(this,lvl,time)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, value :: lvl
      real(WP), value :: time
      ! Default: nothing special needed - rebuild is done in callbacks
   end subroutine amrscalar_on_regrid


end module amrscalar_class
