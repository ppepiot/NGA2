!> AMRVOF solver class
!> Provides Volume-of-Fluid advection for two-phase flow on amrgrid
!> IRL-free implementation using native cutting geometry
module amrvof_class
   use iso_c_binding,    only: c_ptr,c_loc,c_f_pointer,c_char
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use surfmesh_class,   only: surfmesh
   use amrex_amr_module, only: amrex_multifab,amrex_box,amrex_boxarray,amrex_distromap
   implicit none
   private

   ! Expose type
   public :: amrvof

   ! Boundary condition types
   integer, parameter, public :: BC_LIQ    =1  ! All liquid in ghost
   integer, parameter, public :: BC_GAS    =2  ! All gas in ghost
   integer, parameter, public :: BC_REFLECT=3  ! Symmetry (mirror across boundary)
   integer, parameter, public :: BC_USER   =4  ! User-defined callback

   ! Default parameters for volume fraction
   real(WP), parameter, public :: VFlo=1.0e-12_WP    ! Minimum VF value considered
   real(WP), parameter, public :: VFhi=1.0_WP-VFlo   ! Maximum VF value considered
   real(WP), parameter, public :: vol_eps=1.0e-8_WP  ! Volume epsilon for division by zero

   !> AMRVOF solver type
   type, extends(amrsolver) :: amrvof

      ! User-configurable callbacks
      procedure(vof_init_iface   ), pointer, nopass :: user_vof_init   =>null()
      procedure(vof_tagging_iface), pointer, nopass :: user_vof_tagging=>null()
      procedure(vof_bc_iface     ), pointer, nopass :: user_vof_bc     =>null()

      ! Boundary conditions (per face, only used if direction is non-periodic)
      integer, dimension(3) :: lo_bc=BC_REFLECT
      integer, dimension(3) :: hi_bc=BC_REFLECT

      ! Ghost cell count
      integer :: nover=2

      ! VF: cell-centered amrdata at all levels
      type(amrdata) :: VF,VFold

      ! Barycenters: finest-level-only multifabs, 3 components each
      type(amrex_multifab) :: CL,CLold    ! Liquid barycenter
      type(amrex_multifab) :: CG,CGold    ! Gas barycenter

      ! PLIC: finest-level-only multifab, 4 components (nx, ny, nz, d)
      type(amrex_multifab) :: PLIC,PLICold

      ! Tagging parameters
      integer :: regrid_buffer=10  ! Number of cells to buffer around interface for tagging

      ! Monitoring quantities
      real(WP) :: VFmin=0.0_WP
      real(WP) :: VFmax=0.0_WP
      real(WP) :: VFint=0.0_WP
      
      ! Per-rank timing accumulators (reset by get_info)
      real(WP) :: wt_advance=0.0_WP    !< Full advance_vof
      real(WP) :: wt_sl=0.0_WP         !< SL flux+update
      real(WP) :: wt_plic=0.0_WP       !< Full build_plic
      real(WP) :: wt_plicnet=0.0_WP    !< PLICnet reconstruction loop
      real(WP) :: wt_polygon=0.0_WP    !< Polygon extraction loop
      ! Reduced timing (min/max across ranks)
      real(WP) :: wtmax_advance=0.0_WP, wtmin_advance=0.0_WP
      real(WP) :: wtmax_sl=0.0_WP,      wtmin_sl=0.0_WP
      real(WP) :: wtmax_plic=0.0_WP,    wtmin_plic=0.0_WP
      real(WP) :: wtmax_plicnet=0.0_WP, wtmin_plicnet=0.0_WP
      real(WP) :: wtmax_polygon=0.0_WP, wtmin_polygon=0.0_WP
      ! Load distribution diagnostics
      integer :: nmixed_max=0, nmixed_min=0
      
      ! Surface mesh for visualization
      type(surfmesh) :: smesh

      ! Flag to skip registration with amrgrid in case of inheritance
      logical :: skip_registration=.false.

   contains
      ! Type-bound constructor/destructor
      procedure :: initialize
      procedure :: finalize
      ! Lifecycle callbacks
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      procedure :: tagging
      procedure :: get_cost
      ! Utilities
      procedure :: store_old              !< Copy current state to old state
      procedure :: fill                   !< Fill VF/CL/CG/PLIC ghosts+BCs
      ! VOF-PLIC methods
      procedure :: build_plic             !< Reconstruct PLIC from VF and barycenters
      procedure :: build_plicnet          !< PLICnet reconstruction
      procedure :: build_polygons         !< Build polygons from PLIC planes
      procedure :: reset_moments          !< Recompute VF/barycenters from PLIC
      procedure :: build_band             !< Build transport band around interface
      ! Physics methods
      procedure :: advance_vof            !< Advect VF using staggered or collocated velocity
      procedure :: get_vof_cfl            !< Compute advective CFL at finest level
      ! Print solver info
      procedure :: get_info
      procedure :: print=>amrvof_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrvof

   !> Abstract interface for user-overridable on_init callback
   abstract interface
      subroutine vof_init_iface(solver,lvl,time,ba,dm)
         import :: amrvof,WP,amrex_boxarray,amrex_distromap
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine vof_init_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine vof_tagging_iface(solver,lvl,time,tags)
         import :: amrvof,c_ptr,WP
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine vof_tagging_iface
   end interface

   !> Abstract interface for user-defined boundary conditions
   !> User sets VF (and CL/CG/PLIC if associated) in ghost cells
   abstract interface
      subroutine vof_bc_iface(solver,lvl,time,face,bx,pVF,pCL,pCG,pPLIC)
         import :: amrvof,amrex_box,WP
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl                  !< Current level
         real(WP), intent(in) :: time
         integer, intent(in) :: face                 !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         type(amrex_box), intent(in) :: bx           !< Ghost region to fill
         real(WP), dimension(:,:,:,:), contiguous, pointer, intent(inout) :: pVF,pCL,pCG,pPLIC
      end subroutine vof_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrvof type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrvof_on_init(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_init(lvl,time,ba,dm)
      if (associated(this%user_vof_init)) call this%user_vof_init(this,lvl,time,ba,dm)
   end subroutine amrvof_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrvof_on_coarse(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_coarse(lvl,time,ba,dm)
   end subroutine amrvof_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrvof_on_remake(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_remake(lvl,time,ba,dm)
   end subroutine amrvof_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrvof_on_clear(ctx,lvl)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrvof), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrvof_on_clear

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrvof_postregrid(ctx,lbase,time)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrvof), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrvof_postregrid

   !> Dispatch tagging: calls type-bound method then user callback
   subroutine amrvof_tagging(ctx,lvl,time,tags)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrvof), pointer :: this
      call c_f_pointer(ctx,this)
      call this%tagging(lvl,time,tags)
      if (associated(this%user_vof_tagging)) call this%user_vof_tagging(this,lvl,time,tags)
   end subroutine amrvof_tagging

   !> Dispatch cost: calls type-bound method
   subroutine amrvof_get_cost(ctx,lvl,nboxes,costs,ba)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrvof), pointer :: this
      call c_f_pointer(ctx,this)
      call this%get_cost(lvl,nboxes,costs,ba)
   end subroutine amrvof_get_cost

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the VOF solver
   subroutine initialize(this,amr,name)
      use amrdata_class, only: amrex_interp_pc
      implicit none
      class(amrvof), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name
      ! Set name
      if (present(name)) then
         this%name=trim(name)
      else
         this%name='UNNAMED_VOF'
      end if
      ! Store amrgrid pointer
      this%amr=>amr
      ! Initialize VF/VFold as amrdata (all levels, used for tagging + average-down)
      call this%VF%initialize   (amr,name='VF'   ,ncomp=1,ng=this%nover,interp=amrex_interp_pc); this%VF%parent   =>this
      call this%VFold%initialize(amr,name='VFold',ncomp=1,ng=this%nover,interp=amrex_interp_pc); this%VFold%parent=>this
      ! Initialize surface mesh for visualization
      this%smesh%name=trim(this%name)//'_plic'
      ! Register callbacks with amrgrid
      if (.not.this%skip_registration) then
         select type (this)
          type is (amrvof)
            call this%amr%add_on_init   (amrvof_on_init,   c_loc(this))
            call this%amr%add_on_coarse (amrvof_on_coarse, c_loc(this))
            call this%amr%add_on_remake (amrvof_on_remake, c_loc(this))
            call this%amr%add_on_clear  (amrvof_on_clear,  c_loc(this))
            call this%amr%add_tagging   (amrvof_tagging,   c_loc(this))
            call this%amr%add_postregrid(amrvof_postregrid,c_loc(this))
            call this%amr%set_get_cost  (amrvof_get_cost,  c_loc(this))
         end select
      end if
      ! Print solver info
      call this%print()
   end subroutine initialize

   !> Finalize the VOF solver
   subroutine finalize(this)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrvof), intent(inout) :: this
      ! Finalize VF/VFold amrdata
      call this%VF%finalize()
      call this%VFold%finalize()
      ! Destroy finest-level multifabs
      call amrex_multifab_destroy(this%CL)
      call amrex_multifab_destroy(this%CG)
      call amrex_multifab_destroy(this%PLIC)
      call amrex_multifab_destroy(this%CLold)
      call amrex_multifab_destroy(this%CGold)
      call amrex_multifab_destroy(this%PLICold)
      ! Reset surface mesh
      call this%smesh%reset()
      ! Nullify pointers
      nullify(this%amr)
      nullify(this%user_vof_init)
      nullify(this%user_vof_tagging)
      nullify(this%user_vof_bc)
   end subroutine finalize

   ! ============================================================================
   ! LIFECYCLE CALLBACKS
   ! ============================================================================

   !> Override on_init: reset VF/VFold level layout, rebuild finest-level mfabs
   subroutine on_init(this,lvl,time,ba,dm)
      use amrgrid_class, only: mfab_rebuild
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset VF/VFold level layout and zero
      call this%VF%reset_level(lvl,ba,dm); call this%VF%setval(val=0.0_WP,lvl=lvl)
      call this%VFold%reset_level(lvl,ba,dm); call this%VFold%setval(val=0.0_WP,lvl=lvl)
      ! Build finest-level multifabs if we are at the finest level
      if (lvl.eq.this%amr%maxlvl) then
         call mfab_rebuild(this%CL     ,ba,dm,nc=3,ng=this%nover)
         call mfab_rebuild(this%CG     ,ba,dm,nc=3,ng=this%nover)
         call mfab_rebuild(this%PLIC   ,ba,dm,nc=4,ng=this%nover)
         call mfab_rebuild(this%CLold  ,ba,dm,nc=3,ng=this%nover)
         call mfab_rebuild(this%CGold  ,ba,dm,nc=3,ng=this%nover)
         call mfab_rebuild(this%PLICold,ba,dm,nc=4,ng=this%nover)
      end if
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse, rebuild finest-level mfabs
   subroutine on_coarse(this,lvl,time,ba,dm)
      use amrgrid_class, only: mfab_rebuild
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Interpolate VF from coarse level and reset VFold
      call this%VF%on_coarse(this%VF,lvl,time,ba,dm)
      call this%VFold%reset_level(lvl,ba,dm)
      ! Build finest-level multifabs if we are at the finest level
      if (lvl.eq.this%amr%maxlvl) then
         ! Rebuild and zero out
         call mfab_rebuild(this%CL,     ba,dm,nc=3,ng=this%nover)
         call mfab_rebuild(this%CG,     ba,dm,nc=3,ng=this%nover)
         call mfab_rebuild(this%PLIC,   ba,dm,nc=4,ng=this%nover)
         call mfab_rebuild(this%CLold,  ba,dm,nc=3,ng=this%nover)
         call mfab_rebuild(this%CGold,  ba,dm,nc=3,ng=this%nover)
         call mfab_rebuild(this%PLICold,ba,dm,nc=4,ng=this%nover)
         ! Set to trivial values
         trivialize: block
            use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
            real(WP) :: dx,dy,dz
            integer :: i,j,k
            dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
            call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
            do while (mfi%next())
               ! Get pointers to data
               pVF  =>this%VF%mf(lvl)%dataptr(mfi)
               pCL  =>this%CL%dataptr(mfi)
               pCG  =>this%CG%dataptr(mfi)
               pPLIC=>this%PLIC%dataptr(mfi)
               ! Get tilebox
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pCL(i,j,k,:)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pCG(i,j,k,:)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,sign(1.0e10_WP,pVF(i,j,k,1)-0.5_WP)]
               end do; end do; end do
            end do
            call amrex_mfiter_destroy(mfi)
         end block trivialize
      end if
   end subroutine on_coarse

   !> Override on_remake: migrate VF on regrid, rebuild finest-level mfabs and fill with parallel_copy
   subroutine on_remake(this,lvl,time,ba,dm)
      use amrgrid_class, only: mfab_rebuild
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Remake VF from existing+coarse data and reset VFold
      call this%VF%on_remake(this%VF,lvl,time,ba,dm)
      call this%VFold%reset_level(lvl,ba,dm)
      ! Rebuild finest-level multifabs if we are at the finest level
      if (lvl.eq.this%amr%maxlvl) then
         remake_finest: block
            use amrex_amr_module, only: amrex_multifab_build,amrex_multifab_destroy, &
            &                           amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
            type(amrex_multifab) :: CL_new,CG_new,PLIC_new
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
            real(WP) :: dx,dy,dz
            integer :: i,j,k
            ! Build new mfabs
            call amrex_multifab_build(CL_new  ,ba,dm,nc=3,ng=this%nover)
            call amrex_multifab_build(CG_new  ,ba,dm,nc=3,ng=this%nover)
            call amrex_multifab_build(PLIC_new,ba,dm,nc=4,ng=this%nover)
            ! Set to trivial values
            dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
            call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
            do while (mfi%next())
               ! Get pointers to data
               pVF  =>this%VF%mf(lvl)%dataptr(mfi)
               pCL  =>CL_new%dataptr(mfi)
               pCG  =>CG_new%dataptr(mfi)
               pPLIC=>PLIC_new%dataptr(mfi)
               ! Get tilebox
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pCL(i,j,k,:)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pCG(i,j,k,:)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,sign(1.0e10_WP,pVF(i,j,k,1)-0.5_WP)]
               end do; end do; end do
            end do
            call amrex_mfiter_destroy(mfi)
            ! Parallel copy surviving data from old to new
            call CL_new%parallel_copy(this%CL,this%amr%geom(lvl))
            call CG_new%parallel_copy(this%CG,this%amr%geom(lvl))
            call PLIC_new%parallel_copy(this%PLIC,this%amr%geom(lvl))
            ! Destroy old, assign new
            call amrex_multifab_destroy(this%CL  ); call this%CL%move(CL_new)
            call amrex_multifab_destroy(this%CG  ); call this%CG%move(CG_new)
            call amrex_multifab_destroy(this%PLIC); call this%PLIC%move(PLIC_new)
            ! Rebuild old multifabs
            call mfab_rebuild(this%CLold,  ba,dm,nc=3,ng=this%nover)
            call mfab_rebuild(this%CGold,  ba,dm,nc=3,ng=this%nover)
            call mfab_rebuild(this%PLICold,ba,dm,nc=4,ng=this%nover)
         end block remake_finest
      end if
   end subroutine on_remake

   !> Override on_clear: clear VF/VFold at level, destroy finest-level mfabs
   subroutine on_clear(this,lvl)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%VF%clear_level(lvl)
      call this%VFold%clear_level(lvl)
      if (lvl.eq.this%amr%maxlvl) then
         call amrex_multifab_destroy(this%CL)
         call amrex_multifab_destroy(this%CG)
         call amrex_multifab_destroy(this%PLIC)
         call amrex_multifab_destroy(this%CLold)
         call amrex_multifab_destroy(this%CGold)
         call amrex_multifab_destroy(this%PLICold)
      end if
   end subroutine on_clear

   !> Override post_regrid: average down VF
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      integer :: lvl
      call this%VF%average_down(lbase)
      do lvl=this%amr%clvl(),lbase,-1
         call this%fill(lvl,time)
      end do
   end subroutine post_regrid

   !> Tag cells near interface with regrid_buffer layer growth
   subroutine tagging(this,lvl,time,tags)
      use amrex_amr_module, only: amrex_tagboxarray,amrex_mfiter
      use amrgrid_class, only: SETtag
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrex_tagboxarray) :: tba
      type(amrex_multifab) :: band
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand
      integer :: i,j,k
      integer :: eff_buffer
      ! Resolve tagboxarray pointer
      tba=tags
      ! Effective buffer size
      eff_buffer=max(1,this%regrid_buffer/(2**(this%amr%maxlvl-lvl)))
      ! Build band of size eff_buffer
      call this%amr%mfab_build(lvl=lvl,mfab=band,ncomp=1,nover=1)
      call this%build_band(lvl=lvl,VF=this%VF%mf(lvl),band=band,nband=eff_buffer)
      ! Set tags from band
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         tagarr=>tba%dataPtr(mfi)
         pBand =>band%dataptr(mfi)
         ! Get interior tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            if (pBand(i,j,k,1).gt.0.0_WP) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      ! Cleanup
      call this%amr%mfab_destroy(band)
   end subroutine tagging

   !> Estimate per-box costs for load balancing
   !> Cost based on number of mixed cells vs pure cells
   subroutine get_cost(this,lvl,nboxes,costs,ba)
      use iso_c_binding, only: c_associated
      use amrex_amr_module, only: amrex_boxarray,amrex_box,amrex_mfiter,&
      &                           amrex_mfiter_build,amrex_mfiter_destroy,amrex_intersection
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: old_bx,new_bx,isect
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer :: n,i,j,k,ierr
      ! Guard: if VF data doesn't exist yet, return uniform costs
      if (.not.allocated(this%VF%mf)) then; costs=1.0_WP; return; end if
      if (.not.c_associated(this%VF%mf(lvl)%p)) then; costs=1.0_WP; return; end if
      ! Coarser levels don't contribute
      if (lvl.lt.this%amr%maxlvl) then
         costs=1.0_WP; return
      end if
      ! At finest level, count mixed cells per new box from local old data
      costs=0.0_WP
      call amrex_mfiter_build(mfi,this%VF%mf(lvl),tiling=.false.)
      do while (mfi%next())
         old_bx=mfi%tilebox()
         pVF=>this%VF%mf(lvl)%dataptr(mfi)
         do n=1,nboxes
            new_bx=ba%get_box(n-1)  ! 0-indexed
            if (.not.old_bx%intersects(new_bx)) cycle
            isect=amrex_intersection(old_bx,new_bx)
            do k=isect%lo(3),isect%hi(3); do j=isect%lo(2),isect%hi(2); do i=isect%lo(1),isect%hi(1)
               if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) costs(n)=costs(n)+1.0_WP
            end do; end do; end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
      ! Allreduce: sum partial mixed-cell counts across ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE,costs,nboxes,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
   end subroutine get_cost

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Copy current state to old state
   subroutine store_old(this)
      implicit none
      class(amrvof), intent(inout) :: this
      ! Copy VF to VFold
      call this%VFold%copy(src=this%VF)
      ! Copy CL, CG, and PLIC
      call this%CLold%copy  (srcmf=this%CL  ,srccomp=1,dstcomp=1,nc=3,ng=this%nover)
      call this%CGold%copy  (srcmf=this%CG  ,srccomp=1,dstcomp=1,nc=3,ng=this%nover)
      call this%PLICold%copy(srcmf=this%PLIC,srccomp=1,dstcomp=1,nc=4,ng=this%nover)
   end subroutine store_old

   !> Unified ghost fill: VF at all levels, CL/CG/PLIC at maxlvl
   subroutine fill(this,lvl,time)
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in), optional :: lvl
      real(WP), intent(in) :: time
      integer :: l
      if (present(lvl)) then
         call fill_lvl(lvl)
      else
         do l=0,this%amr%clvl()
            call fill_lvl(l)
         end do
      end if
   contains
      !> Fill at a level
      subroutine fill_lvl(lvl)
         use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
         implicit none
         integer, intent(in) :: lvl
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: vbx,bc_bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
         integer :: ig,jg,kg,ilo,ihi,jlo,jhi,klo,khi
         integer :: dlo(3),dhi(3),dir,side,face
         integer :: i1,i2,j1,j2,k1,k2,bc_type
         integer, dimension(3) :: ind
         real(WP) :: Lx,Ly,Lz,bnd(3,2)
         logical :: at_finest
         ! Get domain info
         dlo=this%amr%geom(lvl)%domain%lo
         dhi=this%amr%geom(lvl)%domain%hi
         Lx=this%amr%xhi-this%amr%xlo
         Ly=this%amr%yhi-this%amr%ylo
         Lz=this%amr%zhi-this%amr%zlo
         bnd(:,1)=[this%amr%xlo,this%amr%ylo,this%amr%zlo]
         bnd(:,2)=[this%amr%xhi,this%amr%yhi,this%amr%zhi]
         at_finest=lvl.eq.this%amr%maxlvl
         ! Step 1a: validextrap + fill_boundary for VF
         call this%amr%mfab_validextrap(lvl,this%VF%mf(lvl))
         call this%VF%mf(lvl)%fill_boundary(this%amr%geom(lvl))
         ! Step 1b: same for CL/CG/PLIC at finest
         if (at_finest) then
            call this%amr%mfab_validextrap(lvl,this%CL)
            call this%amr%mfab_validextrap(lvl,this%CG)
            call this%amr%mfab_validextrap(lvl,this%PLIC)
            call this%CL%fill_boundary(this%amr%geom(lvl))
            call this%CG%fill_boundary(this%amr%geom(lvl))
            call this%PLIC%fill_boundary(this%amr%geom(lvl))
         end if
         ! Step 2: periodic shift + BCs + trivial ghosts
         call amrex_mfiter_build(mfi,this%VF%mf(lvl),tiling=.false.)
         do while (mfi%next())
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pCL=>null()
            pCG=>null()
            pPLIC=>null()
            vbx=mfi%validbox()
            ilo=lbound(pVF,1); ihi=ubound(pVF,1)
            jlo=lbound(pVF,2); jhi=ubound(pVF,2)
            klo=lbound(pVF,3); khi=ubound(pVF,3)
            if (at_finest) then
               pCL  =>this%CL%dataptr(mfi)
               pCG  =>this%CG%dataptr(mfi)
               pPLIC=>this%PLIC%dataptr(mfi)
               ! Periodic coordinate shift for CL/CG/PLIC
               if (this%amr%xper) then
                  if (ilo.lt.dlo(1)) then
                     do kg=klo,khi; do jg=jlo,jhi; do ig=ilo,dlo(1)-1
                        pCL(ig,jg,kg,1)=pCL(ig,jg,kg,1)-Lx; pCG(ig,jg,kg,1)=pCG(ig,jg,kg,1)-Lx
                        pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)-pPLIC(ig,jg,kg,1)*Lx
                     end do; end do; end do
                  end if
                  if (ihi.gt.dhi(1)) then
                     do kg=klo,khi; do jg=jlo,jhi; do ig=dhi(1)+1,ihi
                        pCL(ig,jg,kg,1)=pCL(ig,jg,kg,1)+Lx; pCG(ig,jg,kg,1)=pCG(ig,jg,kg,1)+Lx
                        pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)+pPLIC(ig,jg,kg,1)*Lx
                     end do; end do; end do
                  end if
               end if
               if (this%amr%yper) then
                  if (jlo.lt.dlo(2)) then
                     do kg=klo,khi; do jg=jlo,dlo(2)-1; do ig=ilo,ihi
                        pCL(ig,jg,kg,2)=pCL(ig,jg,kg,2)-Ly; pCG(ig,jg,kg,2)=pCG(ig,jg,kg,2)-Ly
                        pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)-pPLIC(ig,jg,kg,2)*Ly
                     end do; end do; end do
                  end if
                  if (jhi.gt.dhi(2)) then
                     do kg=klo,khi; do jg=dhi(2)+1,jhi; do ig=ilo,ihi
                        pCL(ig,jg,kg,2)=pCL(ig,jg,kg,2)+Ly; pCG(ig,jg,kg,2)=pCG(ig,jg,kg,2)+Ly
                        pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)+pPLIC(ig,jg,kg,2)*Ly
                     end do; end do; end do
                  end if
               end if
               if (this%amr%zper) then
                  if (klo.lt.dlo(3)) then
                     do kg=klo,dlo(3)-1; do jg=jlo,jhi; do ig=ilo,ihi
                        pCL(ig,jg,kg,3)=pCL(ig,jg,kg,3)-Lz; pCG(ig,jg,kg,3)=pCG(ig,jg,kg,3)-Lz
                        pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)-pPLIC(ig,jg,kg,3)*Lz
                     end do; end do; end do
                  end if
                  if (khi.gt.dhi(3)) then
                     do kg=dhi(3)+1,khi; do jg=jlo,jhi; do ig=ilo,ihi
                        pCL(ig,jg,kg,3)=pCL(ig,jg,kg,3)+Lz; pCG(ig,jg,kg,3)=pCG(ig,jg,kg,3)+Lz
                        pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)+pPLIC(ig,jg,kg,3)*Lz
                     end do; end do; end do
                  end if
               end if
            end if
            ! Physical BCs
            do dir=1,3; do side=-1,1,2
               ! Skip periodic directions
               if (dir.eq.1.and.this%amr%xper) cycle
               if (dir.eq.2.and.this%amr%yper) cycle
               if (dir.eq.3.and.this%amr%zper) cycle
               ! Get BC type
               if (side.eq.-1) then; bc_type=this%lo_bc(dir)
               else;                 bc_type=this%hi_bc(dir)
               end if
               ! Get ghost cell bounds
               i1=ilo; i2=ihi; j1=jlo; j2=jhi; k1=klo; k2=khi
               select case(dir)
                case(1)
                  if (side.eq.-1) then; i2=dlo(1)-1; if (i2.lt.i1) cycle
                  else;                 i1=dhi(1)+1; if (i2.lt.i1) cycle
                  end if
                case(2)
                  if (side.eq.-1) then; j2=dlo(2)-1; if (j2.lt.j1) cycle
                  else;                 j1=dhi(2)+1; if (j2.lt.j1) cycle
                  end if
                case(3)
                  if (side.eq.-1) then; k2=dlo(3)-1; if (k2.lt.k1) cycle
                  else;                 k1=dhi(3)+1; if (k2.lt.k1) cycle
                  end if
               end select
               ! Apply BCs
               select case(bc_type)
                case(BC_REFLECT)
                  do kg=k1,k2; do jg=j1,j2; do ig=i1,i2
                     ind=[ig,jg,kg]
                     if (side.eq.-1) then; ind(dir)=2*dlo(dir)-1-ind(dir)
                     else;                 ind(dir)=2*dhi(dir)+1-ind(dir)
                     end if
                     pVF(ig,jg,kg,1)=pVF(ind(1),ind(2),ind(3),1)
                     if (at_finest) then
                        pCL(ig,jg,kg,:)=pCL(ind(1),ind(2),ind(3),:)
                        pCG(ig,jg,kg,:)=pCG(ind(1),ind(2),ind(3),:)
                        pCL(ig,jg,kg,dir)=2.0_WP*bnd(dir,(3+side)/2)-pCL(ig,jg,kg,dir)
                        pCG(ig,jg,kg,dir)=2.0_WP*bnd(dir,(3+side)/2)-pCG(ig,jg,kg,dir)
                        pPLIC(ig,jg,kg,:)=pPLIC(ind(1),ind(2),ind(3),:)
                        pPLIC(ig,jg,kg,dir)=-pPLIC(ig,jg,kg,dir)
                        pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)-2.0_WP*pPLIC(ind(1),ind(2),ind(3),dir)*bnd(dir,(3+side)/2)
                     end if
                  end do; end do; end do
                case(BC_LIQ)
                  do kg=k1,k2; do jg=j1,j2; do ig=i1,i2
                     pVF(ig,jg,kg,1)=1.0_WP
                     if (at_finest) then
                        pCL(ig,jg,kg,:)=[this%amr%xlo+(real(ig,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(jg,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(kg,WP)+0.5_WP)*this%amr%dz(lvl)]
                        pCG(ig,jg,kg,:)=[this%amr%xlo+(real(ig,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(jg,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(kg,WP)+0.5_WP)*this%amr%dz(lvl)]
                        pPLIC(ig,jg,kg,:)=[0.0_WP,0.0_WP,0.0_WP,+1.0e10_WP]
                     end if
                  end do; end do; end do
                case(BC_GAS)
                  do kg=k1,k2; do jg=j1,j2; do ig=i1,i2
                     pVF(ig,jg,kg,1)=0.0_WP
                     if (at_finest) then
                        pCL(ig,jg,kg,:)=[this%amr%xlo+(real(ig,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(jg,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(kg,WP)+0.5_WP)*this%amr%dz(lvl)]
                        pCG(ig,jg,kg,:)=[this%amr%xlo+(real(ig,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(jg,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(kg,WP)+0.5_WP)*this%amr%dz(lvl)]
                        pPLIC(ig,jg,kg,:)=[0.0_WP,0.0_WP,0.0_WP,-1.0e10_WP]
                     end if
                  end do; end do; end do
                case(BC_USER)
                  if (associated(this%user_vof_bc)) then
                     bc_bx=amrex_box([i1,j1,k1],[i2,j2,k2])
                     face=2*dir-1+(1+side)/2
                     call this%user_vof_bc(solver=this,lvl=lvl,time=time,face=face,bx=bc_bx,pVF=pVF,pCL=pCL,pCG=pCG,pPLIC=pPLIC)
                  end if
               end select
            end do; end do
            ! Trivialize pure-cell ghosts at finest
            if (at_finest) then
               do kg=klo,khi; do jg=jlo,jhi; do ig=ilo,ihi
                  ! Skip interior cells
                  if (ig.ge.vbx%lo(1).and.ig.le.vbx%hi(1).and.jg.ge.vbx%lo(2).and.jg.le.vbx%hi(2).and.kg.ge.vbx%lo(3).and.kg.le.vbx%hi(3)) cycle
                  ! Check if cell is pure liquid or pure gas
                  if (pVF(ig,jg,kg,1).lt.VFlo.or.pVF(ig,jg,kg,1).gt.VFhi) then
                     pCL(ig,jg,kg,:)=[this%amr%xlo+(real(ig,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(jg,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(kg,WP)+0.5_WP)*this%amr%dz(lvl)]
                     pCG(ig,jg,kg,:)=[this%amr%xlo+(real(ig,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(jg,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(kg,WP)+0.5_WP)*this%amr%dz(lvl)]
                     pPLIC(ig,jg,kg,:)=[0.0_WP,0.0_WP,0.0_WP,sign(1.0e10_WP,pVF(ig,jg,kg,1)-0.5_WP)]
                  end if
               end do; end do; end do
            end if
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine fill_lvl
   end subroutine fill

   ! ============================================================================
   ! VOF-PLIC METHODS
   ! ============================================================================

   !> Build PLIC reconstruction from VF and barycenters using PLICnet
   subroutine build_plic(this,time)
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrvof), intent(inout) :: this
      real(WP), intent(in) :: time
      real(WP) :: t0
      integer :: lvl
      ! Start timer
      t0=MPI_Wtime()

      ! Perform PLICnet reconstruction
      call this%build_plicnet(time)

      ! Build polygons from PLIC
      call this%build_polygons()

      ! Reset moments from PLIC
      call this%reset_moments()

      ! Average down to coarse levels
      call this%VF%average_down()
      do lvl=this%amr%clvl()-1,0,-1
         call this%fill(lvl,time)
      end do

      ! End timer
      this%wt_plic=this%wt_plic+(MPI_Wtime()-t0)
   end subroutine build_plic

   !> PLICnet interface reconstruction
   subroutine build_plicnet(this,time)
      use plicnet, only: get_normal,reflect_moments
      use mathtools, only: normalize
      use amrvof_geometry, only: get_plane_dist
      use amrex_amr_module, only: amrex_mfiter
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrvof), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl,i,j,k,ii,jj,kk,direction,direction2
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      real(WP), dimension(0:188) :: moments
      real(WP), dimension(3) :: normal,center,lo,hi
      real(WP) :: m000,m100,m010,m001,temp,t0
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
      logical :: flip
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      ! Start timer
      t0=MPI_Wtime()
      ! Only build at finest level
      lvl=this%amr%maxlvl
      ! Get cell size at this level
      dx=this%amr%dx(lvl); dxi=1.0_WP/dx
      dy=this%amr%dy(lvl); dyi=1.0_WP/dy
      dz=this%amr%dz(lvl); dzi=1.0_WP/dz
      ! Perform PLICnet reconstruction
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers
         pVF=>this%VF%mf(lvl)%dataptr(mfi)
         pCL=>this%CL%dataptr(mfi)
         pCG=>this%CG%dataptr(mfi)
         pPLIC=>this%PLIC%dataptr(mfi)
         ! Loop over interior cells
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Handle full cells: set trivial plane
            if (pVF(i,j,k,1).lt.VFlo.or.pVF(i,j,k,1).gt.VFhi) then
               pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,sign(1.0e10_WP,pVF(i,j,k,1)-0.5_WP)]
               cycle
            end if
            ! Liquid-gas symmetry
            flip=.false.; if (pVF(i,j,k,1).ge.0.5_WP) flip=.true.
            ! Initialize geometric moments
            m000=0.0_WP; m100=0.0_WP; m010=0.0_WP; m001=0.0_WP
            ! Construct neighborhood of volume moments (3x3x3 stencil)
            if (flip) then
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+0)=1.0_WP-pVF(ii,jj,kk,1)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(pCG(ii,jj,kk,1)-(this%amr%xlo+(real(ii,WP)+0.5_WP)*dx))*dxi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(pCG(ii,jj,kk,2)-(this%amr%ylo+(real(jj,WP)+0.5_WP)*dy))*dyi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(pCG(ii,jj,kk,3)-(this%amr%zlo+(real(kk,WP)+0.5_WP)*dz))*dzi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(pCL(ii,jj,kk,1)-(this%amr%xlo+(real(ii,WP)+0.5_WP)*dx))*dxi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(pCL(ii,jj,kk,2)-(this%amr%ylo+(real(jj,WP)+0.5_WP)*dy))*dyi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(pCL(ii,jj,kk,3)-(this%amr%zlo+(real(kk,WP)+0.5_WP)*dz))*dzi
                  m000=m000+moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                  m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                  m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                  m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
               end do; end do; end do
            else
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+0)=pVF(ii,jj,kk,1)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(pCL(ii,jj,kk,1)-(this%amr%xlo+(real(ii,WP)+0.5_WP)*dx))*dxi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(pCL(ii,jj,kk,2)-(this%amr%ylo+(real(jj,WP)+0.5_WP)*dy))*dyi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(pCL(ii,jj,kk,3)-(this%amr%zlo+(real(kk,WP)+0.5_WP)*dz))*dzi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(pCG(ii,jj,kk,1)-(this%amr%xlo+(real(ii,WP)+0.5_WP)*dx))*dxi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(pCG(ii,jj,kk,2)-(this%amr%ylo+(real(jj,WP)+0.5_WP)*dy))*dyi
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(pCG(ii,jj,kk,3)-(this%amr%zlo+(real(kk,WP)+0.5_WP)*dz))*dzi
                  m000=m000+moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                  m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                  m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                  m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
               end do; end do; end do
            end if
            ! Geometric center of neighborhood
            center=0.0_WP; if (m000.gt.tiny(1.0_WP)) center=[m100,m010,m001]/m000
            ! Apply symmetry (48 symmetries via reflect_moments)
            call reflect_moments(moments,center,direction,direction2)
            ! Get normal from neural network
            call get_normal(moments,normal)
            normal=normalize(normal)
            ! Undo direction2 rotation (axis permutation)
            if (direction2.eq.1) then
               temp=normal(1); normal(1)=normal(2); normal(2)=temp
            else if (direction2.eq.2) then
               temp=normal(2); normal(2)=normal(3); normal(3)=temp
            else if (direction2.eq.3) then
               temp=normal(1); normal(1)=normal(3); normal(3)=temp
            else if (direction2.eq.4) then
               temp=normal(2); normal(2)=normal(3); normal(3)=temp
               temp=normal(1); normal(1)=normal(2); normal(2)=temp
            else if (direction2.eq.5) then
               temp=normal(1); normal(1)=normal(3); normal(3)=temp
               temp=normal(1); normal(1)=normal(2); normal(2)=temp
            end if
            ! Undo direction reflection (octant)
            if (direction.eq.1) then
               normal(1)=-normal(1)
            else if (direction.eq.2) then
               normal(2)=-normal(2)
            else if (direction.eq.3) then
               normal(3)=-normal(3)
            else if (direction.eq.4) then
               normal(1)=-normal(1); normal(2)=-normal(2)
            else if (direction.eq.5) then
               normal(1)=-normal(1); normal(3)=-normal(3)
            else if (direction.eq.6) then
               normal(2)=-normal(2); normal(3)=-normal(3)
            else if (direction.eq.7) then
               normal(1)=-normal(1); normal(2)=-normal(2); normal(3)=-normal(3)
            end if
            ! Undo liquid-gas flip
            if (.not.flip) normal=-normal
            ! Renormalize
            normal=normalize(normal)
            ! Cell bounds
            lo=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hi=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            ! Store PLIC plane
            pPLIC(i,j,k,:)=[normal(1),normal(2),normal(3),get_plane_dist(normal,lo,hi,pVF(i,j,k,1))]
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      ! Fill ghosts (sync + periodic correction + physical BC)
      call this%fill(lvl,time)
      ! End timer
      this%wt_plicnet=this%wt_plicnet+(MPI_Wtime()-t0)
   end subroutine build_plicnet

   !> Build polygons from PLIC planes
   subroutine build_polygons(this)
      use mpi_f08, only: MPI_Wtime
      use amrvof_geometry, only: cut_hex_polygon
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrvof), intent(inout) :: this
      integer :: lvl
      real(WP) :: dx,dy,dz,t0
      integer :: i,j,k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPLIC
      real(WP), dimension(3) :: lo,hi
      real(WP), dimension(4) :: plane
      real(WP), dimension(3,8) :: hex
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx,gbx
      ! Per-FAB polygon storage (allocatable, indexed by cell)
      real(WP), dimension(:,:,:,:,:), allocatable :: polygon_local  ! (3, 6, ilo:ihi, jlo:jhi, klo:khi)
      integer, dimension(:,:,:), allocatable :: poly_nv_local       ! (ilo:ihi, jlo:jhi, klo:khi)
      real(WP), dimension(3,6) :: poly_verts
      integer :: poly_nv
      ! Start timer
      t0=MPI_Wtime()
      ! Get level and cell size
      lvl=this%amr%maxlvl
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      ! Reset polygon storage
      call this%smesh%reset()
      ! Compute new polygons
      call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
      do while (mfi%next())
         ! Get local and grown boxes
         bx=mfi%tilebox()
         gbx=mfi%growntilebox(2)  ! Grown by 2 for 5x5x5 stencil for curvature calculation
         ! Get pointer to PLIC data
         pPLIC=>this%PLIC%dataptr(mfi)
         ! ----- Step A: Allocate per-FAB polygon storage -----
         allocate(polygon_local(1:3,1:6,gbx%lo(1):gbx%hi(1),gbx%lo(2):gbx%hi(2),gbx%lo(3):gbx%hi(3))); polygon_local=0.0_WP
         allocate(poly_nv_local        (gbx%lo(1):gbx%hi(1),gbx%lo(2):gbx%hi(2),gbx%lo(3):gbx%hi(3))); poly_nv_local=0
         ! ----- Step B: Extract polygons (grown box including ghosts) -----
         do k=gbx%lo(3),gbx%hi(3); do j=gbx%lo(2),gbx%hi(2); do i=gbx%lo(1),gbx%hi(1)
            ! Skip cells with no interface
            if (abs(pPLIC(i,j,k,4)).gt.1.0e+9_WP) cycle         
            ! Build hex and plane
            lo=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hi=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            plane=[pPLIC(i,j,k,1),pPLIC(i,j,k,2),pPLIC(i,j,k,3),pPLIC(i,j,k,4)]
            hex(:,1)=[hi(1),lo(2),lo(3)]
            hex(:,2)=[hi(1),hi(2),lo(3)]
            hex(:,3)=[hi(1),hi(2),hi(3)]
            hex(:,4)=[hi(1),lo(2),hi(3)]
            hex(:,5)=[lo(1),lo(2),lo(3)]
            hex(:,6)=[lo(1),hi(2),lo(3)]
            hex(:,7)=[lo(1),hi(2),hi(3)]
            hex(:,8)=[lo(1),lo(2),hi(3)]
            call cut_hex_polygon(hex,plane,poly_nv,poly_verts)
            ! Store in per-FAB array
            poly_nv_local(i,j,k)=poly_nv
            if (poly_nv.ge.3) polygon_local(:,1:poly_nv,i,j,k)=poly_verts(:,1:poly_nv)
         end do; end do; end do
         ! ----- Step C: Compute curvature (valid cells, stencil access) -----
         ! TODO: curvature = f(polygon_local stencil around i,j,k)
         ! For now, skip curvature computation
         ! ----- Step D: Append to smesh (valid cells only) -----
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            poly_nv=poly_nv_local(i,j,k)
            if (poly_nv.ge.3) call this%smesh%add_polygon(polygon_local(:,1:poly_nv,i,j,k),poly_nv)
         end do; end do; end do
         ! ----- Step E: Deallocate per-FAB storage -----
         deallocate(polygon_local,poly_nv_local)
      end do
      call this%amr%mfiter_destroy(mfi)
      ! End timer
      this%wt_polygon=this%wt_polygon+(MPI_Wtime()-t0)
   end subroutine build_polygons

   !> Reset VF and barycenters from PLIC plane to ensure consistency
   !> Computes in valid + ghost cells from PLIC (which is already filled)
   subroutine reset_moments(this)
      use amrvof_geometry, only: cut_hex_vol
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrvof), intent(inout) :: this
      integer :: lvl,i,j,k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
      real(WP), dimension(3,8) :: hex
      real(WP), dimension(4) :: plane
      real(WP) :: vol_liq,vol_gas,cell_vol,dx,dy,dz
      real(WP), dimension(3) :: bary_liq,bary_gas,cell_center
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      ! Only work at finest level
      lvl=this%amr%maxlvl
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      cell_vol=this%amr%cell_vol(lvl)
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pVF=>this%VF%mf(lvl)%dataptr(mfi)
         pCL=>this%CL%dataptr(mfi)
         pCG=>this%CG%dataptr(mfi)
         pPLIC=>this%PLIC%dataptr(mfi)
         ! Loop over tiles grown by nover
         bx=mfi%growntilebox(this%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Cell center
            cell_center=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
            ! Build hex cell
            hex(:,1)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,2)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,3)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,4)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,5)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            hex(:,6)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            hex(:,7)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            hex(:,8)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            ! Get plane from PLIC
            plane=pPLIC(i,j,k,:)
            ! Skip cutting for full cells
            if (abs(plane(4)).ge.1.0e9_WP) then
               if (plane(4).gt.0.0_WP) then
                  pVF(i,j,k,1)=1.0_WP
               else
                  pVF(i,j,k,1)=0.0_WP
               end if
               pCL(i,j,k,1:3)=cell_center
               pCG(i,j,k,1:3)=cell_center
               cycle
            end if
            ! Cut hex by plane
            call cut_hex_vol(hex,plane,vol_liq,vol_gas,bary_liq,bary_gas)
            ! Update VF and barycenters
            pVF(i,j,k,1)=vol_liq/cell_vol
            pCL(i,j,k,1:3)=bary_liq
            pCG(i,j,k,1:3)=bary_gas
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
   end subroutine reset_moments

   !> Build band around interface using provided VF field
   !> Band values: 1=interface cell, 2..nband=extended layers, 0=far field
   !> Caller owns lifecycle of band MultiFab (build/destroy)
   subroutine build_band(this,lvl,VF,band,nband)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter
      use iso_c_binding, only: c_associated
      use messager, only: die
      implicit none
      class(amrvof), intent(in) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(in) :: VF
      type(amrex_multifab), intent(inout) :: band
      integer, intent(in) :: nband
      integer :: i,j,k,layer
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pVF
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      ! Check that band is allocated, with 1 component and ng=1
      if (.not.c_associated(VF%p)) call die('[build_band] VF is not allocated')
      if (.not.c_associated(band%p)) call die('[build_band] band is not allocated')
      if (band%ncomp() .ne.1) call die('[build_band] band must have 1 component')
      if (band%nghost().lt.1) call die('[build_band] band must have at least 1 ghost cell')
      ! Pass 1: Mark interface cells (band=1)
      call band%setval(0.0_WP)
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pVF  =>VF%dataptr(mfi)
         pBand=>band%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Flag all obvious mixture cells
            if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) then
               pBand(i,j,k,1)=1.0_WP
            ! Flag implicit interfaces
            else if ((pVF(i,j,k,1).lt.VFlo.and.maxval([pVF(i-1,j,k,1),pVF(i+1,j,k,1),pVF(i,j-1,k,1),pVF(i,j+1,k,1),pVF(i,j,k-1,1),pVF(i,j,k+1,1)]).gt.VFhi).or.&
            &        (pVF(i,j,k,1).gt.VFhi.and.minval([pVF(i-1,j,k,1),pVF(i+1,j,k,1),pVF(i,j-1,k,1),pVF(i,j+1,k,1),pVF(i,j,k-1,1),pVF(i,j,k+1,1)]).lt.VFlo)) then
               pBand(i,j,k,1)=1.0_WP
            end if
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      ! Synchronize within level
      call band%fill_boundary(this%amr%geom(lvl))
      ! Pass 2: Extend by nband-1 additional layers
      do layer=2,nband
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pBand=>band%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (pBand(i,j,k,1).eq.0.0_WP.and.any(pBand(i-1:i+1,j-1:j+1,k-1:k+1,1).eq.real(layer-1,WP))) pBand(i,j,k,1)=real(layer,WP)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Synchronize within level
         call band%fill_boundary(this%amr%geom(lvl))
      end do
   end subroutine build_band

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Advect VF using velocity field (staggered or collocated, auto-detected from nodality)
   !> User must provide MultiFabs at finest level with >= 2 ghost cells filled
   subroutine advance_vof(this,U,V,W,dt,time)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrvof), intent(inout) :: this
      type(amrex_multifab), intent(in) :: U,V,W
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      type(amrex_multifab) :: band,Vx,Vy,Vz
      logical :: is_staggered
      integer :: lvl
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,vol
      real(WP) :: t0,t1
      ! Start full routine timer
      t0=MPI_Wtime()

      ! Level at which we're working
      lvl=this%amr%maxlvl

      ! Mesh info
      dx=this%amr%dx(lvl); dxi=1.0_WP/dx
      dy=this%amr%dy(lvl); dyi=1.0_WP/dy
      dz=this%amr%dz(lvl); dzi=1.0_WP/dz
      vol=this%amr%cell_vol(lvl)

      ! Check velocity centering and ghost cell requirements
      check_velocity: block
         use messager, only: die
         logical, dimension(3) :: nodal_U,nodal_V,nodal_W
         nodal_U=U%nodal_type()
         nodal_V=V%nodal_type()
         nodal_W=W%nodal_type()
         if (all(nodal_U .eqv. [.true. ,.false.,.false.]) .and. & 
         &   all(nodal_V .eqv. [.false.,.true. ,.false.]) .and. & 
         &   all(nodal_W .eqv. [.false.,.false.,.true. ])) then
            is_staggered=.true.
         else if (.not.any(nodal_U).and..not.any(nodal_V).and..not.any(nodal_W)) then
            is_staggered=.false.
         else
            call die('[advance_vof] velocity must be either staggered (face-centered) or collocated (cell-centered)')
         end if
         if (U%nghost().lt.2.or.V%nghost().lt.2.or.W%nghost().lt.2) then
            call die('[advance_vof] velocity requires >= 2 ghost cells')
         end if
      end block check_velocity
      
      ! Build transport band to localize computation
      call this%amr%mfab_build(lvl=lvl,mfab=band,ncomp=1,nover=1)
      call this%build_band(lvl=lvl,VF=this%VFold%mf(lvl),band=band,nband=2)

      ! Build face-centered volume flux MultiFabs (8 components: Lvol, Gvol, Lbar(3), Gbar(3))
      call this%amr%mfab_build(lvl=lvl,mfab=Vx,ncomp=8,nover=0,atface=[.true. ,.false.,.false.]); call Vx%setval(0.0_WP)
      call this%amr%mfab_build(lvl=lvl,mfab=Vy,ncomp=8,nover=0,atface=[.false.,.true. ,.false.]); call Vy%setval(0.0_WP)
      call this%amr%mfab_build(lvl=lvl,mfab=Vz,ncomp=8,nover=0,atface=[.false.,.false.,.true. ]); call Vz%setval(0.0_WP)

      ! Phase 1: Compute all fluxes
      t1=MPI_Wtime() ! Start SL timer
      compute_fluxes: block
         use amrvof_geometry, only: tet_sign,tet_map,tet_vol,correct_flux_poly,flux_poly_moments,remap_box_staggered,remap_box_collocated,tet2flux_recursive
         integer :: i,j,k,n,nn
         real(WP), dimension(3,9) :: face
         real(WP), dimension(3,4) :: tet
         integer , dimension(3,4) :: ijk
         integer , dimension(3,9) :: fijk
         type(amrex_multifab) :: proj
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pProj
         integer, dimension(3) :: bblo,bbhi
         logical :: bb_pure_liq,bb_pure_gas
         real(WP) :: fvol,vel
         real(WP), dimension(3) :: fbary
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pVx,pVy,pVz,pU,pV,pW,pPLICold
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx
         ! Build node-centered MultiFab for projected positions (3 components)
         call this%amr%mfab_build(lvl=lvl,mfab=proj,ncomp=3,nover=0,atface=[.true.,.true.,.true.])
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get data pointers: PLICold, band, velocity, fluxes
            pPLICold=>this%PLICold%dataptr(mfi)
            pBand=>band%dataptr(mfi)
            pU =>U%dataptr(mfi)
            pV =>V%dataptr(mfi)
            pW =>W%dataptr(mfi)
            pVx=>Vx%dataptr(mfi)
            pVy=>Vy%dataptr(mfi)
            pVz=>Vz%dataptr(mfi)
            ! Remap all tile nodes via RK2
            pProj=>proj%dataptr(mfi)
            if (is_staggered) then; call remap_box_staggered (bx=mfi%nodaltilebox(),dt=dt,xlo=this%amr%xlo,ylo=this%amr%ylo,zlo=this%amr%zlo,dx=dx,dy=dy,dz=dz,U=pU,Ulo=lbound(pU),V=pV,Vlo=lbound(pV),W=pW,Wlo=lbound(pW),proj=pProj,Projlo=lbound(pProj))
            else;                   call remap_box_collocated(bx=mfi%nodaltilebox(),dt=dt,xlo=this%amr%xlo,ylo=this%amr%ylo,zlo=this%amr%zlo,dx=dx,dy=dy,dz=dz,U=pU,Ulo=lbound(pU),V=pV,Vlo=lbound(pV),W=pW,Wlo=lbound(pW),proj=pProj,Projlo=lbound(pProj))
            end if
            ! X-fluxes: loop over nodaltilebox(1)
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i-1:i,j,k,1)).eq.0.0_WP) cycle
               ! Build face velocity
               vel=merge(pU(i,j,k,1),0.5_WP*(pU(i-1,j,k,1)+pU(i,j,k,1)),is_staggered)
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,5)=pProj(i,j  ,k  ,1:3)
               face(:,2)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=pProj(i,j  ,k+1,1:3)
               face(:,3)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,7)=pProj(i,j+1,k+1,1:3)
               face(:,4)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=pProj(i,j+1,k  ,1:3)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dy*dz*vel)
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(1,nn)=merge(i-1,i,vel.gt.0.0_WP); end do
               ! Check polyhedron bounding box
               bblo=[minval(fijk(1,:)),minval(fijk(2,:)),minval(fijk(3,:))]
               bbhi=[maxval(fijk(1,:)),maxval(fijk(2,:)),maxval(fijk(3,:))]
               bb_pure_liq=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).gt.+1.0e9_WP)) bb_pure_liq=.true.
               bb_pure_gas=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).lt.-1.0e9_WP)) bb_pure_gas=.true.
               ! Compute volume flux
               pVx(i,j,k,:)=0.0_WP
               if (bb_pure_liq.or.bb_pure_gas) then
                  ! Lightweight path
                  call flux_poly_moments(face,fvol,fbary)
                  if (bb_pure_liq) then
                     pVx(i,j,k,1)=fvol; pVx(i,j,k,3:5)=fbary
                  else
                     pVx(i,j,k,2)=fvol; pVx(i,j,k,6:8)=fbary
                  end if
               else
                  ! Decompose into tets, cut, and accumulate
                  do n=1,8
                     do nn=1,4; tet(:,nn)=face(:,tet_map(nn,n)); ijk(:,nn)=fijk(:,tet_map(nn,n)); end do
                     ! Per-tet purity check
                     bblo=[minval(ijk(1,:)),minval(ijk(2,:)),minval(ijk(3,:))]
                     bbhi=[maxval(ijk(1,:)),maxval(ijk(2,:)),maxval(ijk(3,:))]
                     bb_pure_liq=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).gt.+1.0e9_WP)) bb_pure_liq=.true.
                     bb_pure_gas=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).lt.-1.0e9_WP)) bb_pure_gas=.true.
                     if (bb_pure_liq.or.bb_pure_gas) then
                        fvol=tet_vol(tet); fbary=0.25_WP*(tet(:,1)+tet(:,2)+tet(:,3)+tet(:,4))
                        if (bb_pure_liq) then
                           pVx(i,j,k,1)=pVx(i,j,k,1)+fvol; pVx(i,j,k,3:5)=pVx(i,j,k,3:5)+fvol*fbary
                        else
                           pVx(i,j,k,2)=pVx(i,j,k,2)+fvol; pVx(i,j,k,6:8)=pVx(i,j,k,6:8)+fvol*fbary
                        end if
                     else
                        pVx(i,j,k,1:8)=pVx(i,j,k,1:8)+tet_sign(tet)*tet2flux_recursive(tet,ijk,pPLICold,lbound(pPLICold),this%amr%xlo,this%amr%ylo,this%amr%zlo,dx,dy,dz)
                     end if
                  end do
               end if
            end do; end do; end do
            ! Y-fluxes: loop over nodaltilebox(2)
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j-1:j,k,1)).eq.0.0_WP) cycle
               ! Build face velocity
               vel=merge(pV(i,j,k,1),0.5_WP*(pV(i,j-1,k,1)+pV(i,j,k,1)),is_staggered)
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,5)=pProj(i+1,j,k+1,1:3)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=pProj(i  ,j,k+1,1:3)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,7)=pProj(i  ,j,k  ,1:3)
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=pProj(i+1,j,k  ,1:3)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dz*dx*vel)
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(2,nn)=merge(j-1,j,vel.gt.0.0_WP); end do
               ! Check polyhedron bounding box
               bblo=[minval(fijk(1,:)),minval(fijk(2,:)),minval(fijk(3,:))]
               bbhi=[maxval(fijk(1,:)),maxval(fijk(2,:)),maxval(fijk(3,:))]
               bb_pure_liq=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).gt.+1.0e9_WP)) bb_pure_liq=.true.
               bb_pure_gas=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).lt.-1.0e9_WP)) bb_pure_gas=.true.
               ! Compute volume flux
               pVy(i,j,k,:)=0.0_WP
               if (bb_pure_liq.or.bb_pure_gas) then
                  ! Lightweight path
                  call flux_poly_moments(face,fvol,fbary)
                  if (bb_pure_liq) then
                     pVy(i,j,k,1)=fvol; pVy(i,j,k,3:5)=fbary
                  else
                     pVy(i,j,k,2)=fvol; pVy(i,j,k,6:8)=fbary
                  end if
               else
                  ! Decompose into tets, cut, and accumulate
                  do n=1,8
                     do nn=1,4; tet(:,nn)=face(:,tet_map(nn,n)); ijk(:,nn)=fijk(:,tet_map(nn,n)); end do
                     ! Per-tet purity check
                     bblo=[minval(ijk(1,:)),minval(ijk(2,:)),minval(ijk(3,:))]
                     bbhi=[maxval(ijk(1,:)),maxval(ijk(2,:)),maxval(ijk(3,:))]
                     bb_pure_liq=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).gt.+1.0e9_WP)) bb_pure_liq=.true.
                     bb_pure_gas=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).lt.-1.0e9_WP)) bb_pure_gas=.true.
                     if (bb_pure_liq.or.bb_pure_gas) then
                        fvol=tet_vol(tet); fbary=0.25_WP*(tet(:,1)+tet(:,2)+tet(:,3)+tet(:,4))
                        if (bb_pure_liq) then
                           pVy(i,j,k,1)=pVy(i,j,k,1)+fvol; pVy(i,j,k,3:5)=pVy(i,j,k,3:5)+fvol*fbary
                        else
                           pVy(i,j,k,2)=pVy(i,j,k,2)+fvol; pVy(i,j,k,6:8)=pVy(i,j,k,6:8)+fvol*fbary
                        end if
                     else
                        pVy(i,j,k,1:8)=pVy(i,j,k,1:8)+tet_sign(tet)*tet2flux_recursive(tet,ijk,pPLICold,lbound(pPLICold),this%amr%xlo,this%amr%ylo,this%amr%zlo,dx,dy,dz)
                     end if
                  end do
               end if
            end do; end do; end do
            ! Z-fluxes: loop over nodaltilebox(3)
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j,k-1:k,1)).eq.0.0_WP) cycle
               ! Build face velocity
               vel=merge(pW(i,j,k,1),0.5_WP*(pW(i,j,k-1,1)+pW(i,j,k,1)),is_staggered)
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,5)=pProj(i+1,j  ,k,1:3)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,6)=pProj(i  ,j  ,k,1:3)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,7)=pProj(i  ,j+1,k,1:3)
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,8)=pProj(i+1,j+1,k,1:3)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dx*dy*vel)
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(3,nn)=merge(k-1,k,vel.gt.0.0_WP); end do
               ! Check polyhedron bounding box
               bblo=[minval(fijk(1,:)),minval(fijk(2,:)),minval(fijk(3,:))]
               bbhi=[maxval(fijk(1,:)),maxval(fijk(2,:)),maxval(fijk(3,:))]
               bb_pure_liq=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).gt.+1.0e9_WP)) bb_pure_liq=.true.
               bb_pure_gas=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).lt.-1.0e9_WP)) bb_pure_gas=.true.
               ! Compute volume flux
               pVz(i,j,k,:)=0.0_WP
               if (bb_pure_liq.or.bb_pure_gas) then
                  ! Lightweight path
                  call flux_poly_moments(face,fvol,fbary)
                  if (bb_pure_liq) then
                     pVz(i,j,k,1)=fvol; pVz(i,j,k,3:5)=fbary
                  else
                     pVz(i,j,k,2)=fvol; pVz(i,j,k,6:8)=fbary
                  end if
               else
                  ! Decompose into tets, cut, and accumulate
                  do n=1,8
                     do nn=1,4; tet(:,nn)=face(:,tet_map(nn,n)); ijk(:,nn)=fijk(:,tet_map(nn,n)); end do
                     ! Per-tet purity check
                     bblo=[minval(ijk(1,:)),minval(ijk(2,:)),minval(ijk(3,:))]
                     bbhi=[maxval(ijk(1,:)),maxval(ijk(2,:)),maxval(ijk(3,:))]
                     bb_pure_liq=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).gt.+1.0e9_WP)) bb_pure_liq=.true.
                     bb_pure_gas=.false.; if (all(pPLICold(bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3),4).lt.-1.0e9_WP)) bb_pure_gas=.true.
                     if (bb_pure_liq.or.bb_pure_gas) then
                        fvol=tet_vol(tet); fbary=0.25_WP*(tet(:,1)+tet(:,2)+tet(:,3)+tet(:,4))
                        if (bb_pure_liq) then
                           pVz(i,j,k,1)=pVz(i,j,k,1)+fvol; pVz(i,j,k,3:5)=pVz(i,j,k,3:5)+fvol*fbary
                        else
                           pVz(i,j,k,2)=pVz(i,j,k,2)+fvol; pVz(i,j,k,6:8)=pVz(i,j,k,6:8)+fvol*fbary
                        end if
                     else
                        pVz(i,j,k,1:8)=pVz(i,j,k,1:8)+tet_sign(tet)*tet2flux_recursive(tet,ijk,pPLICold,lbound(pPLICold),this%amr%xlo,this%amr%ylo,this%amr%zlo,dx,dy,dz)
                     end if
                  end do
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         call this%amr%mfab_destroy(proj)
      end block compute_fluxes

      ! Phase 2: Update VF from fluxes
      update_vf: block
         use amrvof_geometry, only: project_rk2_staggered,project_rk2_collocated
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         integer :: i,j,k
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pVx,pVy,pVz,pU,pV,pW
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pVFold,pCLold,pCGold
         real(WP) :: Lvol_old,Lvol_new,Lvol_flux
         real(WP) :: Gvol_old,Gvol_new,Gvol_flux
         real(WP), dimension(3) :: Lbar_old,Lbar_new,Lbar_flux
         real(WP), dimension(3) :: Gbar_old,Gbar_new,Gbar_flux
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get data pointers
            pBand=>band%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pCL=>this%CL%dataptr(mfi)
            pCG=>this%CG%dataptr(mfi)
            pVFold=>this%VFold%mf(lvl)%dataptr(mfi)
            pCLold=>this%CLold%dataptr(mfi)
            pCGold=>this%CGold%dataptr(mfi)
            pU=>U%dataptr(mfi)
            pV=>V%dataptr(mfi)
            pW=>W%dataptr(mfi)
            pVx=>Vx%dataptr(mfi)
            pVy=>Vy%dataptr(mfi)
            pVz=>Vz%dataptr(mfi)
            ! Loop over tilebox
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Skip if cell not in band
               if (pBand(i,j,k,1).eq.0.0_WP) cycle
               ! Old phasic moments
               Lvol_old=(       pVFold(i,j,k,1))*vol
               Gvol_old=(1.0_WP-pVFold(i,j,k,1))*vol
               Lbar_old=pCLold(i,j,k,1:3)
               Gbar_old=pCGold(i,j,k,1:3)
               ! Net flux (outflow positive)
               Lvol_flux=pVx(i+1,j,k,1)  -pVx(i,j,k,1)  +pVy(i,j+1,k,1)  -pVy(i,j,k,1)  +pVz(i,j,k+1,1)  -pVz(i,j,k,1)
               Gvol_flux=pVx(i+1,j,k,2)  -pVx(i,j,k,2)  +pVy(i,j+1,k,2)  -pVy(i,j,k,2)  +pVz(i,j,k+1,2)  -pVz(i,j,k,2)
               Lbar_flux=pVx(i+1,j,k,3:5)-pVx(i,j,k,3:5)+pVy(i,j+1,k,3:5)-pVy(i,j,k,3:5)+pVz(i,j,k+1,3:5)-pVz(i,j,k,3:5)
               Gbar_flux=pVx(i+1,j,k,6:8)-pVx(i,j,k,6:8)+pVy(i,j+1,k,6:8)-pVy(i,j,k,6:8)+pVz(i,j,k+1,6:8)-pVz(i,j,k,6:8)
               ! New phasic volumes
               Lvol_new=Lvol_old-Lvol_flux
               Gvol_new=Gvol_old-Gvol_flux
               ! New volume fraction and default barycenters
               pVF(i,j,k,1)=Lvol_new/(Lvol_new+Gvol_new)
               pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               ! Clip and update barycenters
               if (pVF(i,j,k,1).lt.VFlo) then
                  pVF(i,j,k,1)=0.0_WP
               else if (pVF(i,j,k,1).gt.VFhi) then
                  pVF(i,j,k,1)=1.0_WP
               else
                   ! Update barycenters from moment conservation and project forward
                   if (Lvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then
                      Lbar_new=(Lbar_old*Lvol_old-Lbar_flux)/Lvol_new
                      if (is_staggered) then; pCL(i,j,k,1:3)=project_rk2_staggered (Lbar_new,dt,pU,lbound(pU),pV,lbound(pV),pW,lbound(pW),this%amr%xlo,this%amr%ylo,this%amr%zlo,dxi,dyi,dzi)
                      else;                   pCL(i,j,k,1:3)=project_rk2_collocated(Lbar_new,dt,pU,lbound(pU),pV,lbound(pV),pW,lbound(pW),this%amr%xlo,this%amr%ylo,this%amr%zlo,dxi,dyi,dzi)
                      end if
                   end if
                   if (Gvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then
                      Gbar_new=(Gbar_old*Gvol_old-Gbar_flux)/Gvol_new
                      if (is_staggered) then; pCG(i,j,k,1:3)=project_rk2_staggered (Gbar_new,dt,pU,lbound(pU),pV,lbound(pV),pW,lbound(pW),this%amr%xlo,this%amr%ylo,this%amr%zlo,dxi,dyi,dzi)
                      else;                   pCG(i,j,k,1:3)=project_rk2_collocated(Gbar_new,dt,pU,lbound(pU),pV,lbound(pV),pW,lbound(pW),this%amr%xlo,this%amr%ylo,this%amr%zlo,dxi,dyi,dzi)
                      end if
                   end if
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end block update_vf
      this%wt_sl=this%wt_sl+(MPI_Wtime()-t1) ! End SL timer

      ! Cleanup flux multifabs
      call this%amr%mfab_destroy(Vx)
      call this%amr%mfab_destroy(Vy)
      call this%amr%mfab_destroy(Vz)

      ! Clean up band multifab
      call this%amr%mfab_destroy(band)
      
      ! Sync and apply BC
      call this%fill(lvl,time)

      ! End full routine timer
      this%wt_advance=this%wt_advance+(MPI_Wtime()-t0)
   end subroutine advance_vof

   !> Compute advective CFL at finest level
   !> Takes external staggered velocity MultiFabs and returns max CFL
   subroutine get_vof_cfl(this,U,V,W,dt,cfl)
      implicit none
      class(amrvof), intent(inout) :: this
      type(amrex_multifab), intent(in) :: U,V,W  !< Staggered velocity at max level
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      real(WP) :: Umax,Vmax,Wmax,CFLx,CFLy,CFLz
      integer :: lvl
      ! Get finest level metrics
      lvl=this%amr%maxlvl
      ! Get max velocity norms
      Umax=U%norm0()
      Vmax=V%norm0()
      Wmax=W%norm0()
      ! Compute directional CFLs
      CFLx=0.0_WP; CFLy=0.0_WP; CFLz=0.0_WP
      if (this%amr%nx.gt.1) CFLx=dt*Umax/this%amr%dx(lvl)
      if (this%amr%ny.gt.1) CFLy=dt*Vmax/this%amr%dy(lvl)
      if (this%amr%nz.gt.1) CFLz=dt*Wmax/this%amr%dz(lvl)
      ! Return max CFL
      cfl=max(CFLx,CFLy,CFLz)
   end subroutine get_vof_cfl

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Get solver information
   subroutine get_info(this)
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX,MPI_MIN,MPI_INTEGER
      implicit none
      class(amrvof), intent(inout) :: this
      integer :: ierr
      
      ! Get min/max at finest level and integral at level 0
      this%VFmin=this%VF%get_min(lvl=this%amr%maxlvl)
      this%VFmax=this%VF%get_max(lvl=this%amr%maxlvl)
      this%VFint=this%VF%get_sum(lvl=0)*this%amr%cell_vol(0)
      
      ! Load distribution diagnostics at finest level
      load_distribution: block
         use amrex_amr_module, only: amrex_mfiter
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         integer :: nmixed
         integer :: i,j,k
         nmixed=0
         call this%amr%mfiter_build(this%amr%maxlvl,mfi)
         do while (mfi%next())
            pVF=>this%VF%mf(this%amr%maxlvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) nmixed=nmixed+1
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         this%nmixed_max=nmixed; call MPI_ALLREDUCE(MPI_IN_PLACE,this%nmixed_max,1,MPI_INTEGER,MPI_MAX,this%amr%comm,ierr)
         this%nmixed_min=nmixed; call MPI_ALLREDUCE(MPI_IN_PLACE,this%nmixed_min,1,MPI_INTEGER,MPI_MIN,this%amr%comm,ierr)
      end block load_distribution

      ! Reduce per-rank timing to min/max across ranks
      call MPI_ALLREDUCE(this%wt_advance, this%wtmax_advance,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_advance, this%wtmin_advance,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_sl,      this%wtmax_sl,     1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_sl,      this%wtmin_sl,     1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_plic,    this%wtmax_plic,   1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_plic,    this%wtmin_plic,   1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_plicnet, this%wtmax_plicnet,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_plicnet, this%wtmin_plicnet,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_polygon, this%wtmax_polygon,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_polygon, this%wtmin_polygon,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      ! Reset per-rank timing accumulators for next interval
      this%wt_advance=0.0_WP; this%wt_sl=0.0_WP; this%wt_plic=0.0_WP; this%wt_plicnet=0.0_WP; this%wt_polygon=0.0_WP

   end subroutine get_info

   !> Print solver info to screen
   subroutine amrvof_print(this)
      use messager, only: log
      use string, only: str_long
      implicit none
      class(amrvof), intent(in) :: this
      character(len=str_long) :: message
      call log("VOF solver: "//trim(this%name))
      write(message,'("  VF range: [",ES12.5,", ",ES12.5,"]")') VFlo,VFhi
      call log(trim(message))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrvof_print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Register checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrvof), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%VF,'VF')
      call io%add_mfab(this%CL,'CL',this%amr%maxlvl)
      call io%add_mfab(this%CG,'CG',this%amr%maxlvl)
      call io%add_mfab(this%PLIC,'PLIC',this%amr%maxlvl)
   end subroutine register_checkpoint

   !> Restore checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrvof), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      call io%read_data(dirname,this%VF,'VF')
      call io%read_mfab(dirname,this%CL,'CL',this%amr%maxlvl)
      call io%read_mfab(dirname,this%CG,'CG',this%amr%maxlvl)
      call io%read_mfab(dirname,this%PLIC,'PLIC',this%amr%maxlvl)
      ! Fill ghosts
      call this%fill(time=time)
      ! Rebuild polygons from restored PLIC
      call this%build_polygons()
   end subroutine restore_checkpoint

end module amrvof_class
