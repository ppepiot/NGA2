!> Data class wrapping AMReX MultiFab
!> Designed to be managed by amrgrid Registry
module amrdata_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrex_amr_module, only: amrex_multifab,amrex_boxarray,amrex_distromap,&
   &                           amrex_multifab_build,amrex_multifab_destroy,amrex_geometry,&
   &                           amrex_interp_cell_cons,amrex_interp_face_linear,amrex_interp_node_bilinear
   implicit none
   private

   public :: amrdata
   public :: amrdata_on_init,amrdata_on_coarse,amrdata_on_remake,amrdata_on_clear,amrdata_fillbc

   ! Special interpolation modes for amrdata
   integer, parameter, public :: amrex_interp_none   = -1  !< Workspace: allocate but don't fill
   integer, parameter, public :: amrex_interp_reinit = -2  !< Reinit: always call user_init on regrid

   ! Forward declaration for interfaces
   type :: amrdata
      ! The underlying AMReX objects (array of MultiFabs, one per level)
      type(amrex_multifab), dimension(:), allocatable :: mf
      ! Pointer to parent amrgrid (set by initialize)
      class(amrgrid), pointer :: amr => null()
      ! Pointer to parent solver (set by solver, for BC context access)
      class(*), pointer :: parent => null()
      ! Metadata
      character(len=str_medium) :: name='UNNAMED_AMRDATA'
      integer :: ncomp=1                                   !< Number of components (default=1)
      integer :: ng=0                                      !< Number of ghost cells (default=0)
      logical :: nodal(3) = [.false., .false., .false.]    !< false=cell, true=vertex in that direction
      integer :: interp=amrex_interp_cell_cons             !< Interpolation method
      integer, dimension(:,:), allocatable :: lo_bc,hi_bc  !< Boundary conditions: lo_bc(3,ncomp), hi_bc(3,ncomp)
      ! Callback pointers (set to defaults in initialize)
      procedure(on_init_iface),   pointer, nopass :: on_init   => null()
      procedure(on_coarse_iface), pointer, nopass :: on_coarse => null()
      procedure(on_remake_iface), pointer, nopass :: on_remake => null()
      procedure(on_clear_iface),  pointer, nopass :: on_clear  => null()
      procedure(fillbc_iface),    pointer, nopass :: fillbc    => null()
      ! User-provided initialization callback
      procedure(on_init_iface),   pointer, nopass :: user_init => null()
   contains
      ! Lifecycle methods
      procedure :: initialize       !< Initialize amrdata with amrgrid and parameters
      procedure :: finalize         !< Finalize the amrdata object
      procedure :: register         !< Register regrid callbacks with amrgrid
      ! Level management
      procedure :: reset_level      !< Regenerates mfab at level
      procedure :: clear_level      !< Clears mfab at level
      ! Fill operations
      procedure :: fill_from_coarse !< Interpolate from coarse level only
      procedure :: fill             !< Fill ghost cells and coarse-fine boundary data
      procedure :: fill_mfab        !< Fill into a target MultiFab
      procedure :: average_down     !< Average from finest to coarsest level
      procedure :: average_downto   !< Average level lvl+1 down to level lvl
      ! Iteration helper
      procedure :: mfiter_build     !< Build MFIter from this data's MultiFab
   end type amrdata

   !> Abstract interface for on_init callback
   abstract interface
      subroutine on_init_iface(this, lvl, time, ba, dm)
         import :: amrdata, amrex_boxarray, amrex_distromap, WP
         class(amrdata), intent(inout) :: this
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine on_init_iface
   end interface

   !> Abstract interface for on_coarse callback
   abstract interface
      subroutine on_coarse_iface(this, lvl, time, ba, dm)
         import :: amrdata, amrex_boxarray, amrex_distromap, WP
         class(amrdata), intent(inout) :: this
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine on_coarse_iface
   end interface

   !> Abstract interface for on_remake callback
   abstract interface
      subroutine on_remake_iface(this, lvl, time, ba, dm)
         import :: amrdata, amrex_boxarray, amrex_distromap, WP
         class(amrdata), intent(inout) :: this
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine on_remake_iface
   end interface

   !> Abstract interface for on_clear callback
   abstract interface
      subroutine on_clear_iface(this, lvl)
         import :: amrdata
         class(amrdata), intent(inout) :: this
         integer, intent(in) :: lvl
      end subroutine on_clear_iface
   end interface

   !> Abstract interface for fillbc callback (matches AMReX order)
   abstract interface
      subroutine fillbc_iface(this, mf, scomp, ncomp, time, geom)
         import :: amrdata, amrex_multifab, amrex_geometry, WP
         class(amrdata), intent(inout) :: this
         type(amrex_multifab), intent(inout) :: mf
         integer, intent(in) :: scomp, ncomp
         real(WP), intent(in) :: time
         type(amrex_geometry), intent(in) :: geom
      end subroutine fillbc_iface
   end interface

contains

   !> Initialize amrdata with amrgrid reference and parameters
   !> This sets the amr pointer, allocates mf array, and configures BCs
   subroutine initialize(this, amr, name, ncomp, ng, nodal, interp)
      use amrex_amr_module, only: amrex_bc_int_dir
      use amrex_amr_module, only: amrex_interp_face_linear,amrex_interp_node_bilinear
      implicit none
      class(amrdata), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in) :: name
      integer, intent(in) :: ncomp
      integer, intent(in) :: ng
      logical, intent(in), optional :: nodal(3)  !< Data location: false=cell, true=vertex
      integer, intent(in), optional :: interp    !< Interpolation method override
      ! Store pointer to amrgrid
      this%amr => amr
      ! Set metadata
      this%name = name
      this%ncomp = ncomp
      this%ng = ng
      ! Set nodal if provided
      if (present(nodal)) this%nodal = nodal
      ! Set interpolation: user override or auto-select based on nodal
      if (present(interp)) then
         this%interp = interp
      else if (any(this%nodal)) then
         if (all(this%nodal)) then
            this%interp = amrex_interp_node_bilinear  ! Node-centered
         else
            this%interp = amrex_interp_face_linear    ! Face-centered
         end if
      end if
      ! Allocate mf array
      if (.not.allocated(this%mf)) allocate(this%mf(0:amr%maxlvl))
      ! Allocate and set default BCs to handle periodicity only
      if (.not.allocated(this%lo_bc)) allocate(this%lo_bc(3,ncomp))
      if (.not.allocated(this%hi_bc)) allocate(this%hi_bc(3,ncomp))
      this%lo_bc = amrex_bc_int_dir
      this%hi_bc = amrex_bc_int_dir
      ! Set default callbacks
      this%on_init   => default_on_init
      this%on_coarse => default_on_coarse
      this%on_remake => default_on_remake
      this%on_clear  => default_on_clear
      this%fillbc    => default_fillbc
   end subroutine initialize

   !> Finalize the amrdata object
   subroutine finalize(this)
      class(amrdata), intent(inout) :: this
      integer :: i
      ! Destroy all MultiFabs
      if (allocated(this%mf)) then
         do i=0,ubound(this%mf,1)
            call amrex_multifab_destroy(this%mf(i))
         end do
         deallocate(this%mf)
      end if
      ! Deallocate BC arrays
      if (allocated(this%lo_bc)) deallocate(this%lo_bc)
      if (allocated(this%hi_bc)) deallocate(this%hi_bc)
      ! Reset pointers
      nullify(this%amr)
      nullify(this%parent)
      nullify(this%on_init)
      nullify(this%on_coarse)
      nullify(this%on_remake)
      nullify(this%on_clear)
      nullify(this%fillbc)
      nullify(this%user_init)
      ! Reset scalars to defaults
      this%name = 'UNNAMED_AMRDATA'
      this%ncomp = 1
      this%ng = 0
      this%nodal = [.false., .false., .false.]
      this%interp = amrex_interp_cell_cons
   end subroutine finalize

   !> Register regrid callbacks with amrgrid
   subroutine register(this)
      use iso_c_binding, only: c_loc
      use messager, only: die
      class(amrdata), target, intent(inout) :: this
      if (.not.associated(this%amr)) call die('[amrdata%register] armdata%initialize must be called first')
      select type (this)
       type is (amrdata)
         call this%amr%add_on_init  (amrdata_on_init,   c_loc(this))
         call this%amr%add_on_coarse(amrdata_on_coarse, c_loc(this))
         call this%amr%add_on_remake(amrdata_on_remake, c_loc(this))
         call this%amr%add_on_clear (amrdata_on_clear,  c_loc(this))
      end select
   end subroutine register

   !> Reset mfab on a level given new BoxArray and DistroMap
   subroutine reset_level(this,lvl,ba,dm)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray),  intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      call amrex_multifab_destroy(this%mf(lvl))
      call amrex_multifab_build(this%mf(lvl),ba,dm,this%ncomp,this%ng,this%nodal)
   end subroutine reset_level

   !> Destroy mfab on a level
   subroutine clear_level(this,lvl)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      call amrex_multifab_destroy(this%mf(lvl))
   end subroutine clear_level

   !> Build MFIter from this data's MultiFab (ensures correct layout for staggered data)
   subroutine mfiter_build(this, lvl, mfi, tiling)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(inout) :: mfi
      logical, intent(in), optional :: tiling
      logical :: use_tiling
      use_tiling = this%amr%default_tiling
      if (present(tiling)) use_tiling = tiling
      call amrex_mfiter_build(mfi, this%mf(lvl), tiling=use_tiling)
   end subroutine mfiter_build

   !> Default fillbc: applies amrex_filcc using lo_bc/hi_bc
   !> NOTE: amrex_filcc is designed for cell-centered data only.
   subroutine default_fillbc(this, mf, scomp, ncomp, time, geom)
      use amrex_amr_module, only: amrex_filcc, amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp, ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer, dimension(4) :: plo, phi
      ! Skip if all directions periodic
      if (this%amr%xper .and. this%amr%yper .and. this%amr%zper) return
      ! Loop over boxes and apply filcc
      call amrex_mfiter_build(mfi, mf, tiling=.false.)
      do while (mfi%next())
         p => mf%dataptr(mfi)
         if (.not.geom%domain%contains(p)) then
            plo = lbound(p); phi = ubound(p)
            call amrex_filcc(p, plo, phi, geom%domain%lo, geom%domain%hi, geom%dx, &
            &   geom%get_physical_location(plo), this%lo_bc, this%hi_bc)
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine default_fillbc

   !> Default on_init: reset level and initialize based on interp mode
   subroutine default_on_init(this, lvl, time, ba, dm)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level
      call this%reset_level(lvl, ba, dm)
      ! Handle different modes
      select case (this%interp)
       case (amrex_interp_none)
         ! Workspace mode: just allocate, don't fill
       case default
         ! Standard or reinit mode: zero out, then call user_init
         call this%mf(lvl)%setval(0.0_WP)
      end select
      ! User-provided initialization (called for all modes except none without user_init)
      if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
   end subroutine default_on_init

   !> Default on_coarse: reset level and fill based on interp mode
   subroutine default_on_coarse(this, lvl, time, ba, dm)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level
      call this%reset_level(lvl, ba, dm)
      ! Handle different modes
      select case (this%interp)
       case (amrex_interp_none)
         ! Workspace mode: just allocate, don't fill
       case (amrex_interp_reinit)
         ! Reinit mode: call user_init instead of interpolating
         if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
       case default
         ! Standard interpolation: fill from coarse
         call this%fill_from_coarse(lvl, time)
      end select
   end subroutine default_on_coarse

   !> Default on_remake: handle level relayout based on interp mode
   subroutine default_on_remake(this, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_multifab_build, amrex_multifab_destroy
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_multifab) :: mf_tmp
      ! Handle different modes
      select case (this%interp)
       case (amrex_interp_none)
         ! Workspace mode: just reallocate
         call this%reset_level(lvl, ba, dm)
       case (amrex_interp_reinit)
         ! Reinit mode: reallocate and call user_init
         call this%reset_level(lvl, ba, dm)
         if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
       case default
         ! Standard interpolation: FillPatch old data into new layout
         ! Build temp MultiFab with new layout (0 ghost cells for FillPatch)
         call amrex_multifab_build(mf_tmp, ba, dm, this%ncomp, 0, this%nodal)
         ! Fill temp from old data via FillPatch
         call this%fill_mfab(mf_tmp, lvl, time)
         ! Reset level and copy from temp
         call this%reset_level(lvl, ba, dm)
         call this%mf(lvl)%copy(mf_tmp, 1, 1, this%ncomp, 0)
         ! Destroy temp
         call amrex_multifab_destroy(mf_tmp)
      end select
   end subroutine default_on_remake

   !> Default on_clear: just clear the level
   subroutine default_on_clear(this, lvl)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%clear_level(lvl)
   end subroutine default_on_clear

   !> Dispatch callback for fillbc
   subroutine amrdata_fillbc(ctx, mf_ptr, scomp, ncomp, time, geom_ptr) bind(c)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_int
      use amrex_amr_module, only: amrex_geometry, amrex_multifab
      type(c_ptr), value, intent(in) :: ctx
      type(c_ptr), value, intent(in) :: mf_ptr
      integer(c_int), value, intent(in) :: scomp
      integer(c_int), value, intent(in) :: ncomp
      real(c_double), value, intent(in) :: time
      type(c_ptr), value, intent(in) :: geom_ptr
      type(amrdata), pointer :: this
      type(amrex_multifab) :: mf
      type(amrex_geometry) :: geom
      call c_f_pointer(ctx, this)
      mf = mf_ptr
      geom = geom_ptr
      ! Call the fillbc callback
      call this%fillbc(this, mf, int(scomp), int(ncomp), real(time, WP), geom)
   end subroutine amrdata_fillbc

   !> Dispatch callback for on_init
   subroutine amrdata_on_init(ctx, lvl, time, ba, dm)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_init(this, lvl, time, ba, dm)
   end subroutine amrdata_on_init

   !> Dispatch callback for on_coarse
   subroutine amrdata_on_coarse(ctx, lvl, time, ba, dm)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_coarse(this, lvl, time, ba, dm)
   end subroutine amrdata_on_coarse

   !> Dispatch callback for on_remake
   subroutine amrdata_on_remake(ctx, lvl, time, ba, dm)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_remake(this, lvl, time, ba, dm)
   end subroutine amrdata_on_remake

   !> Dispatch callback for on_clear
   subroutine amrdata_on_clear(ctx, lvl)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_clear(this, lvl)
   end subroutine amrdata_on_clear

   !> Fill fine level from coarse only (for creating new fine levels)
   subroutine fill_from_coarse(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillcoarsepatch
      implicit none
      class(amrdata), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: data_ctx
      type(c_funptr) :: bc_dispatch_ptr
      ! Only valid for fine levels
      if (lvl .lt. 1) return
      ! Get context pointer
      select type (this)
       type is (amrdata)
         data_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(amrdata_fillbc)
      ! Call C++ wrapper
      call amrmfab_fillcoarsepatch(this%mf(lvl)%p, time, this%mf(lvl-1)%p, &
      &   this%amr%geom(lvl-1)%p, this%amr%geom(lvl)%p, data_ctx, bc_dispatch_ptr, &
      &   1, 1, this%ncomp, this%amr%rref(lvl-1), this%interp, this%lo_bc, this%hi_bc, this%ncomp)
   end subroutine fill_from_coarse

   !> Fill ghost cells and coarse-fine boundary data
   !> For lvl=0: fills ghost cells using boundary conditions
   !> For lvl>0: interpolates from coarser level and fills ghosts
   subroutine fill(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two
      implicit none
      class(amrdata), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      real(WP) :: t_old, t_new
      type(c_ptr) :: data_ctx
      type(c_funptr) :: bc_dispatch_ptr
      t_old = time - 1.0e200_WP
      t_new = time
      ! Get context pointer (this amrdata)
      select type (this)
       type is (amrdata)
         data_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(amrdata_fillbc)
      ! Call appropriate FillPatch (scomp/dcomp use 1-indexed Fortran convention)
      if (lvl .eq. 0) then
         call amrmfab_fillpatch_single(this%mf(0)%p, t_old, this%mf(0)%p, &
         &   t_new, this%mf(0)%p, this%amr%geom(0)%p, data_ctx, bc_dispatch_ptr, &
         &   time, 1, 1, this%ncomp)
      else
         call amrmfab_fillpatch_two(this%mf(lvl)%p, t_old, this%mf(lvl-1)%p, &
         &   t_new, this%mf(lvl-1)%p, this%amr%geom(lvl-1)%p, &
         &   t_old, this%mf(lvl)%p, t_new, this%mf(lvl)%p, this%amr%geom(lvl)%p, &
         &   data_ctx, bc_dispatch_ptr, time, 1, 1, this%ncomp, &
         &   this%amr%rref(lvl-1), this%interp, this%lo_bc, this%hi_bc, this%ncomp)
      end if
   end subroutine fill

   !> Fill into a target MultiFab from this amrdata (for regrid callbacks)
   subroutine fill_mfab(this, dest, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_amr_module, only: amrex_multifab
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two
      implicit none
      class(amrdata), target, intent(inout) :: this
      type(amrex_multifab), intent(inout) :: dest
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      real(WP) :: t_old, t_new
      type(c_ptr) :: data_ctx
      type(c_funptr) :: bc_dispatch_ptr
      t_old = time - 1.0e200_WP
      t_new = time
      ! Get context pointer
      select type (this)
       type is (amrdata)
         data_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(amrdata_fillbc)
      ! Call appropriate FillPatch
      if (lvl .eq. 0) then
         call amrmfab_fillpatch_single(dest%p, t_old, this%mf(0)%p, &
         &   t_new, this%mf(0)%p, this%amr%geom(0)%p, data_ctx, bc_dispatch_ptr, &
         &   time, 1, 1, this%ncomp)
      else
         call amrmfab_fillpatch_two(dest%p, t_old, this%mf(lvl-1)%p, &
         &   t_new, this%mf(lvl-1)%p, this%amr%geom(lvl-1)%p, &
         &   t_old, this%mf(lvl)%p, t_new, this%mf(lvl)%p, this%amr%geom(lvl)%p, &
         &   data_ctx, bc_dispatch_ptr, time, 1, 1, this%ncomp, &
         &   this%amr%rref(lvl-1), this%interp, this%lo_bc, this%hi_bc, this%ncomp)
      end if
   end subroutine fill_mfab

   !> Average down from finest level to coarsest (ensures level consistency)
   subroutine average_down(this)
      use amrex_amr_module, only: amrex_average_down
      implicit none
      class(amrdata), intent(inout) :: this
      integer :: lvl
      if (.not.associated(this%amr)) return
      ! Loop from finest to coarsest
      do lvl = this%amr%clvl()-1, 0, -1
         call amrex_average_down(this%mf(lvl+1), this%mf(lvl), &
         &   this%amr%geom(lvl+1), this%amr%geom(lvl), 1, this%ncomp, this%amr%rref(lvl))
      end do
   end subroutine average_down

   !> Average level lvl+1 down to level lvl
   subroutine average_downto(this, lvl)
      use amrex_amr_module, only: amrex_average_down
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      if (.not.associated(this%amr)) return
      if (lvl.lt.0.or.lvl.ge.this%amr%clvl()) return
      call amrex_average_down(this%mf(lvl+1), this%mf(lvl), &
      &   this%amr%geom(lvl+1), this%amr%geom(lvl), 1, this%ncomp, this%amr%rref(lvl))
   end subroutine average_downto

end module amrdata_class
