!> Data class wrapping AMReX MultiFab
!> Designed to be managed by amrgrid Registry
module amrdata_class
   use iso_c_binding,    only: c_ptr
   use precision,        only: WP
   use string,           only: str_medium
   use amrex_amr_module, only: amrex_multifab,amrex_boxarray,amrex_distromap,&
   &                           amrex_multifab_build,amrex_multifab_destroy,amrex_geometry,&
   &                           amrex_interp_cell_cons
   implicit none
   private

   public :: amrdata

   !> Data object wrapping a MultiFab hierarchy (one MultiFab per level)
   type :: amrdata
      ! The underlying AMReX objects (array of MultiFabs, one per level)
      type(amrex_multifab), dimension(:), allocatable :: mf
      ! Metadata
      character(len=str_medium) :: name='UNNAMED_AMRDATA'
      integer :: ncomp=1                                   !< Number of components
      integer :: ng=0                                      !< Number of ghost cells
      integer :: interp=amrex_interp_cell_cons             !< Interpolation method
      integer, dimension(:,:), allocatable :: lo_bc,hi_bc  !< Boundary conditions: lo_bc(3,ncomp), hi_bc(3,ncomp)
   contains
      procedure :: define
      procedure :: destroy
      procedure :: clear_level
      procedure :: fill_ghosts
      procedure :: fillbc          !< Default physical BC callback (uses amrex_filcc)
      procedure :: get_data_ptr
      procedure :: on_regrid
   end type amrdata

contains

   !> Define the data on a specific grid level (BoxArray + DistroMap)
   !> ncomp and ng should already be set via amrgrid%register
   subroutine define(this,lvl,ba,dm)
      use iso_c_binding, only: c_associated
      use messager, only: die
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Verify allocation (should be handled by amrgrid%register)
      if (.not.allocated(this%mf)) call die('[amrdata define] mf not allocated - use amrgrid%register')
      if (lvl.lt.0.or.lvl.gt.ubound(this%mf,1)) call die('[amrdata define] lvl out of bounds')
      ! Destroy if already allocated at this level
      if (c_associated(this%mf(lvl)%p)) call amrex_multifab_destroy(this%mf(lvl))
      ! Build the MultiFab
      call amrex_multifab_build(this%mf(lvl),ba,dm,this%ncomp,this%ng)
   end subroutine define

   !> Destroy data at a specific level
   subroutine clear_level(this,lvl)
      use iso_c_binding, only: c_associated
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      if (allocated(this%mf)) then
         if (lvl.ge.0.and.lvl.le.ubound(this%mf,1)) then
            if (c_associated(this%mf(lvl)%p)) call amrex_multifab_destroy(this%mf(lvl))
         end if
      end if
   end subroutine clear_level

   !> Destroy the underlying MultiFabs and BC arrays
   subroutine destroy(this)
      use iso_c_binding, only: c_associated
      class(amrdata), intent(inout) :: this
      integer :: i
      if (allocated(this%mf)) then
         do i=0,ubound(this%mf,1)
            if (c_associated(this%mf(i)%p)) call amrex_multifab_destroy(this%mf(i))
         end do
         deallocate(this%mf)
      end if
      if (allocated(this%lo_bc)) deallocate(this%lo_bc)
      if (allocated(this%hi_bc)) deallocate(this%hi_bc)
   end subroutine destroy

   !> Fill ghost cells (boundary conditions) for a specific level
   subroutine fill_ghosts(this,lvl,geom)
      use iso_c_binding, only: c_associated
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_geometry), intent(in) :: geom
      if (allocated(this%mf).and.lvl.ge.0.and.lvl.le.ubound(this%mf,1)) then
         if (c_associated(this%mf(lvl)%p)) call this%mf(lvl)%fill_boundary(geom)
      end if
   end subroutine fill_ghosts

   !> Get access to raw data pointer for a specific block and level
   subroutine get_data_ptr(this,lvl,mfi,ptr)
      use amrex_amr_module, only: amrex_mfiter
      use iso_c_binding, only: c_double
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      real(c_double), dimension(:,:,:,:), pointer, intent(out) :: ptr
      ptr=>this%mf(lvl)%dataPtr(mfi)
   end subroutine get_data_ptr

   !> Hook called by Registry after regridding
   subroutine on_regrid(this,lvl,time)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      ! Default: do nothing (override for fillpatch, etc.)
   end subroutine on_regrid

   !> Default physical BC callback - applies amrex_filcc using lo_bc/hi_bc
   !> This is called by amrgrid%fill() via the C++ FillPatch wrapper
   subroutine fillbc(this,mf_ptr,geom_ptr,time,scomp,ncomp)
      use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build
      use iso_c_binding, only: c_ptr,c_int,c_f_pointer
      class(amrdata), intent(inout) :: this
      type(c_ptr), value :: mf_ptr
      type(c_ptr), value :: geom_ptr
      real(WP), value :: time
      integer, value :: scomp,ncomp
      type(amrex_geometry) :: geom
      type(amrex_multifab) :: mf
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer, dimension(4) :: plo,phi
      ! Convert C pointers to Fortran types
      geom=geom_ptr
      mf=mf_ptr
      ! Loop over boxes and apply filcc
      call amrex_mfiter_build(mfi,mf)
      do while(mfi%next())
         p=>mf%dataptr(mfi)
         ! Check if part of box is outside the domain
         if (.not.geom%domain%contains(p)) then
            plo=lbound(p); phi=ubound(p)
            call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,&
            &   geom%get_physical_location(plo),this%lo_bc,this%hi_bc)
         end if
      end do
   end subroutine fillbc

end module amrdata_class
