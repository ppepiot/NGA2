!> AMR VOF Solver class
!> Provides Volume-of-Fluid advection for two-phase flow with AMReX
!> IRL-free implementation using native cutting geometry
module amrvof_class
   use iso_c_binding,    only: c_ptr, c_null_ptr, c_loc, c_f_pointer
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use amrvof_geometry,  only: VFlo, VFhi
   use amrex_amr_module, only: amrex_multifab, amrex_mfiter, amrex_box, &
   &                           amrex_boxarray, amrex_distromap, amrex_geometry
   implicit none
   private

   ! Expose type and dispatchers
   public :: amrvof, VFlo, VFhi
   public :: amrvof_on_init, amrvof_on_coarse, amrvof_on_remake
   public :: amrvof_on_clear, amrvof_tagging, amrvof_postregrid
   
   ! PLIC boundary condition types
   integer, parameter, public :: BC_LIQ     = 1  !< All liquid in ghost
   integer, parameter, public :: BC_GAS     = 2  !< All gas in ghost
   integer, parameter, public :: BC_REFLECT = 3  !< Symmetry (mirror across boundary)
   integer, parameter, public :: BC_USER    = 4  !< User-defined callback

   !> AMR VOF solver type
   type, extends(amrsolver) :: amrvof
      ! User-configurable callbacks
      procedure(vof_init_iface), pointer, nopass :: user_init => null()
      procedure(vof_tagging_iface), pointer, nopass :: user_tagging => null()
      procedure(plic_bc_iface), pointer, nopass :: user_plic_bc => null()

      ! PLIC boundary conditions (per face, only used if direction is non-periodic)
      integer :: bc_xlo = BC_REFLECT
      integer :: bc_xhi = BC_REFLECT
      integer :: bc_ylo = BC_REFLECT
      integer :: bc_yhi = BC_REFLECT
      integer :: bc_zlo = BC_REFLECT
      integer :: bc_zhi = BC_REFLECT

      ! VOF data (solver owns these - 4 MultiFabs as per plan)
      type(amrdata) :: VF           !< Volume fraction (cell-centered)
      type(amrdata) :: Cliq         !< Liquid barycenter (3 components)
      type(amrdata) :: Cgas         !< Gas barycenter (3 components)
      type(amrdata) :: PLIC         !< PLIC plane (4 components: nx, ny, nz, d)

      ! Old data for time stepping
      type(amrdata) :: VFold
      type(amrdata) :: Cliqold
      type(amrdata) :: Cgasold
      type(amrdata) :: PLICold

      ! Monitoring quantities
      real(WP) :: VFmin = 0.0_WP    !< Minimum VF
      real(WP) :: VFmax = 0.0_WP    !< Maximum VF
      real(WP) :: VFint = 0.0_WP    !< Integral of VF (liquid volume)

   contains
      procedure :: initialize
      procedure :: finalize
      ! Override internal type-bound callbacks from amrsolver
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      ! Deferred from amrsolver base class
      procedure :: get_info
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
      ! VOF-specific procedures
      procedure :: build_plic       !< Reconstruct PLIC from VF and barycenters
      procedure :: advance_vof      !< Advect VF using velocity field
      procedure :: fill_vf_lvl      !< Fill VF ghosts at single level
      procedure :: fill_vf          !< Fill VF ghosts on all levels
      procedure :: average_down_vf  !< Average down VF for C/F consistency
      procedure :: reset_moments    !< Recompute VF/barycenters from PLIC
      procedure :: print => amrvof_print  !< Print solver info
   end type amrvof

   !> Abstract interface for user-overridable on_init callback
   abstract interface
      subroutine vof_init_iface(solver, lvl, time, ba, dm)
         import :: amrvof, WP, amrex_boxarray, amrex_distromap
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine vof_init_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine vof_tagging_iface(solver, lvl, tags, time)
         import :: amrvof, c_ptr, WP
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl
         type(c_ptr), intent(in) :: tags
         real(WP), intent(in) :: time
      end subroutine vof_tagging_iface
   end interface

   !> Abstract interface for user-defined PLIC boundary condition
   abstract interface
      subroutine plic_bc_iface(solver, lvl, dir, side, pPLIC, ilo, ihi, jlo, jhi, klo, khi)
         import :: amrvof, WP
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl       !< Level
         integer, intent(in) :: dir       !< Direction (1=x, 2=y, 3=z)
         integer, intent(in) :: side      !< Side (-1=lo, +1=hi)
         real(WP), dimension(:,:,:,:), contiguous, intent(inout) :: pPLIC
         integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi  !< Ghost region bounds
      end subroutine plic_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrvof type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrvof_on_init(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_init(lvl, time, ba, dm)
      if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
   end subroutine amrvof_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrvof_on_coarse(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_coarse(lvl, time, ba, dm)
   end subroutine amrvof_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrvof_on_remake(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_remake(lvl, time, ba, dm)
   end subroutine amrvof_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrvof_on_clear(ctx, lvl)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_clear(lvl)
   end subroutine amrvof_on_clear

   !> Dispatch tagging: force max refinement at interface
   subroutine amrvof_tagging(ctx, lvl, tags, time)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      ! Default: tag interfacial cells (0 < VF < 1)
      ! TODO: implement interface tagging
      if (associated(this%user_tagging)) call this%user_tagging(this, lvl, tags, time)
   end subroutine amrvof_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrvof_postregrid(ctx, lbase, time)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%post_regrid(lbase, time)
   end subroutine amrvof_postregrid

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the VOF solver
   subroutine initialize(this, amr, name)
      class(amrvof), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Set name
      if (present(name)) then
         this%name = trim(name)
      else
         this%name = 'UNNAMED_VOF'
      end if

      ! Store amrgrid pointer
      this%amr => amr

      ! Initialize VOF data (cell-centered)
      call this%VF%initialize(amr, name='VF', ncomp=1, ng=1)
      call this%Cliq%initialize(amr, name='Cliq', ncomp=3, ng=1)
      call this%Cgas%initialize(amr, name='Cgas', ncomp=3, ng=1)
      call this%PLIC%initialize(amr, name='PLIC', ncomp=4, ng=2)

      ! Initialize old data
      call this%VFold%initialize(amr, name='VFold', ncomp=1, ng=1)
      call this%Cliqold%initialize(amr, name='Cliqold', ncomp=3, ng=1)
      call this%Cgasold%initialize(amr, name='Cgasold', ncomp=3, ng=1)
      call this%PLICold%initialize(amr, name='PLICold', ncomp=4, ng=2)

      ! Set parent pointers for callback context access
      this%VF%parent => this
      this%Cliq%parent => this
      this%Cgas%parent => this
      this%PLIC%parent => this
      this%VFold%parent => this
      this%Cliqold%parent => this
      this%Cgasold%parent => this
      this%PLICold%parent => this

      ! Register all 6 callbacks with amrgrid using concrete dispatchers
      ! NOTE: Standard amrdata interpolation is fine IF the tagger ensures the
      !       interface never escapes the finest level between regrids. Tag with
      !       sufficient buffer: interface cells + CFL*dt*regrid_interval.
      select type (this)
       type is (amrvof)
         call this%amr%add_on_init   (amrvof_on_init,    c_loc(this))
         call this%amr%add_on_coarse (amrvof_on_coarse,  c_loc(this))
         call this%amr%add_on_remake (amrvof_on_remake,  c_loc(this))
         call this%amr%add_on_clear  (amrvof_on_clear,   c_loc(this))
         call this%amr%add_tagging   (amrvof_tagging,    c_loc(this))
         call this%amr%add_postregrid(amrvof_postregrid, c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level layouts
      call this%VF%reset_level(lvl, ba, dm)
      call this%Cliq%reset_level(lvl, ba, dm)
      call this%Cgas%reset_level(lvl, ba, dm)
      call this%PLIC%reset_level(lvl, ba, dm)
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
      ! Set to zero
      call this%VF%setval(val=0.0_WP, lvl=lvl)
      call this%Cliq%setval(val=0.0_WP, lvl=lvl)
      call this%Cgas%setval(val=0.0_WP, lvl=lvl)
      call this%PLIC%setval(val=0.0_WP, lvl=lvl)
      call this%VFold%setval(val=0.0_WP, lvl=lvl)
      call this%Cliqold%setval(val=0.0_WP, lvl=lvl)
      call this%Cgasold%setval(val=0.0_WP, lvl=lvl)
      call this%PLICold%setval(val=0.0_WP, lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse
   subroutine on_coarse(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! VF: allocate then interpolate from coarse
      call this%VF%on_coarse(this%VF, lvl, time, ba, dm)
      call this%Cliq%on_coarse(this%Cliq, lvl, time, ba, dm)
      call this%Cgas%on_coarse(this%Cgas, lvl, time, ba, dm)
      call this%PLIC%on_coarse(this%PLIC, lvl, time, ba, dm)
      ! Old data just needs geometry
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid
   subroutine on_remake(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Remake all data
      call this%VF%on_remake(this%VF, lvl, time, ba, dm)
      call this%Cliq%on_remake(this%Cliq, lvl, time, ba, dm)
      call this%Cgas%on_remake(this%Cgas, lvl, time, ba, dm)
      call this%PLIC%on_remake(this%PLIC, lvl, time, ba, dm)
      ! Old data just needs geometry
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this, lvl)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%VF%clear_level(lvl)
      call this%Cliq%clear_level(lvl)
      call this%Cgas%clear_level(lvl)
      call this%PLIC%clear_level(lvl)
      call this%VFold%clear_level(lvl)
      call this%Cliqold%clear_level(lvl)
      call this%Cgasold%clear_level(lvl)
      call this%PLICold%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this, lbase, time)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      call this%average_down_vf(lbase)
   end subroutine post_regrid

   !> Average down VF/Cliq/Cgas from finest to lbase, with ghost fill at each level
   subroutine average_down_vf(this, lbase)
      class(amrvof), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl, lb
      lb = 0; if (present(lbase)) lb = lbase
      do lvl = this%amr%clvl()-1, lb, -1
         ! Average valid cells from fine level
         call this%VF%average_downto(lvl)
         call this%Cliq%average_downto(lvl)
         call this%Cgas%average_downto(lvl)
         ! Fill ghost cells at this coarse level
         call this%VF%sync_lvl(lvl)
         call this%Cliq%sync_lvl(lvl)
         call this%Cgas%sync_lvl(lvl)
      end do
   end subroutine average_down_vf

   !> Fill VF ghost cells at a single level
   subroutine fill_vf_lvl(this, lvl, time)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      ! TODO: implement ghost cell filling
      call this%VF%sync_lvl(lvl)
   end subroutine fill_vf_lvl

   !> Fill VF ghost cells on all levels
   subroutine fill_vf(this, time)
      class(amrvof), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%fill_vf_lvl(lvl, time)
      end do
   end subroutine fill_vf

   !> Finalize the VOF solver
   subroutine finalize(this)
      class(amrvof), intent(inout) :: this
      call this%VF%finalize()
      call this%Cliq%finalize()
      call this%Cgas%finalize()
      call this%PLIC%finalize()
      call this%VFold%finalize()
      call this%Cliqold%finalize()
      call this%Cgasold%finalize()
      call this%PLICold%finalize()
      nullify(this%amr)
      nullify(this%user_init)
      nullify(this%user_tagging)
   end subroutine finalize

   ! ============================================================================
   ! VOF-SPECIFIC METHODS (STUBS)
   ! ============================================================================

   !> Build PLIC reconstruction from VF and barycenters using PLICnet
   subroutine build_plic(this)
      use plicnet, only: get_normal, reflect_moments
      use mathtools, only: normalize
      use amrvof_geometry, only: get_plane_dist
      class(amrvof), intent(inout) :: this
      integer :: lvl, i, j, k, ii, jj, kk, direction, direction2
      real(WP), dimension(0:188) :: moments
      real(WP), dimension(3) :: normal, center, lo, hi
      real(WP) :: m000, m100, m010, m001, temp, vf_cell
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      real(WP) :: dx, dy, dz
      logical :: flip
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      
      ! Only build at finest level
      lvl = this%amr%clvl()
      
      ! Get cell size at this level
      dx = this%amr%dx(lvl)
      dy = this%amr%dy(lvl)
      dz = this%amr%dz(lvl)
      
      ! Iterate over boxes at this level
      call mfi%build(this%VF%mf(lvl), tiling=this%amr%default_tiling)
      do while (mfi%next())
         bx = mfi%tilebox()
         
         ! Get pointers (with ghost cells for stencil access)
         pVF   => this%VF%dataptr(mfi)
         pCliq => this%Cliq%dataptr(mfi)
         pCgas => this%Cgas%dataptr(mfi)
         pPLIC => this%PLIC%dataptr(mfi)
         
         ! Loop over cells in this box
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  
                  vf_cell = pVF(i,j,k,1)
                  
                  ! Handle full cells: set trivial plane
                  if (vf_cell.lt.VFlo .or. vf_cell.gt.VFhi) then
                     pPLIC(i,j,k,1) = 0.0_WP  ! nx
                     pPLIC(i,j,k,2) = 0.0_WP  ! ny
                     pPLIC(i,j,k,3) = 0.0_WP  ! nz
                     pPLIC(i,j,k,4) = sign(1.0e10_WP, vf_cell - 0.5_WP)  ! d
                     cycle
                  end if
                  
                  ! Liquid-gas symmetry
                  flip = .false.
                  if (vf_cell.ge.0.5_WP) flip = .true.
                  
                  ! Initialize geometric moments
                  m000 = 0.0_WP; m100 = 0.0_WP; m010 = 0.0_WP; m001 = 0.0_WP
                  
                  ! Construct neighborhood of volume moments (3x3x3 stencil)
                  if (flip) then
                     do kk = k-1, k+1
                        do jj = j-1, j+1
                           do ii = i-1, i+1
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))   = 1.0_WP - pVF(ii,jj,kk,1)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1) = (pCgas(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) / dx
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2) = (pCgas(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) / dy
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3) = (pCgas(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) / dz
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4) = (pCliq(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) / dx
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5) = (pCliq(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) / dy
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6) = (pCliq(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) / dz
                              m000 = m000 + moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              m100 = m100 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              m010 = m010 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              m001 = m001 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                           end do
                        end do
                     end do
                  else
                     do kk = k-1, k+1
                        do jj = j-1, j+1
                           do ii = i-1, i+1
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))   = pVF(ii,jj,kk,1)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1) = (pCliq(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) / dx
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2) = (pCliq(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) / dy
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3) = (pCliq(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) / dz
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4) = (pCgas(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) / dx
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5) = (pCgas(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) / dy
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6) = (pCgas(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) / dz
                              m000 = m000 + moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              m100 = m100 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              m010 = m010 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              m001 = m001 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                           end do
                        end do
                     end do
                  end if
                  
                  ! Geometric center of neighborhood
                  if (m000.gt.tiny(1.0_WP)) then
                     center = [m100, m010, m001] / m000
                  else
                     center = 0.0_WP
                  end if
                  
                  ! Apply symmetry (48 symmetries via reflect_moments)
                  call reflect_moments(moments, center, direction, direction2)
                  
                  ! Get normal from neural network
                  call get_normal(moments, normal)
                  normal = normalize(normal)
                  
                  ! Undo direction2 rotation (axis permutation)
                  if (direction2.eq.1) then
                     temp = normal(1); normal(1) = normal(2); normal(2) = temp
                  else if (direction2.eq.2) then
                     temp = normal(2); normal(2) = normal(3); normal(3) = temp
                  else if (direction2.eq.3) then
                     temp = normal(1); normal(1) = normal(3); normal(3) = temp
                  else if (direction2.eq.4) then
                     temp = normal(2); normal(2) = normal(3); normal(3) = temp
                     temp = normal(1); normal(1) = normal(2); normal(2) = temp
                  else if (direction2.eq.5) then
                     temp = normal(1); normal(1) = normal(3); normal(3) = temp
                     temp = normal(1); normal(1) = normal(2); normal(2) = temp
                  end if
                  
                  ! Undo direction reflection (octant)
                  if (direction.eq.1) then
                     normal(1) = -normal(1)
                  else if (direction.eq.2) then
                     normal(2) = -normal(2)
                  else if (direction.eq.3) then
                     normal(3) = -normal(3)
                  else if (direction.eq.4) then
                     normal(1) = -normal(1); normal(2) = -normal(2)
                  else if (direction.eq.5) then
                     normal(1) = -normal(1); normal(3) = -normal(3)
                  else if (direction.eq.6) then
                     normal(2) = -normal(2); normal(3) = -normal(3)
                  else if (direction.eq.7) then
                     normal(1) = -normal(1); normal(2) = -normal(2); normal(3) = -normal(3)
                  end if
                  
                  ! Undo liquid-gas flip
                  if (.not.flip) normal = -normal
                  
                  ! Renormalize
                  normal = normalize(normal)
                  
                  ! Cell bounds
                  lo = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k,WP)*dz]
                  hi = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  
                  ! Store PLIC plane: (nx, ny, nz, d)
                  pPLIC(i,j,k,1) = normal(1)
                  pPLIC(i,j,k,2) = normal(2)
                  pPLIC(i,j,k,3) = normal(3)
                  pPLIC(i,j,k,4) = get_plane_dist(normal, lo, hi, vf_cell)
                  
               end do
            end do
         end do
         
      end do
      call mfi%destroy()
      
      ! Sync PLIC ghost cells at finest level
      call this%PLIC%sync_lvl(lvl)
      
      ! Apply periodic correction to plane distance
      call correct_periodic_plic()
      
      ! Apply physical boundary conditions
      call apply_plic_bc()
      
   contains
      
      !> Correct plane distance d in ghost cells for periodicity
      !> After sync, ghost cells that came from the other side of the domain
      !> have plane distance d that refers to coordinates on that side.
      !> We correct by d ← d ± n·L where L is domain length.
      subroutine correct_periodic_plic()
         use amrex_amr_module, only: amrex_geometry
         type(amrex_mfiter) :: mfi2
         type(amrex_geometry) :: geom
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pP
         integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
         integer :: dlo(3), dhi(3)
         real(WP) :: xL, yL, zL
         
         ! Get geometry and domain bounds at this level
         geom = this%amr%geom(lvl)
         dlo = geom%domain%lo
         dhi = geom%domain%hi
         
         xL = this%amr%xhi - this%amr%xlo
         yL = this%amr%yhi - this%amr%ylo
         zL = this%amr%zhi - this%amr%zlo
         
         call mfi2%build(this%PLIC%mf(lvl), tiling=.false.)
         do while (mfi2%next())
            pP => this%PLIC%dataptr(mfi2)
            ilo = lbound(pP,1); ihi = ubound(pP,1)
            jlo = lbound(pP,2); jhi = ubound(pP,2)
            klo = lbound(pP,3); khi = ubound(pP,3)
            
            ! X-periodic: correct ghosts beyond domain
            if (this%amr%xper) then
               ! Low side: FAB extends below domain
               if (ilo .lt. dlo(1)) then
                  do kg = klo, khi
                     do jg = jlo, jhi
                        do ig = ilo, dlo(1)-1
                           pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,1) * xL
                        end do
                     end do
                  end do
               end if
               ! High side: FAB extends above domain
               if (ihi .gt. dhi(1)) then
                  do kg = klo, khi
                     do jg = jlo, jhi
                        do ig = dhi(1)+1, ihi
                           pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,1) * xL
                        end do
                     end do
                  end do
               end if
            end if
            
            ! Y-periodic: correct ghosts beyond domain
            if (this%amr%yper) then
               if (jlo .lt. dlo(2)) then
                  do kg = klo, khi
                     do jg = jlo, dlo(2)-1
                        do ig = ilo, ihi
                           pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,2) * yL
                        end do
                     end do
                  end do
               end if
               if (jhi .gt. dhi(2)) then
                  do kg = klo, khi
                     do jg = dhi(2)+1, jhi
                        do ig = ilo, ihi
                           pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,2) * yL
                        end do
                     end do
                  end do
               end if
            end if
            
            ! Z-periodic: correct ghosts beyond domain
            if (this%amr%zper) then
               if (klo .lt. dlo(3)) then
                  do kg = klo, dlo(3)-1
                     do jg = jlo, jhi
                        do ig = ilo, ihi
                           pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,3) * zL
                        end do
                     end do
                  end do
               end if
               if (khi .gt. dhi(3)) then
                  do kg = dhi(3)+1, khi
                     do jg = jlo, jhi
                        do ig = ilo, ihi
                           pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,3) * zL
                        end do
                     end do
                  end do
               end if
            end if
            
         end do
         call mfi2%destroy()
         
      end subroutine correct_periodic_plic
      
      !> Apply physical boundary conditions to PLIC based on bc_type
      subroutine apply_plic_bc()
         use amrex_amr_module, only: amrex_geometry
         type(amrex_mfiter) :: mfi3
         type(amrex_geometry) :: geom
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pP
         integer :: ilo, ihi, jlo, jhi, klo, khi
         integer :: dlo(3), dhi(3)
         
         geom = this%amr%geom(lvl)
         dlo = geom%domain%lo
         dhi = geom%domain%hi
         
         call mfi3%build(this%PLIC%mf(lvl), tiling=.false.)
         do while (mfi3%next())
            pP => this%PLIC%dataptr(mfi3)
            ilo = lbound(pP,1); ihi = ubound(pP,1)
            jlo = lbound(pP,2); jhi = ubound(pP,2)
            klo = lbound(pP,3); khi = ubound(pP,3)
            
            ! X-low boundary
            if (.not.this%amr%xper .and. ilo.lt.dlo(1)) then
               call apply_bc_face(pP, 1, -1, this%bc_xlo, ilo, dlo(1)-1, jlo, jhi, klo, khi, dlo(1), this%amr%xlo)
            end if
            ! X-high boundary
            if (.not.this%amr%xper .and. ihi.gt.dhi(1)) then
               call apply_bc_face(pP, 1, +1, this%bc_xhi, dhi(1)+1, ihi, jlo, jhi, klo, khi, dhi(1), this%amr%xhi)
            end if
            
            ! Y-low boundary
            if (.not.this%amr%yper .and. jlo.lt.dlo(2)) then
               call apply_bc_face(pP, 2, -1, this%bc_ylo, ilo, ihi, jlo, dlo(2)-1, klo, khi, dlo(2), this%amr%ylo)
            end if
            ! Y-high boundary
            if (.not.this%amr%yper .and. jhi.gt.dhi(2)) then
               call apply_bc_face(pP, 2, +1, this%bc_yhi, ilo, ihi, dhi(2)+1, jhi, klo, khi, dhi(2), this%amr%yhi)
            end if
            
            ! Z-low boundary
            if (.not.this%amr%zper .and. klo.lt.dlo(3)) then
               call apply_bc_face(pP, 3, -1, this%bc_zlo, ilo, ihi, jlo, jhi, klo, dlo(3)-1, dlo(3), this%amr%zlo)
            end if
            ! Z-high boundary
            if (.not.this%amr%zper .and. khi.gt.dhi(3)) then
               call apply_bc_face(pP, 3, +1, this%bc_zhi, ilo, ihi, jlo, jhi, dhi(3)+1, khi, dhi(3), this%amr%zhi)
            end if
            
         end do
         call mfi3%destroy()
         
      end subroutine apply_plic_bc
      
      !> Apply BC to a single face region
      !> x_bnd is the physical coordinate of the boundary face
      subroutine apply_bc_face(pP, dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd, x_bnd)
         real(WP), dimension(:,:,:,:), contiguous, intent(inout) :: pP
         integer, intent(in) :: dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd
         real(WP), intent(in) :: x_bnd
         integer :: ig, jg, kg, isrc, jsrc, ksrc
         
         select case (bc_type)
         
          case (BC_LIQ)
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               pP(ig,jg,kg,1:3) = 0.0_WP
               pP(ig,jg,kg,4) = 1.0e10_WP  ! All liquid
            end do; end do; end do
            
          case (BC_GAS)
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               pP(ig,jg,kg,1:3) = 0.0_WP
               pP(ig,jg,kg,4) = -1.0e10_WP  ! All gas
            end do; end do; end do
            
          case (BC_REFLECT)
            ! Mirror the interface across the boundary
            ! For plane n·x = d, reflected across x_bnd in direction dir:
            ! n'_dir = -n_dir, d' = d - 2*n_dir*x_bnd
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               ! Find mirror source cell
               isrc = ig; jsrc = jg; ksrc = kg
               if (dir.eq.1) isrc = 2*bnd - ig - side
               if (dir.eq.2) jsrc = 2*bnd - jg - side
               if (dir.eq.3) ksrc = 2*bnd - kg - side
               ! Copy and reflect
               pP(ig,jg,kg,1:4) = pP(isrc,jsrc,ksrc,1:4)
               pP(ig,jg,kg,dir) = -pP(ig,jg,kg,dir)  ! Flip normal component
               pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - 2.0_WP*pP(isrc,jsrc,ksrc,dir)*x_bnd
            end do; end do; end do
            
          case (BC_USER)
            ! Call user callback
            if (associated(this%user_plic_bc)) then
               call this%user_plic_bc(this, lvl, dir, side, pP, i1, i2, j1, j2, k1, k2)
            end if
            
          case default
            ! Unknown BC type - do nothing
            
         end select
         
      end subroutine apply_bc_face
      
   end subroutine build_plic
   
   !> Reset VF and barycenters from PLIC plane to ensure consistency
   !> Computes in valid + ghost cells from PLIC (which is already filled)
   !> Then averages down to coarse levels
   subroutine reset_moments(this)
      use amrvof_geometry, only: cut_hex_vol, VFlo, VFhi
      class(amrvof), intent(inout) :: this
      integer :: lvl, i, j, k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      real(WP), dimension(3,8) :: hex
      real(WP), dimension(4) :: plane
      real(WP) :: vol_liq, vol_gas, cell_vol, dx, dy, dz
      real(WP), dimension(3) :: bary_liq, bary_gas, cell_center
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: ilo, ihi, jlo, jhi, klo, khi
      
      ! Only work at finest level
      lvl = this%amr%clvl()
      
      dx = this%amr%dx(lvl)
      dy = this%amr%dy(lvl)
      dz = this%amr%dz(lvl)
      cell_vol = dx * dy * dz
      
      ! Use tiling=.false. to iterate over full FABs including ghosts
      call mfi%build(this%VF%mf(lvl), tiling=.false.)
      do while (mfi%next())
         bx = mfi%fabbox()  ! Full box including ghosts
         
         pVF   => this%VF%dataptr(mfi)
         pCliq => this%Cliq%dataptr(mfi)
         pCgas => this%Cgas%dataptr(mfi)
         pPLIC => this%PLIC%dataptr(mfi)
         
         ! Get array bounds (same as bx but from array)
         ilo = lbound(pVF,1); ihi = ubound(pVF,1)
         jlo = lbound(pVF,2); jhi = ubound(pVF,2)
         klo = lbound(pVF,3); khi = ubound(pVF,3)
         
         do k = klo, khi
            do j = jlo, jhi
               do i = ilo, ihi
                  
                  ! Cell center
                  cell_center = [this%amr%xlo + (real(i,WP)+0.5_WP)*dx, &
                  &              this%amr%ylo + (real(j,WP)+0.5_WP)*dy, &
                  &              this%amr%zlo + (real(k,WP)+0.5_WP)*dz]
                  
                  ! Build hex cell (8 vertices, standard ordering)
                  hex(:,1) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
                  hex(:,2) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
                  hex(:,3) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
                  hex(:,4) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
                  hex(:,5) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  hex(:,6) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  hex(:,7) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  hex(:,8) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  
                  ! Get plane from PLIC (already filled including ghosts)
                  plane(1:3) = pPLIC(i,j,k,1:3)
                  plane(4)   = pPLIC(i,j,k,4)
                  
                  ! Cut hex by plane
                  call cut_hex_vol(hex, plane, vol_liq, vol_gas, bary_liq, bary_gas)
                  
                  ! Update VF
                  pVF(i,j,k,1) = vol_liq / cell_vol
                  
                  ! Update barycenters
                  pCliq(i,j,k,1:3) = bary_liq
                  pCgas(i,j,k,1:3) = bary_gas
                  
                  ! Clean up edge cases
                  if (pVF(i,j,k,1).lt.VFlo) then
                     pVF(i,j,k,1) = 0.0_WP
                     pCliq(i,j,k,1:3) = cell_center
                     pCgas(i,j,k,1:3) = cell_center
                  end if
                  if (pVF(i,j,k,1).gt.VFhi) then
                     pVF(i,j,k,1) = 1.0_WP
                     pCliq(i,j,k,1:3) = cell_center
                     pCgas(i,j,k,1:3) = cell_center
                  end if
                  
               end do
            end do
         end do
         
      end do
      call mfi%destroy()
      
      ! Average down to coarse levels (uses ghost-capable average_down)
      call this%average_down_vf()
      
   end subroutine reset_moments

   !> Advect VF using provided velocity field
   subroutine advance_vof(this, U, V, W, dt)
      class(amrvof), intent(inout) :: this
      type(amrdata), intent(in) :: U, V, W
      real(WP), intent(in) :: dt
      ! TODO: implement VOF advection using native geometry
   end subroutine advance_vof

   ! ============================================================================
   ! DEFERRED METHODS
   ! ============================================================================

   !> Get solver information
   subroutine get_info(this)
      use parallel, only: MPI_REAL_WP
      use mpi_f08
      class(amrvof), intent(inout) :: this
      integer :: lvl, ierr

      ! Initialize
      this%VFmin = huge(1.0_WP)
      this%VFmax = -huge(1.0_WP)
      this%VFint = 0.0_WP

      ! Loop over levels
      do lvl = 0, this%amr%clvl()
         this%VFmin = min(this%VFmin, this%VF%min(lvl=lvl))
         this%VFmax = max(this%VFmax, this%VF%norm0(lvl=lvl))
      end do

      ! Compute volume integral at level 0
      this%VFint = this%VF%get_sum(lvl=0) * this%amr%dx(0) * this%amr%dy(0) * this%amr%dz(0)

      ! Reduce across MPI ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%VFint, 1, MPI_REAL_WP, MPI_SUM, this%amr%comm, ierr)
   end subroutine get_info

   !> Register checkpoint
   subroutine register_checkpoint(this, io)
      use amrio_class, only: amrio
      class(amrvof), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_multifab(this%VF, 'VF')
      call io%add_multifab(this%Cliq, 'Cliq')
      call io%add_multifab(this%Cgas, 'Cgas')
      call io%add_multifab(this%PLIC, 'PLIC')
   end subroutine register_checkpoint

   !> Restore checkpoint
   subroutine restore_checkpoint(this, io, dirname)
      use amrio_class, only: amrio
      class(amrvof), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      call io%read_multifab(this%VF, 'VF', dirname)
      call io%read_multifab(this%Cliq, 'Cliq', dirname)
      call io%read_multifab(this%Cgas, 'Cgas', dirname)
      call io%read_multifab(this%PLIC, 'PLIC', dirname)
   end subroutine restore_checkpoint

   !> Print solver info to screen
   subroutine amrvof_print(this)
      use messager, only: log
      use string, only: str_long
      class(amrvof), intent(in) :: this
      character(len=str_long) :: message
      call log("VOF solver: "//trim(this%name))
      write(message,'("  VF range: [",ES12.5,", ",ES12.5,"]")') VFlo, VFhi
      call log(trim(message))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrvof_print

end module amrvof_class
