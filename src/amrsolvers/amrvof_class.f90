!> AMR VOF Solver class
!> Provides Volume-of-Fluid advection for two-phase flow with AMReX
!> IRL-free implementation using native cutting geometry
module amrvof_class
   use iso_c_binding,    only: c_ptr, c_null_ptr, c_loc, c_f_pointer, c_char
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use amrvof_geometry,  only: VFlo, VFhi
   use surfmesh_class,   only: surfmesh
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
      procedure(vof_bc_iface), pointer, nopass :: user_vof_bc => null()

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
      
      ! Surface mesh for visualization (finest level polygons)
      type(surfmesh) :: smesh

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
      procedure :: build_plic             !< Reconstruct PLIC from VF and barycenters
      procedure :: advance_vof_stag       !< Advect VF using staggered velocity
      procedure :: advance_vof_col        !< Advect VF using collocated velocity
      procedure :: fill_moments           !< Fill VF/Cliq/Cgas ghosts (sync + BC)
      procedure :: sync_moments_lvl       !< Sync VF/Cliq/Cgas ghosts at level + fix periodic barycenters
      procedure :: sync_moments           !< Sync VF/Cliq/Cgas ghosts all levels + fix periodic barycenters
      procedure :: average_down_moments   !< Average down VF/Cliq/Cgas to coarse levels
      procedure :: reset_moments          !< Recompute VF/barycenters from PLIC
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

   !> Abstract interface for user-defined VOF boundary condition
   !> User must set VF, Cliq, Cgas, and PLIC consistently in ghost cells
   abstract interface
      subroutine vof_bc_iface(solver, lvl, time, dir, side, pVF, pCliq, pCgas, pPLIC, ilo, ihi, jlo, jhi, klo, khi)
         import :: amrvof, WP
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl       !< Level
         real(WP), intent(in) :: time     !< Current simulation time
         integer, intent(in) :: dir       !< Direction (1=x, 2=y, 3=z)
         integer, intent(in) :: side      !< Side (-1=lo, +1=hi)
         real(WP), dimension(:,:,:,:), contiguous, intent(inout) :: pVF, pCliq, pCgas, pPLIC
         integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi  !< Ghost region bounds
      end subroutine vof_bc_iface
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

   !> Dispatch tagging: tag cells near interface (3x3x3 stencil check)
   subroutine amrvof_tagging(ctx, lvl, tags, time)
      use amrex_amr_module, only: amrex_tagboxarray, amrex_mfiter, amrex_box
      use amrgrid_class, only: SETtag
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
      type(amrvof), pointer :: this
      type(amrex_tagboxarray) :: tba
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      logical :: near_interface
      integer :: i, j, k, ii, jj, kk
      real(WP) :: vf_center, vf_neighbor
      
      call c_f_pointer(ctx, this)
      
      ! Tag cells at interface using 3x3x3 stencil
      tba = tags
      call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tba%dataPtr(mfi)
         pVF => this%VF%mf(lvl)%dataptr(mfi)
         
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  vf_center = pVF(i,j,k,1)
                  near_interface = .false.
                  
                  ! Check if cell itself is interfacial
                  if (vf_center .gt. VFlo .and. vf_center .lt. VFhi) then
                     near_interface = .true.
                  else
                     ! Pure cell - check 3x3x3 neighborhood for different phase
                     stencil: do kk = k-1, k+1
                        do jj = j-1, j+1
                           do ii = i-1, i+1
                              vf_neighbor = pVF(ii,jj,kk,1)
                              ! If center is pure gas and neighbor has liquid, or vice versa
                              if (vf_center .lt. VFlo .and. vf_neighbor .gt. VFlo) then
                                 near_interface = .true.
                                 exit stencil
                              end if
                              if (vf_center .gt. VFhi .and. vf_neighbor .lt. VFhi) then
                                 near_interface = .true.
                                 exit stencil
                              end if
                           end do
                        end do
                     end do stencil
                  end if
                  
                  if (near_interface) tagarr(i,j,k,1) = SETtag
               end do
            end do
         end do
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Call user tagging if provided
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
      
      ! Initialize surface mesh for visualization
      this%smesh%name = trim(this%name)//'_plic'

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
      ! Call user init to set VF
      if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse
   subroutine on_coarse(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Allocate and interpolate moments from coarse
      call this%VF%on_coarse(this%VF, lvl, time, ba, dm)
      call this%Cliq%on_coarse(this%Cliq, lvl, time, ba, dm)
      call this%Cgas%on_coarse(this%Cgas, lvl, time, ba, dm)
      ! PLIC: just allocate, then set to trivial planes (new cells are pure, away from interface)
      call this%PLIC%reset_level(lvl, ba, dm)
      call set_trivial_plic()
      ! Old data just needs geometry
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
      ! Fill moment ghosts
      call this%fill_moments(lvl, time)
   contains
      !> Set PLIC to trivial planes based on VF
      subroutine set_trivial_plic()
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
         integer :: i, j, k
         call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
         do while (mfi%next())
            bx = mfi%tilebox()
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pPLIC(i,j,k,1:3) = 0.0_WP
                     pPLIC(i,j,k,4) = sign(1.0e10_WP, pVF(i,j,k,1) - 0.5_WP)
                  end do
               end do
            end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end subroutine set_trivial_plic
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid
   subroutine on_remake(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Remake all data (copies existing + fills from coarse for new areas)
      call this%VF%on_remake(this%VF, lvl, time, ba, dm)
      call this%Cliq%on_remake(this%Cliq, lvl, time, ba, dm)
      call this%Cgas%on_remake(this%Cgas, lvl, time, ba, dm)
      call this%PLIC%on_remake(this%PLIC, lvl, time, ba, dm)
      ! Fix PLIC for pure cells (new cells from coarse have interpolated PLIC which is wrong)
      call fix_pure_plic()
      ! Old data just needs geometry
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
      ! Fill moment ghosts
      call this%fill_moments(lvl, time)
   contains
      !> Set PLIC to trivial for pure cells
      subroutine fix_pure_plic()
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
         integer :: i, j, k
         real(WP) :: vf
         call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
         do while (mfi%next())
            bx = mfi%tilebox()
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     vf = pVF(i,j,k,1)
                     if (vf.lt.VFlo .or. vf.gt.VFhi) then
                        pPLIC(i,j,k,1:3) = 0.0_WP
                        pPLIC(i,j,k,4) = sign(1.0e10_WP, vf - 0.5_WP)
                     end if
                  end do
               end do
            end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end subroutine fix_pure_plic
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
      call this%average_down_moments(lbase)
   end subroutine post_regrid

   !> Average down VF/Cliq/Cgas from finest to lbase, then sync ghost cells
   subroutine average_down_moments(this, lbase)
      use amrex_interface, only: amrmfab_average_down_cell
      class(amrvof), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl, lb
      lb = 0; if (present(lbase)) lb = lbase
      ! Average valid cells from fine to coarse
      do lvl = this%amr%clvl()-1, lb, -1
         call amrmfab_average_down_cell(fmf=this%VF%mf(lvl+1)  , cmf=this%VF%mf(lvl)  , rr=this%amr%rref(lvl), cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_cell(fmf=this%Cliq%mf(lvl+1), cmf=this%Cliq%mf(lvl), rr=this%amr%rref(lvl), cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_cell(fmf=this%Cgas%mf(lvl+1), cmf=this%Cgas%mf(lvl), rr=this%amr%rref(lvl), cgeom=this%amr%geom(lvl))
      end do
      ! Sync ghost cells on all levels + fix periodic barycenters
      call this%sync_moments()
   end subroutine average_down_moments

   !> Sync VF/Cliq/Cgas ghosts on all levels + fix periodic barycenters
   subroutine sync_moments(this)
      class(amrvof), intent(inout) :: this
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%sync_moments_lvl(lvl)
      end do
   end subroutine sync_moments

   !> Sync VF/Cliq/Cgas ghosts at level + fix periodic barycenters
   subroutine sync_moments_lvl(this, lvl)
      use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_geometry
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pCliq, pCgas
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL
      
      ! Sync ghosts (periodic + MPI exchange)
      call this%VF%sync_lvl(lvl)
      call this%Cliq%sync_lvl(lvl)
      call this%Cgas%sync_lvl(lvl)
      
      ! Fix barycenter positions in periodic ghost cells
      geom = this%amr%geom(lvl)
      dlo = geom%domain%lo
      dhi = geom%domain%hi
      xL = this%amr%xhi - this%amr%xlo
      yL = this%amr%yhi - this%amr%ylo
      zL = this%amr%zhi - this%amr%zlo
      
      call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
      do while (mfi%next())
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         ilo = lbound(pCliq,1); ihi = ubound(pCliq,1)
         jlo = lbound(pCliq,2); jhi = ubound(pCliq,2)
         klo = lbound(pCliq,3); khi = ubound(pCliq,3)
         ! X-periodic
         if (this%amr%xper) then
            if (ilo .lt. dlo(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = ilo, dlo(1)-1
                  pCliq(ig,jg,kg,1) = pCliq(ig,jg,kg,1) - xL
                  pCgas(ig,jg,kg,1) = pCgas(ig,jg,kg,1) - xL
               end do; end do; end do
            end if
            if (ihi .gt. dhi(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = dhi(1)+1, ihi
                  pCliq(ig,jg,kg,1) = pCliq(ig,jg,kg,1) + xL
                  pCgas(ig,jg,kg,1) = pCgas(ig,jg,kg,1) + xL
               end do; end do; end do
            end if
         end if
         ! Y-periodic
         if (this%amr%yper) then
            if (jlo .lt. dlo(2)) then
               do kg = klo, khi; do jg = jlo, dlo(2)-1; do ig = ilo, ihi
                  pCliq(ig,jg,kg,2) = pCliq(ig,jg,kg,2) - yL
                  pCgas(ig,jg,kg,2) = pCgas(ig,jg,kg,2) - yL
               end do; end do; end do
            end if
            if (jhi .gt. dhi(2)) then
               do kg = klo, khi; do jg = dhi(2)+1, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,2) = pCliq(ig,jg,kg,2) + yL
                  pCgas(ig,jg,kg,2) = pCgas(ig,jg,kg,2) + yL
               end do; end do; end do
            end if
         end if
         ! Z-periodic
         if (this%amr%zper) then
            if (klo .lt. dlo(3)) then
               do kg = klo, dlo(3)-1; do jg = jlo, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,3) = pCliq(ig,jg,kg,3) - zL
                  pCgas(ig,jg,kg,3) = pCgas(ig,jg,kg,3) - zL
               end do; end do; end do
            end if
            if (khi .gt. dhi(3)) then
               do kg = dhi(3)+1, khi; do jg = jlo, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,3) = pCliq(ig,jg,kg,3) + zL
                  pCgas(ig,jg,kg,3) = pCgas(ig,jg,kg,3) + zL
               end do; end do; end do
            end if
         end if
      end do
      call this%amr%mfiter_destroy(mfi)
   end subroutine sync_moments_lvl

   !> Fill VF/Cliq/Cgas ghosts at a level (sync + periodic correction + physical BC)
   !> Also sets PLIC for BC_LIQ/BC_GAS/BC_USER. BC_REFLECT PLIC is deferred to build_plic.
   subroutine fill_moments(this, lvl, time)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      
      ! Sync VF/Cliq/Cgas ghosts + fix periodic barycenters
      call this%sync_moments_lvl(lvl)
      
      ! Sync PLIC ghosts + fix periodic PLIC
      call this%PLIC%sync_lvl(lvl)
      call correct_periodic_plic()
      
      ! Apply physical BC (sets VF/Cliq/Cgas for all types; PLIC for BC_LIQ/BC_GAS/BC_USER)
      call apply_vof_bc()
      
   contains
      
      !> Correct PLIC plane distance in periodic ghost cells
      subroutine correct_periodic_plic()
         use amrex_amr_module, only: amrex_geometry
         type(amrex_mfiter) :: mfi
         type(amrex_geometry) :: geom
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pP
         integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
         integer :: dlo(3), dhi(3)
         real(WP) :: xL, yL, zL
         
         geom = this%amr%geom(lvl)
         dlo = geom%domain%lo
         dhi = geom%domain%hi
         
         xL = this%amr%xhi - this%amr%xlo
         yL = this%amr%yhi - this%amr%ylo
         zL = this%amr%zhi - this%amr%zlo
         
         call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
         do while (mfi%next())
            pP => this%PLIC%mf(lvl)%dataptr(mfi)
            ilo = lbound(pP,1); ihi = ubound(pP,1)
            jlo = lbound(pP,2); jhi = ubound(pP,2)
            klo = lbound(pP,3); khi = ubound(pP,3)
            
            ! X-periodic
            if (this%amr%xper) then
               if (ilo .lt. dlo(1)) then
                  do kg = klo, khi; do jg = jlo, jhi; do ig = ilo, dlo(1)-1
                     pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,1)*xL
                  end do; end do; end do
               end if
               if (ihi .gt. dhi(1)) then
                  do kg = klo, khi; do jg = jlo, jhi; do ig = dhi(1)+1, ihi
                     pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,1)*xL
                  end do; end do; end do
               end if
            end if
            
            ! Y-periodic
            if (this%amr%yper) then
               if (jlo .lt. dlo(2)) then
                  do kg = klo, khi; do jg = jlo, dlo(2)-1; do ig = ilo, ihi
                     pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,2)*yL
                  end do; end do; end do
               end if
               if (jhi .gt. dhi(2)) then
                  do kg = klo, khi; do jg = dhi(2)+1, jhi; do ig = ilo, ihi
                     pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,2)*yL
                  end do; end do; end do
               end if
            end if
            
            ! Z-periodic
            if (this%amr%zper) then
               if (klo .lt. dlo(3)) then
                  do kg = klo, dlo(3)-1; do jg = jlo, jhi; do ig = ilo, ihi
                     pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,3)*zL
                  end do; end do; end do
               end if
               if (khi .gt. dhi(3)) then
                  do kg = dhi(3)+1, khi; do jg = jlo, jhi; do ig = ilo, ihi
                     pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,3)*zL
                  end do; end do; end do
               end if
            end if
            
         end do
         call this%amr%mfiter_destroy(mfi)
      end subroutine correct_periodic_plic
      
      !> Apply physical BC to VF/Cliq/Cgas/PLIC at domain boundaries
      subroutine apply_vof_bc()
         use amrex_amr_module, only: amrex_geometry
         type(amrex_mfiter) :: mfi
         type(amrex_geometry) :: geom
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
         integer :: ilo, ihi, jlo, jhi, klo, khi
         integer :: dlo(3), dhi(3)
         real(WP) :: dx, dy, dz, x_bnd
         
         geom = this%amr%geom(lvl)
         dlo = geom%domain%lo
         dhi = geom%domain%hi
         dx = this%amr%dx(lvl)
         dy = this%amr%dy(lvl)
         dz = this%amr%dz(lvl)
         
         call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
         do while (mfi%next())
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
            pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            ilo = lbound(pVF,1); ihi = ubound(pVF,1)
            jlo = lbound(pVF,2); jhi = ubound(pVF,2)
            klo = lbound(pVF,3); khi = ubound(pVF,3)
            
            ! X-boundaries
            if (.not.this%amr%xper) then
               if (ilo.lt.dlo(1)) call apply_vof_bc_face(pVF, pCliq, pCgas, pPLIC, 1, -1, this%bc_xlo, &
               &   ilo, dlo(1)-1, jlo, jhi, klo, khi, dlo(1), dx, dy, dz)
               if (ihi.gt.dhi(1)) call apply_vof_bc_face(pVF, pCliq, pCgas, pPLIC, 1, +1, this%bc_xhi, &
               &   dhi(1)+1, ihi, jlo, jhi, klo, khi, dhi(1), dx, dy, dz)
            end if
            
            ! Y-boundaries
            if (.not.this%amr%yper) then
               if (jlo.lt.dlo(2)) call apply_vof_bc_face(pVF, pCliq, pCgas, pPLIC, 2, -1, this%bc_ylo, &
               &   ilo, ihi, jlo, dlo(2)-1, klo, khi, dlo(2), dx, dy, dz)
               if (jhi.gt.dhi(2)) call apply_vof_bc_face(pVF, pCliq, pCgas, pPLIC, 2, +1, this%bc_yhi, &
               &   ilo, ihi, dhi(2)+1, jhi, klo, khi, dhi(2), dx, dy, dz)
            end if
            
            ! Z-boundaries
            if (.not.this%amr%zper) then
               if (klo.lt.dlo(3)) call apply_vof_bc_face(pVF, pCliq, pCgas, pPLIC, 3, -1, this%bc_zlo, &
               &   ilo, ihi, jlo, jhi, klo, dlo(3)-1, dlo(3), dx, dy, dz)
               if (khi.gt.dhi(3)) call apply_vof_bc_face(pVF, pCliq, pCgas, pPLIC, 3, +1, this%bc_zhi, &
               &   ilo, ihi, jlo, jhi, dhi(3)+1, khi, dhi(3), dx, dy, dz)
            end if
            
         end do
         call this%amr%mfiter_destroy(mfi)
      end subroutine apply_vof_bc
      
      !> Apply BC to VF/Cliq/Cgas on a single face; PLIC for BC_LIQ/BC_GAS/BC_USER only
      !> BC_REFLECT PLIC is handled separately in build_plic via apply_reflect_plic
      subroutine apply_vof_bc_face(pVF, pCliq, pCgas, pPLIC, dir, side, bc_type, &
      &   i1, i2, j1, j2, k1, k2, bnd, dx, dy, dz)
         real(WP), dimension(:,:,:,:), contiguous, intent(inout) :: pVF, pCliq, pCgas, pPLIC
         integer, intent(in) :: dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd
         real(WP), intent(in) :: dx, dy, dz
         integer :: ig, jg, kg, isrc, jsrc, ksrc
         real(WP), dimension(3) :: center
         
         select case (bc_type)
         
          case (BC_LIQ)
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               center = [this%amr%xlo + (real(ig,WP)+0.5_WP)*dx, &
               &         this%amr%ylo + (real(jg,WP)+0.5_WP)*dy, &
               &         this%amr%zlo + (real(kg,WP)+0.5_WP)*dz]
               pVF(ig,jg,kg,1) = 1.0_WP
               pCliq(ig,jg,kg,1:3) = center
               pCgas(ig,jg,kg,1:3) = center
               pPLIC(ig,jg,kg,1:3) = 0.0_WP
               pPLIC(ig,jg,kg,4) = 1.0e10_WP
            end do; end do; end do
            
          case (BC_GAS)
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               center = [this%amr%xlo + (real(ig,WP)+0.5_WP)*dx, &
               &         this%amr%ylo + (real(jg,WP)+0.5_WP)*dy, &
               &         this%amr%zlo + (real(kg,WP)+0.5_WP)*dz]
               pVF(ig,jg,kg,1) = 0.0_WP
               pCliq(ig,jg,kg,1:3) = center
               pCgas(ig,jg,kg,1:3) = center
               pPLIC(ig,jg,kg,1:3) = 0.0_WP
               pPLIC(ig,jg,kg,4) = -1.0e10_WP
            end do; end do; end do
            
          case (BC_REFLECT)
            ! Mirror VF/Cliq/Cgas only; PLIC handled by build_plic
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               isrc = ig; jsrc = jg; ksrc = kg
               if (dir.eq.1) isrc = 2*bnd - ig - side
               if (dir.eq.2) jsrc = 2*bnd - jg - side
               if (dir.eq.3) ksrc = 2*bnd - kg - side
               pVF(ig,jg,kg,1) = pVF(isrc,jsrc,ksrc,1)
               pCliq(ig,jg,kg,1:3) = pCliq(isrc,jsrc,ksrc,1:3)
               pCgas(ig,jg,kg,1:3) = pCgas(isrc,jsrc,ksrc,1:3)
            end do; end do; end do
            
          case (BC_USER)
            if (associated(this%user_vof_bc)) then
               call this%user_vof_bc(this, lvl, time, dir, side, pVF, pCliq, pCgas, pPLIC, i1, i2, j1, j2, k1, k2)
            end if
            
          case default
            ! Do nothing
            
         end select
         
      end subroutine apply_vof_bc_face
      
   end subroutine fill_moments

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
      call this%smesh%reset()
      nullify(this%amr)
      nullify(this%user_init)
      nullify(this%user_tagging)
   end subroutine finalize

   ! ============================================================================
   ! VOF-SPECIFIC METHODS (STUBS)
   ! ============================================================================

   !> Build PLIC reconstruction from VF and barycenters using PLICnet
   !> Also extracts PLIC polygons and accumulates them to smesh
   subroutine build_plic(this)
      use plicnet, only: get_normal, reflect_moments
      use mathtools, only: normalize
      use amrvof_geometry, only: get_plane_dist, cut_hex_polygon
      class(amrvof), intent(inout) :: this
      integer :: lvl
      real(WP) :: dx, dy, dz
      
      ! Only build at finest level
      lvl = this%amr%clvl()
      
      ! Get cell size at this level
      dx = this%amr%dx(lvl)
      dy = this%amr%dy(lvl)
      dz = this%amr%dz(lvl)
      
      ! ========== Pass 1: Compute PLIC planes ==========
      plic_reconstruction: block
         integer :: i, j, k, ii, jj, kk, direction, direction2
         real(WP), dimension(0:188) :: moments
         real(WP), dimension(3) :: normal, center, lo, hi
         real(WP) :: m000, m100, m010, m001, temp, vf_cell
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
         logical :: flip
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
         
            ! Get pointers (with ghost cells for stencil access)
            pVF   => this%VF%mf(lvl)%dataptr(mfi)
            pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
            pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
         
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
         call this%amr%mfiter_destroy(mfi)
      end block plic_reconstruction
      
      ! Sync PLIC ghost cells at finest level
      call this%PLIC%sync_lvl(lvl)
      
      ! Apply periodic correction to plane distance
      call correct_periodic_plic()
      
      ! Apply physical boundary conditions
      call apply_plic_bc()
      
      ! ========== Pass 2: Per-FAB polygon extraction ==========
      call this%smesh%reset()
      
      polygon_extraction: block
         integer :: i, j, k, ii, jj, kk
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
         real(WP), dimension(3) :: lo, hi
         real(WP), dimension(4) :: plane
         real(WP), dimension(3,8) :: hex
         real(WP) :: vf_cell
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx, gbx
         
         ! Per-FAB polygon storage (allocatable, indexed by cell)
         real(WP), dimension(:,:,:,:,:), allocatable :: polygon_local  ! (3, 6, ilo:ihi, jlo:jhi, klo:khi)
         integer, dimension(:,:,:), allocatable :: poly_nv_local       ! (ilo:ihi, jlo:jhi, klo:khi)
         real(WP), dimension(3,6) :: poly_verts
         integer :: poly_nv
         integer :: glo(3), ghi(3)
         
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            gbx = mfi%growntilebox(2)  ! Grown by 2 for 5x5x5 stencil
            pVF   => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            
            ! Get bounds for per-FAB allocation
            glo = gbx%lo
            ghi = gbx%hi
            
            ! ----- Step A: Allocate per-FAB polygon storage -----
            allocate(polygon_local(1:3, 1:6, glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3)))
            allocate(poly_nv_local(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3)))
            polygon_local = 0.0_WP
            poly_nv_local = 0
            
            ! ----- Step B: Extract polygons (grown box including ghosts) -----
            do kk = glo(3), ghi(3)
               do jj = glo(2), ghi(2)
                  do ii = glo(1), ghi(1)
                     vf_cell = pVF(ii,jj,kk,1)
                     if (vf_cell.lt.VFlo .or. vf_cell.gt.VFhi) cycle
                     
                     ! Build hex and plane
                     lo = [this%amr%xlo + real(ii,WP)*dx, this%amr%ylo + real(jj,WP)*dy, this%amr%zlo + real(kk,WP)*dz]
                     hi = [this%amr%xlo + real(ii+1,WP)*dx, this%amr%ylo + real(jj+1,WP)*dy, this%amr%zlo + real(kk+1,WP)*dz]
                     plane = [pPLIC(ii,jj,kk,1), pPLIC(ii,jj,kk,2), pPLIC(ii,jj,kk,3), pPLIC(ii,jj,kk,4)]
                     hex(:,1) = [hi(1), lo(2), lo(3)]
                     hex(:,2) = [hi(1), hi(2), lo(3)]
                     hex(:,3) = [hi(1), hi(2), hi(3)]
                     hex(:,4) = [hi(1), lo(2), hi(3)]
                     hex(:,5) = [lo(1), lo(2), lo(3)]
                     hex(:,6) = [lo(1), hi(2), lo(3)]
                     hex(:,7) = [lo(1), hi(2), hi(3)]
                     hex(:,8) = [lo(1), lo(2), hi(3)]
                     
                     call cut_hex_polygon(hex, plane, poly_nv, poly_verts)
                     
                     ! Store in per-FAB array
                     poly_nv_local(ii,jj,kk) = poly_nv
                     if (poly_nv.ge.3) then
                        polygon_local(:,1:poly_nv,ii,jj,kk) = poly_verts(:,1:poly_nv)
                     end if
                  end do
               end do
            end do
            
            ! ----- Step C: Compute curvature (valid cells, stencil access) -----
            ! TODO: curvature = f(polygon_local stencil around i,j,k)
            ! For now, skip curvature computation
            
            ! ----- Step D: Append to smesh (valid cells only) -----
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     poly_nv = poly_nv_local(i,j,k)
                     if (poly_nv.ge.3) then
                        call this%smesh%add_polygon(polygon_local(:,1:poly_nv,i,j,k), poly_nv)
                     end if
                  end do
               end do
            end do
            
            ! ----- Step E: Deallocate per-FAB storage -----
            deallocate(polygon_local, poly_nv_local)
            
         end do
         call this%amr%mfiter_destroy(mfi)
      end block polygon_extraction
      
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
         
         call this%amr%mfiter_build(lvl, mfi2, tiling=.false.)
         do while (mfi2%next())
            pP => this%PLIC%mf(lvl)%dataptr(mfi2)
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
         call this%amr%mfiter_destroy(mfi2)
         
      end subroutine correct_periodic_plic
      
      !> Apply BC_REFLECT for PLIC (BC_LIQ/BC_GAS/BC_USER already handled in fill_moments)
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
         
         call this%amr%mfiter_build(lvl, mfi3, tiling=.false.)
         do while (mfi3%next())
            pP => this%PLIC%mf(lvl)%dataptr(mfi3)
            ilo = lbound(pP,1); ihi = ubound(pP,1)
            jlo = lbound(pP,2); jhi = ubound(pP,2)
            klo = lbound(pP,3); khi = ubound(pP,3)
            
            ! X-boundaries - only BC_REFLECT needs handling here
            if (.not.this%amr%xper) then
               if (ilo.lt.dlo(1) .and. this%bc_xlo.eq.BC_REFLECT) then
                  call apply_reflect_plic(pP, 1, -1, ilo, dlo(1)-1, jlo, jhi, klo, khi, dlo(1), this%amr%xlo)
               end if
               if (ihi.gt.dhi(1) .and. this%bc_xhi.eq.BC_REFLECT) then
                  call apply_reflect_plic(pP, 1, +1, dhi(1)+1, ihi, jlo, jhi, klo, khi, dhi(1), this%amr%xhi)
               end if
            end if
            
            ! Y-boundaries
            if (.not.this%amr%yper) then
               if (jlo.lt.dlo(2) .and. this%bc_ylo.eq.BC_REFLECT) then
                  call apply_reflect_plic(pP, 2, -1, ilo, ihi, jlo, dlo(2)-1, klo, khi, dlo(2), this%amr%ylo)
               end if
               if (jhi.gt.dhi(2) .and. this%bc_yhi.eq.BC_REFLECT) then
                  call apply_reflect_plic(pP, 2, +1, ilo, ihi, dhi(2)+1, jhi, klo, khi, dhi(2), this%amr%yhi)
               end if
            end if
            
            ! Z-boundaries
            if (.not.this%amr%zper) then
               if (klo.lt.dlo(3) .and. this%bc_zlo.eq.BC_REFLECT) then
                  call apply_reflect_plic(pP, 3, -1, ilo, ihi, jlo, jhi, klo, dlo(3)-1, dlo(3), this%amr%zlo)
               end if
               if (khi.gt.dhi(3) .and. this%bc_zhi.eq.BC_REFLECT) then
                  call apply_reflect_plic(pP, 3, +1, ilo, ihi, jlo, jhi, dhi(3)+1, khi, dhi(3), this%amr%zhi)
               end if
            end if
            
         end do
         call this%amr%mfiter_destroy(mfi3)
         
      end subroutine apply_plic_bc
      
      !> Apply BC_REFLECT for PLIC on a single face
      subroutine apply_reflect_plic(pP, dir, side, i1, i2, j1, j2, k1, k2, bnd, x_bnd)
         real(WP), dimension(:,:,:,:), contiguous, intent(inout) :: pP
         integer, intent(in) :: dir, side, i1, i2, j1, j2, k1, k2, bnd
         real(WP), intent(in) :: x_bnd
         integer :: ig, jg, kg, isrc, jsrc, ksrc
         
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
         
      end subroutine apply_reflect_plic
      
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
      call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
      do while (mfi%next())
         bx = mfi%fabbox()  ! Full box including ghosts
         
         pVF   => this%VF%mf(lvl)%dataptr(mfi)
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
         
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
                  
                  ! Skip cutting for full cells (trivial PLIC with large distance)
                  if (abs(plane(4)) .ge. 1.0e9_WP) then
                     if (plane(4) .gt. 0.0_WP) then
                        pVF(i,j,k,1) = 1.0_WP  ! All liquid
                     else
                        pVF(i,j,k,1) = 0.0_WP  ! All gas
                     end if
                     pCliq(i,j,k,1:3) = cell_center
                     pCgas(i,j,k,1:3) = cell_center
                     cycle
                  end if
                  
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
      call this%amr%mfiter_destroy(mfi)
      
      ! Average down to coarse levels (uses ghost-capable average_down)
      call this%average_down_moments()
      
   end subroutine reset_moments

   !> Advect VF using staggered velocity field (U at x-faces, V at y-faces, W at z-faces)
   !> User must provide MultiFabs at finest level with >= 2 ghost cells filled
   subroutine advance_vof_stag(this, U, V, W, dt, time)
      use amrvof_geometry, only: cut_tet_vol
      use amrex_amr_module, only: amrex_multifab
      class(amrvof), intent(inout) :: this
      type(amrex_multifab), intent(in) :: U, V, W
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      integer :: lvl, i, j, k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFold, pCliqold, pCgasold, pPLICold
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      real(WP), dimension(:,:,:,:), allocatable :: Fx, Fy, Fz  ! Volume fluxes (Lvol, Gvol, Lbar(3), Gbar(3))
      real(WP) :: Lvol_old, Gvol_old, Lvol_new, Gvol_new, Lvol_flux, Gvol_flux
      real(WP), dimension(3) :: Lbar_old, Gbar_old, Lbar_new, Gbar_new, Lbar_flux, Gbar_flux
      real(WP) :: dx, dy, dz, vol, face_vel
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: ilo, ihi, jlo, jhi, klo, khi
      ! Shared variables for internal functions
      real(WP), dimension(3,9) :: face     ! 4 current + 4 projected + 1 center
      real(WP) :: dxi, dyi, dzi
      integer :: ilo_, ihi_, jlo_, jhi_, klo_, khi_  ! Velocity bounds
      ! 8-tet decomposition of 9-point polyhedron (from noirl)
      integer, dimension(4,8), parameter :: tet_map = reshape([ &
         7, 4, 3, 6, &
         6, 3, 2, 4, &
         6, 2, 1, 4, &
         7, 8, 4, 6, &
         6, 5, 8, 4, &
         6, 5, 4, 1, &
         5, 6, 8, 9, &
         6, 7, 8, 9], shape(tet_map))
      
      ! Only advect at finest level
      lvl = this%amr%clvl()
      
      ! Require >= 2 ghost cells for velocity (RK2 backtracking needs neighbors)
      if (U%nghost() .lt. 2 .or. V%nghost() .lt. 2 .or. W%nghost() .lt. 2) then
         error stop 'advance_vof_stag: velocity requires >= 2 ghost cells'
      end if
      
      ! Get cell size
      dx = this%amr%dx(lvl)
      dy = this%amr%dy(lvl)
      dz = this%amr%dz(lvl)
      vol = dx * dy * dz
      
      ! Iterate over boxes
      call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
      do while (mfi%next())
         bx = mfi%tilebox()
         
         ! Get pointers
         pVF => this%VF%mf(lvl)%dataptr(mfi)
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
         pVFold => this%VFold%mf(lvl)%dataptr(mfi)
         pCliqold => this%Cliqold%mf(lvl)%dataptr(mfi)
         pCgasold => this%Cgasold%mf(lvl)%dataptr(mfi)
         pPLICold => this%PLICold%mf(lvl)%dataptr(mfi)
         pU => U%dataptr(mfi)
         pV => V%dataptr(mfi)
         pW => W%dataptr(mfi)
         
         ilo = lbound(pVF,1); ihi = ubound(pVF,1)
         jlo = lbound(pVF,2); jhi = ubound(pVF,2)
         klo = lbound(pVF,3); khi = ubound(pVF,3)
         
         ! Allocate flux arrays (8 components: Lvol, Gvol, Lbar(3), Gbar(3))
         allocate(Fx(8,ilo:ihi+1,jlo:jhi,klo:khi)); Fx = 0.0_WP
         allocate(Fy(8,ilo:ihi,jlo:jhi+1,klo:khi)); Fy = 0.0_WP
         allocate(Fz(8,ilo:ihi,jlo:jhi,klo:khi+1)); Fz = 0.0_WP
         
         ! Compute fluxes at each face
         do k = bx%lo(3), bx%hi(3)+1
            do j = bx%lo(2), bx%hi(2)+1
               do i = bx%lo(1), bx%hi(1)+1
                  
                  ! X-flux at face (i,j,k) - U is at x-faces
                  if (i.le.bx%hi(1)+1 .and. j.le.bx%hi(2) .and. k.le.bx%hi(3)) then
                     face_vel = pU(i,j,k,1)  ! Staggered: U directly at face i
                     call compute_flux_stag(i, j, k, 1, face_vel, dt, dx, dy, dz, &
                     &   pVFold, pCliqold, pCgasold, pPLICold, pU, pV, pW, Fx(:,i,j,k))
                  end if
                  
                  ! Y-flux at face (i,j,k) - V is at y-faces
                  if (i.le.bx%hi(1) .and. j.le.bx%hi(2)+1 .and. k.le.bx%hi(3)) then
                     face_vel = pV(i,j,k,1)  ! Staggered: V directly at face j
                     call compute_flux_stag(i, j, k, 2, face_vel, dt, dx, dy, dz, &
                     &   pVFold, pCliqold, pCgasold, pPLICold, pU, pV, pW, Fy(:,i,j,k))
                  end if
                  
                  ! Z-flux at face (i,j,k) - W is at z-faces
                  if (i.le.bx%hi(1) .and. j.le.bx%hi(2) .and. k.le.bx%hi(3)+1) then
                     face_vel = pW(i,j,k,1)  ! Staggered: W directly at face k
                     call compute_flux_stag(i, j, k, 3, face_vel, dt, dx, dy, dz, &
                     &   pVFold, pCliqold, pCgasold, pPLICold, pU, pV, pW, Fz(:,i,j,k))
                  end if
                  
               end do
            end do
         end do
         
         ! Update volume moments
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  ! Old volumes
                  Lvol_old = pVFold(i,j,k,1) * vol
                  Gvol_old = (1.0_WP - pVFold(i,j,k,1)) * vol
                  Lbar_old = pCliqold(i,j,k,1:3)
                  Gbar_old = pCgasold(i,j,k,1:3)
                  
                  ! Net flux (outflow positive)
                  Lvol_flux = Fx(1,i+1,j,k) - Fx(1,i,j,k) + Fy(1,i,j+1,k) - Fy(1,i,j,k) + Fz(1,i,j,k+1) - Fz(1,i,j,k)
                  Gvol_flux = Fx(2,i+1,j,k) - Fx(2,i,j,k) + Fy(2,i,j+1,k) - Fy(2,i,j,k) + Fz(2,i,j,k+1) - Fz(2,i,j,k)
                  Lbar_flux = Fx(3:5,i+1,j,k) - Fx(3:5,i,j,k) + Fy(3:5,i,j+1,k) - Fy(3:5,i,j,k) + Fz(3:5,i,j,k+1) - Fz(3:5,i,j,k)
                  Gbar_flux = Fx(6:8,i+1,j,k) - Fx(6:8,i,j,k) + Fy(6:8,i,j+1,k) - Fy(6:8,i,j,k) + Fz(6:8,i,j,k+1) - Fz(6:8,i,j,k)
                  
                  ! New volumes
                  Lvol_new = Lvol_old - Lvol_flux
                  Gvol_new = Gvol_old - Gvol_flux
                  
                  ! New VF
                  pVF(i,j,k,1) = Lvol_new / (Lvol_new + Gvol_new)
                  
                  ! Clip and update barycenters
                  if (pVF(i,j,k,1) .lt. VFlo) then
                     pVF(i,j,k,1) = 0.0_WP
                     pCliq(i,j,k,1:3) = [this%amr%xlo + (real(i,WP)+0.5_WP)*dx, &
                     &                   this%amr%ylo + (real(j,WP)+0.5_WP)*dy, &
                     &                   this%amr%zlo + (real(k,WP)+0.5_WP)*dz]
                     pCgas(i,j,k,1:3) = pCliq(i,j,k,1:3)
                  else if (pVF(i,j,k,1) .gt. VFhi) then
                     pVF(i,j,k,1) = 1.0_WP
                     pCliq(i,j,k,1:3) = [this%amr%xlo + (real(i,WP)+0.5_WP)*dx, &
                     &                   this%amr%ylo + (real(j,WP)+0.5_WP)*dy, &
                     &                   this%amr%zlo + (real(k,WP)+0.5_WP)*dz]
                     pCgas(i,j,k,1:3) = pCliq(i,j,k,1:3)
                  else
                     ! Update barycenters from moment conservation
                     Lbar_new = (Lbar_old * Lvol_old - Lbar_flux) / Lvol_new
                     Gbar_new = (Gbar_old * Gvol_old - Gbar_flux) / Gvol_new
                     pCliq(i,j,k,1:3) = Lbar_new
                     pCgas(i,j,k,1:3) = Gbar_new
                  end if
               end do
            end do
         end do
         
         deallocate(Fx, Fy, Fz)
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Sync and apply BC
      call this%fill_moments(lvl, time)
      
   contains
      
      !> Compute flux through a face using full semi-Lagrangian RK2 (noirl-style)
      subroutine compute_flux_stag(i, j, k, dir, face_vel, dt, dx, dy, dz, &
      &   pVFold, pCliqold, pCgasold, pPLICold, pU, pV, pW, flux)
         integer, intent(in) :: i, j, k, dir
         real(WP), intent(in) :: face_vel, dt, dx, dy, dz
         real(WP), dimension(:,:,:,:), contiguous, intent(in) :: pVFold, pCliqold, pCgasold, pPLICold
         real(WP), dimension(:,:,:,:), contiguous, intent(in) :: pU, pV, pW
         real(WP), dimension(8), intent(out) :: flux
         real(WP), dimension(3,4) :: tetra    ! Single tet for cutting
         real(WP), dimension(3) :: normal, a, b, c
         real(WP) :: d, vol_expected, signed_vol, vol_liq, vol_gas
         real(WP), dimension(3) :: Lbar, Gbar, tet_Lbar, tet_Gbar
         integer :: src_i, src_j, src_k, ntet, nn
         real(WP) :: total_Lvol, total_Gvol, tet_sign_val
         real(WP), dimension(3) :: total_Lbar, total_Gbar
         
         flux = 0.0_WP
         
         ! Precompute inverse cell sizes and velocity array bounds
         dxi = 1.0_WP / dx; dyi = 1.0_WP / dy; dzi = 1.0_WP / dz
         ilo_ = lbound(pU,1); ihi_ = ubound(pU,1) - 1
         jlo_ = lbound(pV,2); jhi_ = ubound(pV,2) - 1
         klo_ = lbound(pW,3); khi_ = ubound(pW,3) - 1
         
         ! Determine upwind cell
         if (dir .eq. 1) then
            if (face_vel .ge. 0.0_WP) then
               src_i = i-1; src_j = j; src_k = k
            else
               src_i = i; src_j = j; src_k = k
            end if
         else if (dir .eq. 2) then
            if (face_vel .ge. 0.0_WP) then
               src_i = i; src_j = j-1; src_k = k
            else
               src_i = i; src_j = j; src_k = k
            end if
         else
            if (face_vel .ge. 0.0_WP) then
               src_i = i; src_j = j; src_k = k-1
            else
               src_i = i; src_j = j; src_k = k
            end if
         end if
         
         ! Skip pure cells - fast path
         if (pVFold(src_i,src_j,src_k,1) .lt. VFlo .or. pVFold(src_i,src_j,src_k,1) .gt. VFhi) then
            vol_expected = abs(face_vel) * dt
            if (dir .eq. 1) vol_expected = vol_expected * dy * dz
            if (dir .eq. 2) vol_expected = vol_expected * dx * dz
            if (dir .eq. 3) vol_expected = vol_expected * dx * dy
            if (pVFold(src_i,src_j,src_k,1) .lt. VFlo) then
               flux(2) = sign(vol_expected, face_vel)  ! Gas volume
            else
               flux(1) = sign(vol_expected, face_vel)  ! Liquid volume
            end if
            return
         end if
         
         ! Build face vertices (1-4 at current position, 5-8 projected back)
         if (dir .eq. 1) then
            face(:,1) = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
            face(:,2) = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
            face(:,3) = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
            face(:,4) = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
            vol_expected = -face_vel * dt * dy * dz
         else if (dir .eq. 2) then
            face(:,1) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
            face(:,2) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
            face(:,3) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
            face(:,4) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
            vol_expected = -face_vel * dt * dx * dz
         else
            face(:,1) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k,WP)*dz]
            face(:,2) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k,WP)*dz]
            face(:,3) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k,WP)*dz]
            face(:,4) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k,WP)*dz]
            vol_expected = -face_vel * dt * dx * dy
         end if
         
         ! Project vertices back in time with RK2
         face(:,5) = project(face(:,1), -dt)
         face(:,6) = project(face(:,2), -dt)
         face(:,7) = project(face(:,3), -dt)
         face(:,8) = project(face(:,4), -dt)
         
         ! Center point (will be adjusted by volume_correct)
         face(:,9) = 0.25_WP * (face(:,5) + face(:,6) + face(:,7) + face(:,8))
         
         ! Apply volume correction for exact flux
         if (dir .eq. 1) then
            call volume_correct_x(vol_expected)
         else if (dir .eq. 2) then
            call volume_correct_y(vol_expected)
         else
            call volume_correct_z(vol_expected)
         end if
         
         ! Get PLIC from upwind cell
         normal = pPLICold(src_i,src_j,src_k,1:3)
         d = pPLICold(src_i,src_j,src_k,4)
         
         ! Cut each of the 8 tets and accumulate flux
         total_Lvol = 0.0_WP; total_Gvol = 0.0_WP
         total_Lbar = 0.0_WP; total_Gbar = 0.0_WP
         
         do ntet = 1, 8
            ! Build tet from face vertices
            do nn = 1, 4
               tetra(:,nn) = face(:,tet_map(nn,ntet))
            end do
            
            ! Compute signed tet volume
            a = tetra(:,1) - tetra(:,4)
            b = tetra(:,2) - tetra(:,4)
            c = tetra(:,3) - tetra(:,4)
            signed_vol = (a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
            tet_sign_val = sign(1.0_WP, -signed_vol)
            
            ! Cut tet with PLIC
            call cut_tet_vol(tetra, [normal(1), normal(2), normal(3), d], vol_liq, vol_gas, tet_Lbar, tet_Gbar)
            
            ! Accumulate with proper sign
            total_Lvol = total_Lvol + tet_sign_val * vol_liq
            total_Gvol = total_Gvol + tet_sign_val * vol_gas
            total_Lbar = total_Lbar + tet_sign_val * vol_liq * tet_Lbar
            total_Gbar = total_Gbar + tet_sign_val * vol_gas * tet_Gbar
         end do
         
         ! Normalize barycenters
         if (abs(total_Lvol) .gt. epsilon(1.0_WP)) total_Lbar = total_Lbar / total_Lvol
         if (abs(total_Gvol) .gt. epsilon(1.0_WP)) total_Gbar = total_Gbar / total_Gvol
         
         ! Output flux (sign already handled in tet cutting)
         flux(1) = total_Lvol
         flux(2) = total_Gvol
         flux(3:5) = total_Lbar * total_Lvol
         flux(6:8) = total_Gbar * total_Gvol
         
      end subroutine compute_flux_stag
      
      !> RK2 vertex projection back in time
      function project(p1, mydt) result(p2)
         real(WP), dimension(3), intent(in) :: p1
         real(WP), intent(in) :: mydt
         real(WP), dimension(3) :: p2, vel1, vel_mid
         vel1 = interp_velocity_stag(p1)
         p2 = p1 + mydt * vel1
         vel_mid = interp_velocity_stag(0.5_WP * (p1 + p2))
         p2 = p1 + mydt * vel_mid
      end function project
      
      !> Trilinear interpolation of staggered velocity with index clipping
      function interp_velocity_stag(pos) result(vel)
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel
         integer :: ipc, jpc, kpc, ipu, jpv, kpw
         real(WP) :: wxc1, wyc1, wzc1, wxc2, wyc2, wzc2
         real(WP) :: wxu1, wyv1, wzw1, wxu2, wyv2, wzw2
         
         ! Cell-centered indices and weights
         ipc = max(ilo_, min(ihi_, floor((pos(1) - this%amr%xlo) * dxi)))
         jpc = max(jlo_, min(jhi_, floor((pos(2) - this%amr%ylo) * dyi)))
         kpc = max(klo_, min(khi_, floor((pos(3) - this%amr%zlo) * dzi)))
         wxc1 = max(0.0_WP, min(1.0_WP, (pos(1) - (this%amr%xlo + (real(ipc,WP)+0.5_WP)*dx)) * dxi + 0.5_WP))
         wyc1 = max(0.0_WP, min(1.0_WP, (pos(2) - (this%amr%ylo + (real(jpc,WP)+0.5_WP)*dy)) * dyi + 0.5_WP))
         wzc1 = max(0.0_WP, min(1.0_WP, (pos(3) - (this%amr%zlo + (real(kpc,WP)+0.5_WP)*dz)) * dzi + 0.5_WP))
         wxc2 = 1.0_WP - wxc1; wyc2 = 1.0_WP - wyc1; wzc2 = 1.0_WP - wzc1
         
         ! Face-centered indices and weights
         ipu = max(ilo_, min(ihi_, floor((pos(1) - this%amr%xlo) * dxi)))
         jpv = max(jlo_, min(jhi_, floor((pos(2) - this%amr%ylo) * dyi)))
         kpw = max(klo_, min(khi_, floor((pos(3) - this%amr%zlo) * dzi)))
         wxu1 = max(0.0_WP, min(1.0_WP, (pos(1) - (this%amr%xlo + real(ipu,WP)*dx)) * dxi))
         wyv1 = max(0.0_WP, min(1.0_WP, (pos(2) - (this%amr%ylo + real(jpv,WP)*dy)) * dyi))
         wzw1 = max(0.0_WP, min(1.0_WP, (pos(3) - (this%amr%zlo + real(kpw,WP)*dz)) * dzi))
         wxu2 = 1.0_WP - wxu1; wyv2 = 1.0_WP - wyv1; wzw2 = 1.0_WP - wzw1
         
         ! Trilinear interpolation of U (staggered in x)
         vel(1) = wzc1*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc+1,1)+wxu2*pU(ipu,jpc+1,kpc+1,1)) + &
         &             wyc2*(wxu1*pU(ipu+1,jpc  ,kpc+1,1)+wxu2*pU(ipu,jpc  ,kpc+1,1))) + &
         &        wzc2*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc  ,1)+wxu2*pU(ipu,jpc+1,kpc  ,1)) + &
         &             wyc2*(wxu1*pU(ipu+1,jpc  ,kpc  ,1)+wxu2*pU(ipu,jpc  ,kpc  ,1)))
         ! Trilinear interpolation of V (staggered in y)
         vel(2) = wzc1*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc+1,1)+wxc2*pV(ipc,jpv+1,kpc+1,1)) + &
         &             wyv2*(wxc1*pV(ipc+1,jpv  ,kpc+1,1)+wxc2*pV(ipc,jpv  ,kpc+1,1))) + &
         &        wzc2*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc  ,1)+wxc2*pV(ipc,jpv+1,kpc  ,1)) + &
         &             wyv2*(wxc1*pV(ipc+1,jpv  ,kpc  ,1)+wxc2*pV(ipc,jpv  ,kpc  ,1)))
         ! Trilinear interpolation of W (staggered in z)
         vel(3) = wzw1*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw+1,1)+wxc2*pW(ipc,jpc+1,kpw+1,1)) + &
         &             wyc2*(wxc1*pW(ipc+1,jpc  ,kpw+1,1)+wxc2*pW(ipc,jpc  ,kpw+1,1))) + &
         &        wzw2*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw  ,1)+wxc2*pW(ipc,jpc+1,kpw  ,1)) + &
         &             wyc2*(wxc1*pW(ipc+1,jpc  ,kpw  ,1)+wxc2*pW(ipc,jpc  ,kpw  ,1)))
      end function interp_velocity_stag
      
      !> Volume correction for X-flux (analytical from noirl)
      subroutine volume_correct_x(volume)
         real(WP), intent(inout) :: volume
         real(WP), dimension(3) :: va, vb, vc
         integer :: ntet_
         ! Compute current volume mismatch using first 6 tets
         do ntet_ = 1, 6
            va = face(:,tet_map(1,ntet_)) - face(:,tet_map(4,ntet_))
            vb = face(:,tet_map(2,ntet_)) - face(:,tet_map(4,ntet_))
            vc = face(:,tet_map(3,ntet_)) - face(:,tet_map(4,ntet_))
            volume = volume + (va(1)*(vb(2)*vc(3)-vc(2)*vb(3)) - va(2)*(vb(1)*vc(3)-vc(1)*vb(3)) + va(3)*(vb(1)*vc(2)-vc(1)*vb(2))) / 6.0_WP
         end do
         ! Analytical correction for x-component of 9th point
         face(1,9) = (-6.0_WP*volume + face(1,5)*((face(2,8)-face(2,9))*(face(3,6)-face(3,9))-(face(2,6)-face(2,9))*(face(3,8)-face(3,9))) + &
         &   face(2,5)*((face(3,8)-face(3,9))*face(1,6)-(face(3,6)-face(3,9))*face(1,8)) + &
         &   face(2,9)*((face(3,6)-face(3,9))*face(1,8)-(face(3,8)-face(3,9))*face(1,6)) + &
         &   face(3,5)*((face(2,6)-face(2,9))*face(1,8)-(face(2,8)-face(2,9))*face(1,6)) + &
         &   face(3,9)*((face(2,8)-face(2,9))*face(1,6)-(face(2,6)-face(2,9))*face(1,8)) + &
         &   face(1,6)*((face(2,8)-face(2,9))*(face(3,7)-face(3,9))-(face(2,7)-face(2,9))*(face(3,8)-face(3,9))) + &
         &   face(2,6)*((face(3,8)-face(3,9))*face(1,7)-(face(3,7)-face(3,9))*face(1,8)) + &
         &   face(2,9)*((face(3,7)-face(3,9))*face(1,8)-(face(3,8)-face(3,9))*face(1,7)) + &
         &   face(3,6)*((face(2,7)-face(2,9))*face(1,8)-(face(2,8)-face(2,9))*face(1,7)) + &
         &   face(3,9)*((face(2,8)-face(2,9))*face(1,7)-(face(2,7)-face(2,9))*face(1,8))) / &
         &  (-(face(2,6)-face(2,9))*(face(3,8)-face(3,9))+(face(2,8)-face(2,9))*(face(3,6)-face(3,9)) - &
         &   face(2,5)*(face(3,6)-face(3,9))+face(2,5)*(face(3,8)-face(3,9))+face(2,9)*(face(3,6)-face(3,9))-face(2,9)*(face(3,8)-face(3,9)) - &
         &   face(3,5)*(face(2,8)-face(2,9))+face(3,5)*(face(2,6)-face(2,9))+face(3,9)*(face(2,8)-face(2,9))-face(3,9)*(face(2,6)-face(2,9)) - &
         &   (face(2,7)-face(2,9))*(face(3,8)-face(3,9))+(face(2,8)-face(2,9))*(face(3,7)-face(3,9)) - &
         &   face(2,6)*(face(3,7)-face(3,9))+face(2,6)*(face(3,8)-face(3,9))+face(2,9)*(face(3,7)-face(3,9))-face(2,9)*(face(3,8)-face(3,9)) - &
         &   face(3,6)*(face(2,8)-face(2,9))+face(3,6)*(face(2,7)-face(2,9))+face(3,9)*(face(2,8)-face(2,9))-face(3,9)*(face(2,7)-face(2,9)))
         face(2,9) = 0.25_WP * sum(face(2,5:8))
         face(3,9) = 0.25_WP * sum(face(3,5:8))
      end subroutine volume_correct_x
      
      !> Volume correction for Y-flux (analytical from noirl)
      subroutine volume_correct_y(volume)
         real(WP), intent(inout) :: volume
         real(WP), dimension(3) :: va, vb, vc
         integer :: ntet_
         ! Compute current volume mismatch (negative sign for Y)
         do ntet_ = 1, 6
            va = face(:,tet_map(1,ntet_)) - face(:,tet_map(4,ntet_))
            vb = face(:,tet_map(2,ntet_)) - face(:,tet_map(4,ntet_))
            vc = face(:,tet_map(3,ntet_)) - face(:,tet_map(4,ntet_))
            volume = volume - (va(1)*(vb(2)*vc(3)-vc(2)*vb(3)) - va(2)*(vb(1)*vc(3)-vc(1)*vb(3)) + va(3)*(vb(1)*vc(2)-vc(1)*vb(2))) / 6.0_WP
         end do
         ! Analytical correction for y-component of 9th point
         face(1,9) = 0.25_WP * sum(face(1,5:8))
         face(2,9) = (6.0_WP*volume + face(1,5)*((face(3,6)-face(3,9))*face(2,8)-(face(3,8)-face(3,9))*face(2,6)) + &
         &   face(1,9)*((face(3,8)-face(3,9))*face(2,6)-(face(3,6)-face(3,9))*face(2,8)) + &
         &   face(2,5)*((face(3,8)-face(3,9))*(face(1,6)-face(1,9))-(face(3,6)-face(3,9))*(face(1,8)-face(1,9))) + &
         &   face(3,5)*((face(1,8)-face(1,9))*face(2,6)-(face(1,6)-face(1,9))*face(2,8)) + &
         &   face(3,9)*((face(1,6)-face(1,9))*face(2,8)-(face(1,8)-face(1,9))*face(2,6)) + &
         &   face(1,6)*((face(3,7)-face(3,9))*face(2,8)-(face(3,8)-face(3,9))*face(2,7)) + &
         &   face(1,9)*((face(3,8)-face(3,9))*face(2,7)-(face(3,7)-face(3,9))*face(2,8)) + &
         &   face(2,6)*((face(3,8)-face(3,9))*(face(1,7)-face(1,9))-(face(3,7)-face(3,9))*(face(1,8)-face(1,9))) + &
         &   face(3,6)*((face(1,8)-face(1,9))*face(2,7)-(face(1,7)-face(1,9))*face(2,8)) + &
         &   face(3,9)*((face(1,7)-face(1,9))*face(2,8)-(face(1,8)-face(1,9))*face(2,7))) / &
         &  (face(1,5)*((face(3,6)-face(3,9))-(face(3,8)-face(3,9))) + &
         &   face(1,9)*((face(3,8)-face(3,9))-(face(3,6)-face(3,9))) + &
         &   ((face(3,8)-face(3,9))*(face(1,6)-face(1,9))-(face(3,6)-face(3,9))*(face(1,8)-face(1,9))) + &
         &   face(3,5)*((face(1,8)-face(1,9))-(face(1,6)-face(1,9)))+face(3,9)*((face(1,6)-face(1,9))-(face(1,8)-face(1,9))) + &
         &   face(1,9)*((face(3,8)-face(3,9))-(face(3,7)-face(3,9))) + &
         &   ((face(3,8)-face(3,9))*(face(1,7)-face(1,9))-(face(3,7)-face(3,9))*(face(1,8)-face(1,9))) + &
         &   face(3,6)*((face(1,8)-face(1,9))-(face(1,7)-face(1,9)))+face(3,9)*((face(1,7)-face(1,9))-(face(1,8)-face(1,9))))
         face(3,9) = 0.25_WP * sum(face(3,5:8))
      end subroutine volume_correct_y
      
      !> Volume correction for Z-flux (analytical from noirl)
      subroutine volume_correct_z(volume)
         real(WP), intent(inout) :: volume
         real(WP), dimension(3) :: va, vb, vc
         integer :: ntet_
         ! Compute current volume mismatch
         do ntet_ = 1, 6
            va = face(:,tet_map(1,ntet_)) - face(:,tet_map(4,ntet_))
            vb = face(:,tet_map(2,ntet_)) - face(:,tet_map(4,ntet_))
            vc = face(:,tet_map(3,ntet_)) - face(:,tet_map(4,ntet_))
            volume = volume + (va(1)*(vb(2)*vc(3)-vc(2)*vb(3)) - va(2)*(vb(1)*vc(3)-vc(1)*vb(3)) + va(3)*(vb(1)*vc(2)-vc(1)*vb(2))) / 6.0_WP
         end do
         ! Analytical correction for z-component
         face(1,9) = 0.25_WP * sum(face(1,5:8))
         face(2,9) = 0.25_WP * sum(face(2,5:8))
         face(3,9) = (6.0_WP*volume + &
         &   face(1,5)*face(2,6)*face(3,8) - face(1,5)*face(3,6)*face(2,8) - &
         &   face(2,5)*face(1,6)*face(3,8) + face(2,5)*face(3,6)*face(1,8) + &
         &   face(3,5)*face(1,6)*face(2,8) - face(3,5)*face(2,6)*face(1,8) + &
         &   face(1,5)*face(3,6)*face(2,5) - face(2,5)*face(3,6)*face(1,5) - &
         &   face(3,5)*face(1,6)*face(2,5) + face(3,5)*face(2,6)*face(1,5) + &
         &   face(1,6)*face(2,7)*face(3,8) - face(1,6)*face(3,7)*face(2,8) - &
         &   face(2,6)*face(1,7)*face(3,8) + face(2,6)*face(3,7)*face(1,8) + &
         &   face(3,6)*face(1,7)*face(2,8) - face(3,6)*face(2,7)*face(1,8) - &
         &   face(1,5)*face(3,8)*face(2,5) + face(2,5)*face(3,8)*face(1,5) + &
         &   face(3,5)*face(1,8)*face(2,5) - face(3,5)*face(2,8)*face(1,5) + &
         &   face(1,6)*face(3,7)*face(2,5) - face(2,6)*face(3,7)*face(1,5) - &
         &   face(3,6)*face(1,7)*face(2,5) + face(3,6)*face(2,7)*face(1,5) + &
         &   face(1,7)*face(3,8)*face(2,5) - face(2,7)*face(3,8)*face(1,5) - &
         &   face(3,7)*face(1,8)*face(2,5) + face(3,7)*face(2,8)*face(1,5)) / &
         &  (face(1,5)*face(2,6) - face(2,5)*face(1,6) - face(1,5)*face(2,8) + face(2,5)*face(1,8) + &
         &   face(1,6)*face(2,7) - face(2,6)*face(1,7) + face(1,7)*face(2,8) - face(2,7)*face(1,8))
      end subroutine volume_correct_z
      
   end subroutine advance_vof_stag

   !> Advect VF using collocated velocity field (U, V, W all at cell centers)
   !> User must provide MultiFabs at finest level with >= 2 ghost cells filled
   subroutine advance_vof_col(this, U, V, W, dt, time)
      use amrvof_geometry, only: cut_tet_vol
      use amrex_amr_module, only: amrex_multifab
      class(amrvof), intent(inout) :: this
      type(amrex_multifab), intent(in) :: U, V, W
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      integer :: lvl, i, j, k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFold, pCliqold, pCgasold, pPLICold
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      real(WP), dimension(:,:,:,:), allocatable :: Fx, Fy, Fz
      real(WP) :: Lvol_old, Gvol_old, Lvol_new, Gvol_new, Lvol_flux, Gvol_flux
      real(WP), dimension(3) :: Lbar_old, Gbar_old, Lbar_new, Gbar_new, Lbar_flux, Gbar_flux
      real(WP) :: dx, dy, dz, vol, face_vel
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: ilo, ihi, jlo, jhi, klo, khi
      ! Shared variables for internal functions
      real(WP), dimension(3,9) :: face_col    ! 4 current + 4 projected + 1 center
      real(WP) :: dxi, dyi, dzi
      integer :: ilo_, ihi_, jlo_, jhi_, klo_, khi_  ! Velocity bounds
      ! 8-tet decomposition of 9-point polyhedron (from noirl)
      integer, dimension(4,8), parameter :: tet_map = reshape([ &
         7, 4, 3, 6, &
         6, 3, 2, 4, &
         6, 2, 1, 4, &
         7, 8, 4, 6, &
         6, 5, 8, 4, &
         6, 5, 4, 1, &
         5, 6, 8, 9, &
         6, 7, 8, 9], shape(tet_map))
      
      ! Only advect at finest level
      lvl = this%amr%clvl()
      
      ! Require >= 2 ghost cells for velocity (RK2 backtracking needs neighbors)
      if (U%nghost() .lt. 2 .or. V%nghost() .lt. 2 .or. W%nghost() .lt. 2) then
         error stop 'advance_vof_col: velocity requires >= 2 ghost cells'
      end if
      
      dx = this%amr%dx(lvl)
      dy = this%amr%dy(lvl)
      dz = this%amr%dz(lvl)
      vol = dx * dy * dz
      
      ! Iterate over boxes
      call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
      do while (mfi%next())
         bx = mfi%tilebox()
         
         pVF => this%VF%mf(lvl)%dataptr(mfi)
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
         pVFold => this%VFold%mf(lvl)%dataptr(mfi)
         pCliqold => this%Cliqold%mf(lvl)%dataptr(mfi)
         pCgasold => this%Cgasold%mf(lvl)%dataptr(mfi)
         pPLICold => this%PLICold%mf(lvl)%dataptr(mfi)
         pU => U%dataptr(mfi)
         pV => V%dataptr(mfi)
         pW => W%dataptr(mfi)
         
         ilo = lbound(pVF,1); ihi = ubound(pVF,1)
         jlo = lbound(pVF,2); jhi = ubound(pVF,2)
         klo = lbound(pVF,3); khi = ubound(pVF,3)
         
         allocate(Fx(8,ilo:ihi+1,jlo:jhi,klo:khi)); Fx = 0.0_WP
         allocate(Fy(8,ilo:ihi,jlo:jhi+1,klo:khi)); Fy = 0.0_WP
         allocate(Fz(8,ilo:ihi,jlo:jhi,klo:khi+1)); Fz = 0.0_WP
         
         ! Compute fluxes using collocated velocities
         do k = bx%lo(3), bx%hi(3)+1
            do j = bx%lo(2), bx%hi(2)+1
               do i = bx%lo(1), bx%hi(1)+1
                  
                  ! X-flux: average cell-centered U
                  if (i.le.bx%hi(1)+1 .and. j.le.bx%hi(2) .and. k.le.bx%hi(3)) then
                     face_vel = 0.5_WP*(pU(i-1,j,k,1) + pU(i,j,k,1))
                     call compute_flux_col(i, j, k, 1, face_vel, dt, dx, dy, dz, &
                     &   pVFold, pCliqold, pCgasold, pPLICold, pU, pV, pW, Fx(:,i,j,k))
                  end if
                  
                  ! Y-flux: average cell-centered V
                  if (i.le.bx%hi(1) .and. j.le.bx%hi(2)+1 .and. k.le.bx%hi(3)) then
                     face_vel = 0.5_WP*(pV(i,j-1,k,1) + pV(i,j,k,1))
                     call compute_flux_col(i, j, k, 2, face_vel, dt, dx, dy, dz, &
                     &   pVFold, pCliqold, pCgasold, pPLICold, pU, pV, pW, Fy(:,i,j,k))
                  end if
                  
                  ! Z-flux: average cell-centered W
                  if (i.le.bx%hi(1) .and. j.le.bx%hi(2) .and. k.le.bx%hi(3)+1) then
                     face_vel = 0.5_WP*(pW(i,j,k-1,1) + pW(i,j,k,1))
                     call compute_flux_col(i, j, k, 3, face_vel, dt, dx, dy, dz, &
                     &   pVFold, pCliqold, pCgasold, pPLICold, pU, pV, pW, Fz(:,i,j,k))
                  end if
                  
               end do
            end do
         end do
         
         ! Update volume moments (same as staggered version)
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  Lvol_old = pVFold(i,j,k,1) * vol
                  Gvol_old = (1.0_WP - pVFold(i,j,k,1)) * vol
                  Lbar_old = pCliqold(i,j,k,1:3)
                  Gbar_old = pCgasold(i,j,k,1:3)
                  
                  Lvol_flux = Fx(1,i+1,j,k) - Fx(1,i,j,k) + Fy(1,i,j+1,k) - Fy(1,i,j,k) + Fz(1,i,j,k+1) - Fz(1,i,j,k)
                  Gvol_flux = Fx(2,i+1,j,k) - Fx(2,i,j,k) + Fy(2,i,j+1,k) - Fy(2,i,j,k) + Fz(2,i,j,k+1) - Fz(2,i,j,k)
                  Lbar_flux = Fx(3:5,i+1,j,k) - Fx(3:5,i,j,k) + Fy(3:5,i,j+1,k) - Fy(3:5,i,j,k) + Fz(3:5,i,j,k+1) - Fz(3:5,i,j,k)
                  Gbar_flux = Fx(6:8,i+1,j,k) - Fx(6:8,i,j,k) + Fy(6:8,i,j+1,k) - Fy(6:8,i,j,k) + Fz(6:8,i,j,k+1) - Fz(6:8,i,j,k)
                  
                  Lvol_new = Lvol_old - Lvol_flux
                  Gvol_new = Gvol_old - Gvol_flux
                  
                  pVF(i,j,k,1) = Lvol_new / (Lvol_new + Gvol_new)
                  
                  if (pVF(i,j,k,1) .lt. VFlo) then
                     pVF(i,j,k,1) = 0.0_WP
                     pCliq(i,j,k,1:3) = [this%amr%xlo + (real(i,WP)+0.5_WP)*dx, &
                     &                   this%amr%ylo + (real(j,WP)+0.5_WP)*dy, &
                     &                   this%amr%zlo + (real(k,WP)+0.5_WP)*dz]
                     pCgas(i,j,k,1:3) = pCliq(i,j,k,1:3)
                  else if (pVF(i,j,k,1) .gt. VFhi) then
                     pVF(i,j,k,1) = 1.0_WP
                     pCliq(i,j,k,1:3) = [this%amr%xlo + (real(i,WP)+0.5_WP)*dx, &
                     &                   this%amr%ylo + (real(j,WP)+0.5_WP)*dy, &
                     &                   this%amr%zlo + (real(k,WP)+0.5_WP)*dz]
                     pCgas(i,j,k,1:3) = pCliq(i,j,k,1:3)
                  else
                     Lbar_new = (Lbar_old * Lvol_old - Lbar_flux) / Lvol_new
                     Gbar_new = (Gbar_old * Gvol_old - Gbar_flux) / Gvol_new
                     pCliq(i,j,k,1:3) = Lbar_new
                     pCgas(i,j,k,1:3) = Gbar_new
                  end if
               end do
            end do
         end do
         
         deallocate(Fx, Fy, Fz)
      end do
      call this%amr%mfiter_destroy(mfi)
      
      call this%fill_moments(lvl, time)
      
   contains
      
      !> Compute flux through a face using full semi-Lagrangian RK2 with collocated velocities
      subroutine compute_flux_col(i, j, k, dir, face_vel, dt, dx, dy, dz, &
      &   pVFold, pCliqold, pCgasold, pPLICold, pU, pV, pW, flux)
         integer, intent(in) :: i, j, k, dir
         real(WP), intent(in) :: face_vel, dt, dx, dy, dz
         real(WP), dimension(:,:,:,:), contiguous, intent(in) :: pVFold, pCliqold, pCgasold, pPLICold
         real(WP), dimension(:,:,:,:), contiguous, intent(in) :: pU, pV, pW
         real(WP), dimension(8), intent(out) :: flux
         ! Use face_col from parent scope for flux polyhedron
         real(WP), dimension(3,4) :: tetra
         real(WP), dimension(3) :: normal, a, b, c
         real(WP) :: d, vol_expected, signed_vol, vol_liq, vol_gas
         real(WP), dimension(3) :: Lbar, Gbar, tet_Lbar, tet_Gbar
         integer :: src_i, src_j, src_k, ntet, nn
         integer :: ilo_, ihi_, jlo_, jhi_, klo_, khi_
         real(WP) :: dxi, dyi, dzi
         real(WP) :: total_Lvol, total_Gvol, tet_sign_val
         real(WP), dimension(3) :: total_Lbar, total_Gbar
         integer, dimension(4,8), parameter :: tet_map = reshape([ &
            7, 4, 3, 6, 6, 3, 2, 4, 6, 2, 1, 4, 7, 8, 4, 6, &
            6, 5, 8, 4, 6, 5, 4, 1, 5, 6, 8, 9, 6, 7, 8, 9], shape(tet_map))
         
         flux = 0.0_WP
         
         dxi = 1.0_WP / dx; dyi = 1.0_WP / dy; dzi = 1.0_WP / dz
         ilo_ = lbound(pU,1); ihi_ = ubound(pU,1) - 1
         jlo_ = lbound(pU,2); jhi_ = ubound(pU,2) - 1
         klo_ = lbound(pU,3); khi_ = ubound(pU,3) - 1
         
         ! Determine upwind cell
         if (dir .eq. 1) then
            if (face_vel .ge. 0.0_WP) then; src_i = i-1; src_j = j; src_k = k
            else; src_i = i; src_j = j; src_k = k; end if
         else if (dir .eq. 2) then
            if (face_vel .ge. 0.0_WP) then; src_i = i; src_j = j-1; src_k = k
            else; src_i = i; src_j = j; src_k = k; end if
         else
            if (face_vel .ge. 0.0_WP) then; src_i = i; src_j = j; src_k = k-1
            else; src_i = i; src_j = j; src_k = k; end if
         end if
         
         ! Skip pure cells
         if (pVFold(src_i,src_j,src_k,1) .lt. VFlo .or. pVFold(src_i,src_j,src_k,1) .gt. VFhi) then
            vol_expected = abs(face_vel) * dt
            if (dir .eq. 1) vol_expected = vol_expected * dy * dz
            if (dir .eq. 2) vol_expected = vol_expected * dx * dz
            if (dir .eq. 3) vol_expected = vol_expected * dx * dy
            if (pVFold(src_i,src_j,src_k,1) .lt. VFlo) then
               flux(2) = sign(vol_expected, face_vel)
            else
               flux(1) = sign(vol_expected, face_vel)
            end if
            return
         end if
         
         ! Build face vertices
         if (dir .eq. 1) then
            face_col(:,1) = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
            face_col(:,2) = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
            face_col(:,3) = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
            face_col(:,4) = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
            vol_expected = -face_vel * dt * dy * dz
         else if (dir .eq. 2) then
            face_col(:,1) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
            face_col(:,2) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
            face_col(:,3) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
            face_col(:,4) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
            vol_expected = -face_vel * dt * dx * dz
         else
            face_col(:,1) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k,WP)*dz]
            face_col(:,2) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k,WP)*dz]
            face_col(:,3) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k,WP)*dz]
            face_col(:,4) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k,WP)*dz]
            vol_expected = -face_vel * dt * dx * dy
         end if
         
         ! Project vertices back with RK2
         face_col(:,5) = project_col(face_col(:,1), -dt)
         face_col(:,6) = project_col(face_col(:,2), -dt)
         face_col(:,7) = project_col(face_col(:,3), -dt)
         face_col(:,8) = project_col(face_col(:,4), -dt)
         face_col(:,9) = 0.25_WP * (face_col(:,5) + face_col(:,6) + face_col(:,7) + face_col(:,8))
         
         ! Volume correction
         if (dir .eq. 1) then; call volume_correct_x_col(vol_expected)
         else if (dir .eq. 2) then; call volume_correct_y_col(vol_expected)
         else; call volume_correct_z_col(vol_expected); end if
         
         normal = pPLICold(src_i,src_j,src_k,1:3)
         d = pPLICold(src_i,src_j,src_k,4)
         
         ! Cut each of 8 tets
         total_Lvol = 0.0_WP; total_Gvol = 0.0_WP
         total_Lbar = 0.0_WP; total_Gbar = 0.0_WP
         
         do ntet = 1, 8
            do nn = 1, 4; tetra(:,nn) = face_col(:,tet_map(nn,ntet)); end do
            a = tetra(:,1) - tetra(:,4); b = tetra(:,2) - tetra(:,4); c = tetra(:,3) - tetra(:,4)
            signed_vol = (a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
            tet_sign_val = sign(1.0_WP, -signed_vol)
            call cut_tet_vol(tetra, [normal(1), normal(2), normal(3), d], vol_liq, vol_gas, tet_Lbar, tet_Gbar)
            total_Lvol = total_Lvol + tet_sign_val * vol_liq
            total_Gvol = total_Gvol + tet_sign_val * vol_gas
            total_Lbar = total_Lbar + tet_sign_val * vol_liq * tet_Lbar
            total_Gbar = total_Gbar + tet_sign_val * vol_gas * tet_Gbar
         end do
         
         if (abs(total_Lvol) .gt. epsilon(1.0_WP)) total_Lbar = total_Lbar / total_Lvol
         if (abs(total_Gvol) .gt. epsilon(1.0_WP)) total_Gbar = total_Gbar / total_Gvol
         
         flux(1) = total_Lvol; flux(2) = total_Gvol
         flux(3:5) = total_Lbar * total_Lvol; flux(6:8) = total_Gbar * total_Gvol
         
      end subroutine compute_flux_col
      
      !> RK2 vertex projection back in time (collocated version)
      function project_col(p1, mydt) result(p2)
         real(WP), dimension(3), intent(in) :: p1
         real(WP), intent(in) :: mydt
         real(WP), dimension(3) :: p2, v1, vm
         v1 = interp_velocity_col(p1); p2 = p1 + mydt * v1
         vm = interp_velocity_col(0.5_WP * (p1 + p2)); p2 = p1 + mydt * vm
      end function project_col
      
      !> Trilinear interpolation of collocated (cell-centered) velocity
      function interp_velocity_col(pos) result(vel)
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel
         integer :: ip, jp, kp
         real(WP) :: wx1, wy1, wz1, wx2, wy2, wz2
         ! Cell-centered indices
         ip = max(ilo_, min(ihi_, floor((pos(1) - this%amr%xlo) * dxi)))
         jp = max(jlo_, min(jhi_, floor((pos(2) - this%amr%ylo) * dyi)))
         kp = max(klo_, min(khi_, floor((pos(3) - this%amr%zlo) * dzi)))
         wx1 = max(0.0_WP, min(1.0_WP, (pos(1) - (this%amr%xlo + (real(ip,WP)+0.5_WP)*dx)) * dxi + 0.5_WP))
         wy1 = max(0.0_WP, min(1.0_WP, (pos(2) - (this%amr%ylo + (real(jp,WP)+0.5_WP)*dy)) * dyi + 0.5_WP))
         wz1 = max(0.0_WP, min(1.0_WP, (pos(3) - (this%amr%zlo + (real(kp,WP)+0.5_WP)*dz)) * dzi + 0.5_WP))
         wx2 = 1.0_WP - wx1; wy2 = 1.0_WP - wy1; wz2 = 1.0_WP - wz1
         vel(1) = wz1*(wy1*(wx1*pU(ip+1,jp+1,kp+1,1)+wx2*pU(ip,jp+1,kp+1,1))+wy2*(wx1*pU(ip+1,jp,kp+1,1)+wx2*pU(ip,jp,kp+1,1))) + &
         &        wz2*(wy1*(wx1*pU(ip+1,jp+1,kp  ,1)+wx2*pU(ip,jp+1,kp  ,1))+wy2*(wx1*pU(ip+1,jp,kp  ,1)+wx2*pU(ip,jp,kp  ,1)))
         vel(2) = wz1*(wy1*(wx1*pV(ip+1,jp+1,kp+1,1)+wx2*pV(ip,jp+1,kp+1,1))+wy2*(wx1*pV(ip+1,jp,kp+1,1)+wx2*pV(ip,jp,kp+1,1))) + &
         &        wz2*(wy1*(wx1*pV(ip+1,jp+1,kp  ,1)+wx2*pV(ip,jp+1,kp  ,1))+wy2*(wx1*pV(ip+1,jp,kp  ,1)+wx2*pV(ip,jp,kp  ,1)))
         vel(3) = wz1*(wy1*(wx1*pW(ip+1,jp+1,kp+1,1)+wx2*pW(ip,jp+1,kp+1,1))+wy2*(wx1*pW(ip+1,jp,kp+1,1)+wx2*pW(ip,jp,kp+1,1))) + &
         &        wz2*(wy1*(wx1*pW(ip+1,jp+1,kp  ,1)+wx2*pW(ip,jp+1,kp  ,1))+wy2*(wx1*pW(ip+1,jp,kp  ,1)+wx2*pW(ip,jp,kp  ,1)))
      end function interp_velocity_col
      
      subroutine volume_correct_x_col(volume)
         real(WP), intent(inout) :: volume
         real(WP), dimension(3) :: va, vb, vc
         integer :: nt
         do nt = 1, 6
            va = face_col(:,tet_map(1,nt)) - face_col(:,tet_map(4,nt))
            vb = face_col(:,tet_map(2,nt)) - face_col(:,tet_map(4,nt))
            vc = face_col(:,tet_map(3,nt)) - face_col(:,tet_map(4,nt))
            volume = volume + (va(1)*(vb(2)*vc(3)-vc(2)*vb(3)) - va(2)*(vb(1)*vc(3)-vc(1)*vb(3)) + va(3)*(vb(1)*vc(2)-vc(1)*vb(2))) / 6.0_WP
         end do
         face_col(1,9) = (-6.0_WP*volume + face_col(1,5)*((face_col(2,8)-face_col(2,9))*(face_col(3,6)-face_col(3,9))-(face_col(2,6)-face_col(2,9))*(face_col(3,8)-face_col(3,9))) + &
         &   face_col(2,5)*((face_col(3,8)-face_col(3,9))*face_col(1,6)-(face_col(3,6)-face_col(3,9))*face_col(1,8)) + &
         &   face_col(2,9)*((face_col(3,6)-face_col(3,9))*face_col(1,8)-(face_col(3,8)-face_col(3,9))*face_col(1,6)) + &
         &   face_col(3,5)*((face_col(2,6)-face_col(2,9))*face_col(1,8)-(face_col(2,8)-face_col(2,9))*face_col(1,6)) + &
         &   face_col(3,9)*((face_col(2,8)-face_col(2,9))*face_col(1,6)-(face_col(2,6)-face_col(2,9))*face_col(1,8)) + &
         &   face_col(1,6)*((face_col(2,8)-face_col(2,9))*(face_col(3,7)-face_col(3,9))-(face_col(2,7)-face_col(2,9))*(face_col(3,8)-face_col(3,9))) + &
         &   face_col(2,6)*((face_col(3,8)-face_col(3,9))*face_col(1,7)-(face_col(3,7)-face_col(3,9))*face_col(1,8)) + &
         &   face_col(2,9)*((face_col(3,7)-face_col(3,9))*face_col(1,8)-(face_col(3,8)-face_col(3,9))*face_col(1,7)) + &
         &   face_col(3,6)*((face_col(2,7)-face_col(2,9))*face_col(1,8)-(face_col(2,8)-face_col(2,9))*face_col(1,7)) + &
         &   face_col(3,9)*((face_col(2,8)-face_col(2,9))*face_col(1,7)-(face_col(2,7)-face_col(2,9))*face_col(1,8))) / &
         &  (-(face_col(2,6)-face_col(2,9))*(face_col(3,8)-face_col(3,9))+(face_col(2,8)-face_col(2,9))*(face_col(3,6)-face_col(3,9)) - &
         &   face_col(2,5)*(face_col(3,6)-face_col(3,9))+face_col(2,5)*(face_col(3,8)-face_col(3,9))+face_col(2,9)*(face_col(3,6)-face_col(3,9))-face_col(2,9)*(face_col(3,8)-face_col(3,9)) - &
         &   face_col(3,5)*(face_col(2,8)-face_col(2,9))+face_col(3,5)*(face_col(2,6)-face_col(2,9))+face_col(3,9)*(face_col(2,8)-face_col(2,9))-face_col(3,9)*(face_col(2,6)-face_col(2,9)) - &
         &   (face_col(2,7)-face_col(2,9))*(face_col(3,8)-face_col(3,9))+(face_col(2,8)-face_col(2,9))*(face_col(3,7)-face_col(3,9)) - &
         &   face_col(2,6)*(face_col(3,7)-face_col(3,9))+face_col(2,6)*(face_col(3,8)-face_col(3,9))+face_col(2,9)*(face_col(3,7)-face_col(3,9))-face_col(2,9)*(face_col(3,8)-face_col(3,9)) - &
         &   face_col(3,6)*(face_col(2,8)-face_col(2,9))+face_col(3,6)*(face_col(2,7)-face_col(2,9))+face_col(3,9)*(face_col(2,8)-face_col(2,9))-face_col(3,9)*(face_col(2,7)-face_col(2,9)))
         face_col(2,9) = 0.25_WP * sum(face_col(2,5:8)); face_col(3,9) = 0.25_WP * sum(face_col(3,5:8))
      end subroutine volume_correct_x_col
      
      subroutine volume_correct_y_col(volume)
         real(WP), intent(inout) :: volume
         real(WP), dimension(3) :: va, vb, vc
         integer :: nt
         do nt = 1, 6
            va = face_col(:,tet_map(1,nt)) - face_col(:,tet_map(4,nt))
            vb = face_col(:,tet_map(2,nt)) - face_col(:,tet_map(4,nt))
            vc = face_col(:,tet_map(3,nt)) - face_col(:,tet_map(4,nt))
            volume = volume - (va(1)*(vb(2)*vc(3)-vc(2)*vb(3)) - va(2)*(vb(1)*vc(3)-vc(1)*vb(3)) + va(3)*(vb(1)*vc(2)-vc(1)*vb(2))) / 6.0_WP
         end do
         face_col(1,9) = 0.25_WP * sum(face_col(1,5:8))
         face_col(2,9) = (6.0_WP*volume + face_col(1,5)*((face_col(3,6)-face_col(3,9))*face_col(2,8)-(face_col(3,8)-face_col(3,9))*face_col(2,6)) + &
         &   face_col(1,9)*((face_col(3,8)-face_col(3,9))*face_col(2,6)-(face_col(3,6)-face_col(3,9))*face_col(2,8)) + &
         &   face_col(2,5)*((face_col(3,8)-face_col(3,9))*(face_col(1,6)-face_col(1,9))-(face_col(3,6)-face_col(3,9))*(face_col(1,8)-face_col(1,9))) + &
         &   face_col(3,5)*((face_col(1,8)-face_col(1,9))*face_col(2,6)-(face_col(1,6)-face_col(1,9))*face_col(2,8)) + &
         &   face_col(3,9)*((face_col(1,6)-face_col(1,9))*face_col(2,8)-(face_col(1,8)-face_col(1,9))*face_col(2,6)) + &
         &   face_col(1,6)*((face_col(3,7)-face_col(3,9))*face_col(2,8)-(face_col(3,8)-face_col(3,9))*face_col(2,7)) + &
         &   face_col(1,9)*((face_col(3,8)-face_col(3,9))*face_col(2,7)-(face_col(3,7)-face_col(3,9))*face_col(2,8)) + &
         &   face_col(2,6)*((face_col(3,8)-face_col(3,9))*(face_col(1,7)-face_col(1,9))-(face_col(3,7)-face_col(3,9))*(face_col(1,8)-face_col(1,9))) + &
         &   face_col(3,6)*((face_col(1,8)-face_col(1,9))*face_col(2,7)-(face_col(1,7)-face_col(1,9))*face_col(2,8)) + &
         &   face_col(3,9)*((face_col(1,7)-face_col(1,9))*face_col(2,8)-(face_col(1,8)-face_col(1,9))*face_col(2,7))) / &
         &  (face_col(1,5)*((face_col(3,6)-face_col(3,9))-(face_col(3,8)-face_col(3,9))) + &
         &   face_col(1,9)*((face_col(3,8)-face_col(3,9))-(face_col(3,6)-face_col(3,9))) + &
         &   ((face_col(3,8)-face_col(3,9))*(face_col(1,6)-face_col(1,9))-(face_col(3,6)-face_col(3,9))*(face_col(1,8)-face_col(1,9))) + &
         &   face_col(3,5)*((face_col(1,8)-face_col(1,9))-(face_col(1,6)-face_col(1,9)))+face_col(3,9)*((face_col(1,6)-face_col(1,9))-(face_col(1,8)-face_col(1,9))) + &
         &   face_col(1,9)*((face_col(3,8)-face_col(3,9))-(face_col(3,7)-face_col(3,9))) + &
         &   ((face_col(3,8)-face_col(3,9))*(face_col(1,7)-face_col(1,9))-(face_col(3,7)-face_col(3,9))*(face_col(1,8)-face_col(1,9))) + &
         &   face_col(3,6)*((face_col(1,8)-face_col(1,9))-(face_col(1,7)-face_col(1,9)))+face_col(3,9)*((face_col(1,7)-face_col(1,9))-(face_col(1,8)-face_col(1,9))))
         face_col(3,9) = 0.25_WP * sum(face_col(3,5:8))
      end subroutine volume_correct_y_col
      
      subroutine volume_correct_z_col(volume)
         real(WP), intent(inout) :: volume
         real(WP), dimension(3) :: va, vb, vc
         integer :: nt
         do nt = 1, 6
            va = face_col(:,tet_map(1,nt)) - face_col(:,tet_map(4,nt))
            vb = face_col(:,tet_map(2,nt)) - face_col(:,tet_map(4,nt))
            vc = face_col(:,tet_map(3,nt)) - face_col(:,tet_map(4,nt))
            volume = volume + (va(1)*(vb(2)*vc(3)-vc(2)*vb(3)) - va(2)*(vb(1)*vc(3)-vc(1)*vb(3)) + va(3)*(vb(1)*vc(2)-vc(1)*vb(2))) / 6.0_WP
         end do
         face_col(1,9) = 0.25_WP * sum(face_col(1,5:8)); face_col(2,9) = 0.25_WP * sum(face_col(2,5:8))
         face_col(3,9) = (6.0_WP*volume + &
         &   face_col(1,5)*face_col(2,6)*face_col(3,8) - face_col(1,5)*face_col(3,6)*face_col(2,8) - &
         &   face_col(2,5)*face_col(1,6)*face_col(3,8) + face_col(2,5)*face_col(3,6)*face_col(1,8) + &
         &   face_col(3,5)*face_col(1,6)*face_col(2,8) - face_col(3,5)*face_col(2,6)*face_col(1,8) + &
         &   face_col(1,5)*face_col(3,6)*face_col(2,5) - face_col(2,5)*face_col(3,6)*face_col(1,5) - &
         &   face_col(3,5)*face_col(1,6)*face_col(2,5) + face_col(3,5)*face_col(2,6)*face_col(1,5) + &
         &   face_col(1,6)*face_col(2,7)*face_col(3,8) - face_col(1,6)*face_col(3,7)*face_col(2,8) - &
         &   face_col(2,6)*face_col(1,7)*face_col(3,8) + face_col(2,6)*face_col(3,7)*face_col(1,8) + &
         &   face_col(3,6)*face_col(1,7)*face_col(2,8) - face_col(3,6)*face_col(2,7)*face_col(1,8) - &
         &   face_col(1,5)*face_col(3,8)*face_col(2,5) + face_col(2,5)*face_col(3,8)*face_col(1,5) + &
         &   face_col(3,5)*face_col(1,8)*face_col(2,5) - face_col(3,5)*face_col(2,8)*face_col(1,5) + &
         &   face_col(1,6)*face_col(3,7)*face_col(2,5) - face_col(2,6)*face_col(3,7)*face_col(1,5) - &
         &   face_col(3,6)*face_col(1,7)*face_col(2,5) + face_col(3,6)*face_col(2,7)*face_col(1,5) + &
         &   face_col(1,7)*face_col(3,8)*face_col(2,5) - face_col(2,7)*face_col(3,8)*face_col(1,5) - &
         &   face_col(3,7)*face_col(1,8)*face_col(2,5) + face_col(3,7)*face_col(2,8)*face_col(1,5)) / &
         &  (face_col(1,5)*face_col(2,6) - face_col(2,5)*face_col(1,6) - face_col(1,5)*face_col(2,8) + face_col(2,5)*face_col(1,8) + &
         &   face_col(1,6)*face_col(2,7) - face_col(2,6)*face_col(1,7) + face_col(1,7)*face_col(2,8) - face_col(2,7)*face_col(1,8))
      end subroutine volume_correct_z_col
      
      
   end subroutine advance_vof_col

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
         this%VFmin = min(this%VFmin, this%VF%get_min(lvl=lvl))
         this%VFmax = max(this%VFmax, this%VF%get_max(lvl=lvl))
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
      call io%add_data(this%VF, 'VF')
      call io%add_data(this%Cliq, 'Cliq')
      call io%add_data(this%Cgas, 'Cgas')
      call io%add_data(this%PLIC, 'PLIC')
   end subroutine register_checkpoint

   !> Restore checkpoint
   subroutine restore_checkpoint(this, io, dirname)
      use amrio_class, only: amrio
      class(amrvof), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      call io%read_data(dirname, this%VF, 'VF')
      call io%read_data(dirname, this%Cliq, 'Cliq')
      call io%read_data(dirname, this%Cgas, 'Cgas')
      call io%read_data(dirname, this%PLIC, 'PLIC')
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
