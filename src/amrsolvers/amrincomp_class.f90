!> AMR Incompressible solver class
!> Provides projection and basic operations for constant-density flow
module amrincomp_class
   use iso_c_binding,    only: c_ptr, c_null_ptr, c_loc, c_f_pointer
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use amrmg_class,      only: amrmg, amrmg_cstcoef
   use amrex_amr_module, only: amrex_multifab, amrex_mfiter, amrex_box, &
   &                           amrex_boxarray, amrex_distromap, amrex_geometry, &
   &                           amrex_interp_face_divfree
   implicit none
   private

   ! Expose type and dispatchers
   public :: amrincomp
   public :: amrincomp_on_init, amrincomp_on_coarse, amrincomp_on_remake
   public :: amrincomp_on_clear, amrincomp_tagging, amrincomp_postregrid

   !> AMR Incompressible solver type
   type, extends(amrsolver) :: amrincomp
      ! User-configurable callbacks
      procedure(incomp_init_iface), pointer, nopass :: user_init => null()
      procedure(incomp_tagging_iface), pointer, nopass :: user_tagging => null()

      ! Flow data (solver owns these)
      type(amrdata) :: U, V, W          !< Current face velocities
      type(amrdata) :: Uold, Vold, Wold !< Old face velocities
      type(amrdata) :: P                !< Pressure (cell-centered)
      type(amrdata) :: div              !< Divergence (cell-centered)

      ! Pressure solver
      type(amrmg) :: psolver

      ! Physical properties
      real(WP) :: rho  = 1.0_WP   !< Density (constant)
      real(WP) :: visc = 0.0_WP   !< Dynamic viscosity

      ! Monitoring quantities
      real(WP) :: divmax = 0.0_WP !< Maximum divergence

      ! Velocity interpolation method (for C/F interface fills)
      integer :: interp = amrex_interp_face_divfree

   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: get_div      !< Compute divergence (assumes velocity ghosts filled)
      procedure :: get_pgrad    !< Compute pressure gradient (assumes pressure ghosts filled)
      ! Override internal type-bound callbacks from amrsolver
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      procedure :: fill_velocity_lvl         !< Fill velocity ghosts at single level
      procedure :: fill_velocity             !< Fill velocity ghosts on all levels
      procedure :: fill_velocity_from_coarse !< Fill velocity from coarse
      procedure :: sync_velocity_lvl         !< Sync velocity ghosts at single level
      procedure :: sync_velocity             !< Sync velocity ghosts on all levels
      procedure :: fill_velocity_mfab        !< Fill dest MultiFabs for regridding
      procedure :: average_down_velocity     !< Average down MAC velocity for C/F consistency
      procedure :: average_down_velocity_to  !< Average down MAC velocity for single level
      procedure :: average_down_pressure     !< Average down pressure for C/F consistency
      procedure :: average_down_pressure_to  !< Average down pressure for single level
      ! Deferred from amrsolver base class
      procedure :: get_info
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
      ! Physics procedures
      procedure :: get_dmomdt                 !< Compute momentum advection RHS
   end type amrincomp

   !> Abstract interface for user-overridable on_init callback
   abstract interface
      subroutine incomp_init_iface(solver, lvl, time, ba, dm)
         import :: amrincomp, WP, amrex_boxarray, amrex_distromap
         class(amrincomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine incomp_init_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine incomp_tagging_iface(solver, lvl, tags, time)
         import :: amrincomp, c_ptr, WP
         class(amrincomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         type(c_ptr), intent(in) :: tags
         real(WP), intent(in) :: time
      end subroutine incomp_tagging_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrincomp type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrincomp_on_init(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_init(lvl, time, ba, dm)
      if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
   end subroutine amrincomp_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrincomp_on_coarse(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_coarse(lvl, time, ba, dm)
   end subroutine amrincomp_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrincomp_on_remake(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_remake(lvl, time, ba, dm)
   end subroutine amrincomp_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrincomp_on_clear(ctx, lvl)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_clear(lvl)
   end subroutine amrincomp_on_clear

   !> Dispatch tagging: calls user callback if set
   subroutine amrincomp_tagging(ctx, lvl, tags, time)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      if (associated(this%user_tagging)) call this%user_tagging(this, lvl, tags, time)
   end subroutine amrincomp_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrincomp_postregrid(ctx, lbase, time)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx, this)
      call this%post_regrid(lbase, time)
   end subroutine amrincomp_postregrid

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the incompressible solver
   subroutine initialize(this, amr, name)
      class(amrincomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Set name
      if (present(name)) then
         this%name = trim(name)
      else
         this%name = 'UNNAMED_INCOMP'
      end if

      ! Store amrgrid pointer
      this%amr => amr

      ! Initialize staggered velocity (do we need ghosts?)
      call this%U%initialize(amr, name='U', ncomp=1, ng=1, nodal=[.true., .false., .false.])
      call this%V%initialize(amr, name='V', ncomp=1, ng=1, nodal=[.false., .true., .false.])
      call this%W%initialize(amr, name='W', ncomp=1, ng=1, nodal=[.false., .false., .true.])
      call this%Uold%initialize(amr, name='Uold', ncomp=1, ng=1, nodal=[.true., .false., .false.])
      call this%Vold%initialize(amr, name='Vold', ncomp=1, ng=1, nodal=[.false., .true., .false.])
      call this%Wold%initialize(amr, name='Wold', ncomp=1, ng=1, nodal=[.false., .false., .true.])

      ! Initialize pressure (cell-centered) (do we need ghosts?)
      call this%P%initialize(amr, name='P', ncomp=1, ng=1)

      ! Initialize divergence (cell-centered, no ghosts needed)
      call this%div%initialize(amr, name='div', ncomp=1, ng=0)

      ! Set parent pointers for BC context access
      this%U%parent => this; this%V%parent => this; this%W%parent => this
      this%Uold%parent => this; this%Vold%parent => this; this%Wold%parent => this
      this%P%parent => this
      this%div%parent => this

      ! Set velocity fillbc callbacks to shared handler
      this%U%fillbc => velocity_fillbc
      this%V%fillbc => velocity_fillbc
      this%W%fillbc => velocity_fillbc

      ! Initialize pressure solver
      call this%psolver%initialize(amr, type=amrmg_cstcoef)

      ! Register all 6 callbacks with amrgrid using concrete dispatchers
      select type (this)
       type is (amrincomp)
         call this%amr%add_on_init   (amrincomp_on_init,    c_loc(this))
         call this%amr%add_on_coarse (amrincomp_on_coarse,  c_loc(this))
         call this%amr%add_on_remake (amrincomp_on_remake,  c_loc(this))
         call this%amr%add_on_clear  (amrincomp_on_clear,   c_loc(this))
         call this%amr%add_tagging   (amrincomp_tagging,    c_loc(this))
         call this%amr%add_postregrid(amrincomp_postregrid, c_loc(this))
      end select

   end subroutine initialize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this, lvl, time, ba, dm)
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level layouts
      call this%U%reset_level(lvl, ba, dm)
      call this%V%reset_level(lvl, ba, dm)
      call this%W%reset_level(lvl, ba, dm)
      call this%Uold%reset_level(lvl, ba, dm)
      call this%Vold%reset_level(lvl, ba, dm)
      call this%Wold%reset_level(lvl, ba, dm)
      call this%P%reset_level(lvl, ba, dm)
      call this%div%reset_level(lvl, ba, dm)
      ! Set to zero
      call this%U%setval(val=0.0_WP, lvl=lvl)
      call this%V%setval(val=0.0_WP, lvl=lvl)
      call this%W%setval(val=0.0_WP, lvl=lvl)
      call this%Uold%setval(val=0.0_WP, lvl=lvl)
      call this%Vold%setval(val=0.0_WP, lvl=lvl)
      call this%Wold%setval(val=0.0_WP, lvl=lvl)
      call this%P%setval(val=0.0_WP, lvl=lvl)
      call this%div%setval(val=0.0_WP, lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using divergence-free interpolation
   subroutine on_coarse(this, lvl, time, ba, dm)
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Velocity: allocate then fill with divergence-free interpolation
      call this%U%reset_level(lvl, ba, dm)
      call this%V%reset_level(lvl, ba, dm)
      call this%W%reset_level(lvl, ba, dm)
      call this%fill_velocity_from_coarse(lvl, time)
      ! Old velocity just needs geometry
      call this%Uold%reset_level(lvl, ba, dm)
      call this%Vold%reset_level(lvl, ba, dm)
      call this%Wold%reset_level(lvl, ba, dm)
      ! Pressure on coarse
      call this%P%on_coarse(this%P, lvl, time, ba, dm)
      ! Divergence just needs geometry
      call this%div%reset_level(lvl, ba, dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using divergence-free interpolation
   subroutine on_remake(this, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_multifab_build, amrex_multifab_destroy
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_multifab) :: Utmp, Vtmp, Wtmp
      ! Build temp MultiFabs with new layout (0 ghost cells for FillPatch)
      call amrex_multifab_build(Utmp, ba, dm, 1, 0, this%U%nodal)
      call amrex_multifab_build(Vtmp, ba, dm, 1, 0, this%V%nodal)
      call amrex_multifab_build(Wtmp, ba, dm, 1, 0, this%W%nodal)
      ! Fill temps from old data via coupled FillPatch
      call this%fill_velocity_mfab(Utmp, Vtmp, Wtmp, lvl, time)
      ! Reset levels and copy from temps
      call this%U%reset_level(lvl, ba, dm)
      call this%V%reset_level(lvl, ba, dm)
      call this%W%reset_level(lvl, ba, dm)
      call this%U%mf(lvl)%copy(Utmp, 1, 1, 1, 0)
      call this%V%mf(lvl)%copy(Vtmp, 1, 1, 1, 0)
      call this%W%mf(lvl)%copy(Wtmp, 1, 1, 1, 0)
      ! Destroy temps
      call amrex_multifab_destroy(Utmp)
      call amrex_multifab_destroy(Vtmp)
      call amrex_multifab_destroy(Wtmp)
      ! Old velocity just needs geometry
      call this%Uold%reset_level(lvl, ba, dm)
      call this%Vold%reset_level(lvl, ba, dm)
      call this%Wold%reset_level(lvl, ba, dm)
      ! Pressure remake
      call this%P%on_remake(this%P, lvl, time, ba, dm)
      ! Divergence just needs geometry
      call this%div%reset_level(lvl, ba, dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this, lvl)
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%clear_level(lvl)
      call this%V%clear_level(lvl)
      call this%W%clear_level(lvl)
      call this%Uold%clear_level(lvl)
      call this%Vold%clear_level(lvl)
      call this%Wold%clear_level(lvl)
      call this%P%clear_level(lvl)
      call this%div%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this, lbase, time)
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      call this%average_down_velocity(lbase)
      call this%average_down_pressure(lbase)
   end subroutine post_regrid

   !> Average down MAC velocity for a single level (lvl+1 -> lvl)
   !> Uses amrdata infrastructure which handles face-centered averaging correctly
   subroutine average_down_velocity_to(this, lvl)
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%average_downto(lvl)
      call this%V%average_downto(lvl)
      call this%W%average_downto(lvl)
   end subroutine average_down_velocity_to

   !> Average down MAC velocity from finest to lbase
   !> Simply calls average_down_velocity_to in a loop
   !> @param lbase Optional: lowest level to average down to (default 0)
   subroutine average_down_velocity(this, lbase)
      class(amrincomp), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl, lb
      lb = 0; if (present(lbase)) lb = lbase
      do lvl = this%amr%clvl()-1, lb, -1
         call this%average_down_velocity_to(lvl)
      end do
   end subroutine average_down_velocity

   !> Fill velocity ghost cells at a single level using divergence-free interpolation
   subroutine fill_velocity_lvl(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two_faces
      use amrdata_class, only: amrdata_fillbc
      class(amrincomp), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: ctx_u, ctx_v, ctx_w
      type(c_funptr) :: bc_dispatch
      integer :: rr, lo_bc(9), hi_bc(9)
      real(WP) :: t_old, t_new

      ! Get contexts for each velocity component
      ctx_u = c_loc(this%U); ctx_v = c_loc(this%V); ctx_w = c_loc(this%W)
      bc_dispatch = c_funloc(amrdata_fillbc)
      t_old = time - 1.0e200_WP
      t_new = time

      if (lvl .eq. 0) then
         ! Level 0: single-level fill (just physical BCs)
         call amrmfab_fillpatch_single(this%U%mf(0), time, this%U%mf(0), &
         &   time, this%U%mf(0), this%amr%geom(0), ctx_u, bc_dispatch, time, 1, 1, 1)
         call amrmfab_fillpatch_single(this%V%mf(0), time, this%V%mf(0), &
         &   time, this%V%mf(0), this%amr%geom(0), ctx_v, bc_dispatch, time, 1, 1, 1)
         call amrmfab_fillpatch_single(this%W%mf(0), time, this%W%mf(0), &
         &   time, this%W%mf(0), this%amr%geom(0), ctx_w, bc_dispatch, time, 1, 1, 1)
      else
         ! Build combined BC array: [U_x,U_y,U_z, V_x,V_y,V_z, W_x,W_y,W_z]
         lo_bc(1:3) = this%U%lo_bc(:,1)
         lo_bc(4:6) = this%V%lo_bc(:,1)
         lo_bc(7:9) = this%W%lo_bc(:,1)
         hi_bc(1:3) = this%U%hi_bc(:,1)
         hi_bc(4:6) = this%V%hi_bc(:,1)
         hi_bc(7:9) = this%W%hi_bc(:,1)
         rr = this%amr%rref(lvl-1)
         ! Call 3-component divfree FillPatch from two levels
         call amrmfab_fillpatch_two_faces( &
         &   this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), time, &
         &   t_old, this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
         &   t_new, this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
         &   this%amr%geom(lvl-1), &
         &   t_old, this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
         &   t_new, this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
         &   this%amr%geom(lvl), &
         &   ctx_u, ctx_v, ctx_w, bc_dispatch, bc_dispatch, bc_dispatch, &
         &   1, 1, 1, rr, this%interp, lo_bc, hi_bc)
      end if
   end subroutine fill_velocity_lvl

   !> Fill velocity ghost cells on all levels
   subroutine fill_velocity(this, time)
      class(amrincomp), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%fill_velocity_lvl(lvl, time)
      end do
   end subroutine fill_velocity

   !> Fill new fine level velocity from coarse using divergence-free interpolation
   !> Used during MakeNewLevelFromCoarse (creation of new fine levels)
   subroutine fill_velocity_from_coarse(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillcoarsepatch_faces
      use amrdata_class, only: amrdata_fillbc
      class(amrincomp), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: ctx_u, ctx_v, ctx_w
      type(c_funptr) :: bc_dispatch
      integer :: rr, lo_bc(9), hi_bc(9)
      ! Get contexts for each velocity component
      ctx_u = c_loc(this%U); ctx_v = c_loc(this%V); ctx_w = c_loc(this%W)
      bc_dispatch = c_funloc(amrdata_fillbc)
      ! Build combined BC array: [U_x,U_y,U_z, V_x,V_y,V_z, W_x,W_y,W_z]
      lo_bc(1:3) = this%U%lo_bc(:,1)
      lo_bc(4:6) = this%V%lo_bc(:,1)
      lo_bc(7:9) = this%W%lo_bc(:,1)
      hi_bc(1:3) = this%U%hi_bc(:,1)
      hi_bc(4:6) = this%V%hi_bc(:,1)
      hi_bc(7:9) = this%W%hi_bc(:,1)
      rr = this%amr%rref(lvl-1)
      ! Call 3-component divfree FillCoarsePatch
      call amrmfab_fillcoarsepatch_faces( &
      &   this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), time, &
      &   this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
      &   this%amr%geom(lvl-1), this%amr%geom(lvl), &
      &   ctx_u, ctx_v, ctx_w, bc_dispatch, bc_dispatch, bc_dispatch, &
      &   1, 1, 1, rr, this%interp, lo_bc, hi_bc)
   end subroutine fill_velocity_from_coarse

   !> Sync velocity ghost cells at a single level (lightweight, no C/F interpolation)
   subroutine sync_velocity_lvl(this, lvl)
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%sync_lvl(lvl)
      call this%V%sync_lvl(lvl)
      call this%W%sync_lvl(lvl)
   end subroutine sync_velocity_lvl

   !> Sync velocity ghost cells on all levels
   subroutine sync_velocity(this)
      class(amrincomp), intent(inout) :: this
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%sync_velocity_lvl(lvl)
      end do
   end subroutine sync_velocity

   !> Fill destination MultiFabs with velocity using divergence-free interpolation
   !> Used during regridding (on_remake) to fill new layout MultiFabs
   subroutine fill_velocity_mfab(this, Udest, Vdest, Wdest, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two_faces
      use amrdata_class, only: amrdata_fillbc
      class(amrincomp), target, intent(inout) :: this
      type(amrex_multifab), intent(inout) :: Udest, Vdest, Wdest
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: ctx_u, ctx_v, ctx_w
      type(c_funptr) :: bc_dispatch
      integer :: rr, lo_bc(9), hi_bc(9)
      real(WP) :: t_old, t_new

      ! Get contexts for each velocity component
      ctx_u = c_loc(this%U); ctx_v = c_loc(this%V); ctx_w = c_loc(this%W)
      bc_dispatch = c_funloc(amrdata_fillbc)
      t_old = time - 1.0e200_WP
      t_new = time

      if (lvl .eq. 0) then
         ! Level 0: single-level fill (just physical BCs)
         call amrmfab_fillpatch_single(Udest, t_old, this%U%mf(0), &
         &   t_new, this%U%mf(0), this%amr%geom(0), ctx_u, bc_dispatch, time, 1, 1, 1)
         call amrmfab_fillpatch_single(Vdest, t_old, this%V%mf(0), &
         &   t_new, this%V%mf(0), this%amr%geom(0), ctx_v, bc_dispatch, time, 1, 1, 1)
         call amrmfab_fillpatch_single(Wdest, t_old, this%W%mf(0), &
         &   t_new, this%W%mf(0), this%amr%geom(0), ctx_w, bc_dispatch, time, 1, 1, 1)
      else
         ! Build combined BC array
         lo_bc(1:3) = this%U%lo_bc(:,1)
         lo_bc(4:6) = this%V%lo_bc(:,1)
         lo_bc(7:9) = this%W%lo_bc(:,1)
         hi_bc(1:3) = this%U%hi_bc(:,1)
         hi_bc(4:6) = this%V%hi_bc(:,1)
         hi_bc(7:9) = this%W%hi_bc(:,1)
         rr = this%amr%rref(lvl-1)
         ! Call 3-component divfree FillPatch
         call amrmfab_fillpatch_two_faces( &
         &   Udest, Vdest, Wdest, time, &
         &   t_old, this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
         &   t_new, this%U%mf(lvl-1), this%V%mf(lvl-1), this%W%mf(lvl-1), &
         &   this%amr%geom(lvl-1), &
         &   t_old, this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
         &   t_new, this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
         &   this%amr%geom(lvl), &
         &   ctx_u, ctx_v, ctx_w, bc_dispatch, bc_dispatch, bc_dispatch, &
         &   1, 1, 1, rr, this%interp, lo_bc, hi_bc)
      end if
   end subroutine fill_velocity_mfab

   !> Average down pressure for a single level (lvl+1 -> lvl)
   !> Uses amrdata infrastructure which handles cell-centered averaging correctly
   subroutine average_down_pressure_to(this, lvl)
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%P%average_downto(lvl)
   end subroutine average_down_pressure_to

   !> Average down pressure from finest to lbase
   !> Simply calls average_down_pressure_to in a loop
   !> @param lbase Optional: lowest level to average down to (default 0)
   subroutine average_down_pressure(this, lbase)
      class(amrincomp), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl, lb
      lb = 0; if (present(lbase)) lb = lbase
      do lvl = this%amr%clvl()-1, lb, -1
         call this%average_down_pressure_to(lvl)
      end do
   end subroutine average_down_pressure

   !> Finalize the incompressible solver
   subroutine finalize(this)
      class(amrincomp), intent(inout) :: this
      call this%U%finalize()
      call this%V%finalize()
      call this%W%finalize()
      call this%Uold%finalize()
      call this%Vold%finalize()
      call this%Wold%finalize()
      call this%P%finalize()
      call this%div%finalize()
      call this%psolver%finalize()
      nullify(this%amr)
      nullify(this%user_init)
      nullify(this%user_tagging)
   end subroutine finalize

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Compute divergence of velocity into internal div field, update divmax
   subroutine get_div(this)
      use amrex_interface, only: amrmfab_compute_divergence
      class(amrincomp), intent(inout) :: this
      integer :: lvl
      ! Use our wrapper to amrex's 2nd order staggered divergence
      do lvl = 0, this%amr%clvl()
         call amrmfab_compute_divergence(this%div%mf(lvl), &
            this%U%mf(lvl), this%V%mf(lvl), this%W%mf(lvl), &
            this%amr%geom(lvl))
      end do
      ! Update divmax
      this%divmax = 0.0_WP
      do lvl = 0, this%amr%clvl()
         this%divmax = max(this%divmax, this%div%norm0(lvl=lvl))
      end do
   end subroutine get_div

   !> Compute pressure gradient into user-provided face amrdata
   subroutine get_pgrad(this, dPdx, dPdy, dPdz)
      class(amrincomp), intent(inout) :: this
      type(amrdata), intent(inout) :: dPdx, dPdy, dPdz
      integer :: lvl, i, j, k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP, pGx, pGy, pGz
      real(WP) :: dxi, dyi, dzi

      do lvl = 0, this%amr%clvl()
         dxi = 1.0_WP / this%amr%dx(lvl)
         dyi = 1.0_WP / this%amr%dy(lvl)
         dzi = 1.0_WP / this%amr%dz(lvl)

         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pP  => this%P%mf(lvl)%dataptr(mfi)
            pGx => dPdx%mf(lvl)%dataptr(mfi)
            pGy => dPdy%mf(lvl)%dataptr(mfi)
            pGz => dPdz%mf(lvl)%dataptr(mfi)

            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pGx(i,j,k,1) = (pP(i,j,k,1) - pP(i-1,j,k,1)) * dxi
                     pGy(i,j,k,1) = (pP(i,j,k,1) - pP(i,j-1,k,1)) * dyi
                     pGz(i,j,k,1) = (pP(i,j,k,1) - pP(i,j,k-1,1)) * dzi
                  end do
               end do
            end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_pgrad

   !> Get solver information (max velocity, divergence, etc.)
   subroutine get_info(this)
      class(amrincomp), intent(inout) :: this
      ! Compute divmax (call get_div to update it)
      call this%get_div()
      ! TODO: Compute Umax, Vmax, Wmax
   end subroutine get_info

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this, io)
      use amrio_class, only: amrio
      class(amrincomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%U, 'U')
      call io%add_data(this%V, 'V')
      call io%add_data(this%W, 'W')
      call io%add_data(this%P, 'P')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this, io, dirname)
      use amrio_class, only: amrio
      class(amrincomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      call io%read_data(dirname, this%U, 'U')
      call io%read_data(dirname, this%V, 'V')
      call io%read_data(dirname, this%W, 'W')
      call io%read_data(dirname, this%P, 'P')
   end subroutine restore_checkpoint

   !> Compute momentum advection RHS for all levels
   !> drhoUdt = -∇·(ρuu), drhoVdt = -∇·(ρuv), drhoWdt = -∇·(ρuw)
   !> Uses flux averaging at C/F interfaces for conservation
   subroutine get_dmomdt(this, U, V, W, drhoUdt, drhoVdt, drhoWdt)
      use amrex_amr_module, only: amrex_multifab, amrex_multifab_destroy, amrex_mfiter, amrex_box
      class(amrincomp), intent(inout) :: this
      class(amrdata), intent(inout) :: U, V, W                          !< Velocity state (face-centered)
      class(amrdata), intent(inout) :: drhoUdt, drhoVdt, drhoWdt        !< Output: momentum RHS (face-centered)
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: FUx, FUy, FUz
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: FVx, FVy, FVz
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: FWx, FWy, FWz
      type(amrex_multifab) :: Ufill, Vfill, Wfill
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: lvl

      ! Phase 1: Compute fluxes on all levels
      do lvl = 0, this%amr%clvl()
         ! FUx, FVy, FWz: cell-centered
         call this%amr%mfab_build(lvl, FUx(lvl), ncomp=1, nover=1, atface=[.false.,.false.,.false.])
         call this%amr%mfab_build(lvl, FVy(lvl), ncomp=1, nover=1, atface=[.false.,.false.,.false.])
         call this%amr%mfab_build(lvl, FWz(lvl), ncomp=1, nover=1, atface=[.false.,.false.,.false.])
         ! FUy, FVx: xy-edge
         call this%amr%mfab_build(lvl, FUy(lvl), ncomp=1, nover=1, atface=[.true. ,.true. ,.false.])
         call this%amr%mfab_build(lvl, FVx(lvl), ncomp=1, nover=1, atface=[.true. ,.true. ,.false.])
         ! FUz, FWx: xz-edge
         call this%amr%mfab_build(lvl, FUz(lvl), ncomp=1, nover=1, atface=[.true. ,.false.,.true. ])
         call this%amr%mfab_build(lvl, FWx(lvl), ncomp=1, nover=1, atface=[.true. ,.false.,.true. ])
         ! FVz, FWy: yz-edge
         call this%amr%mfab_build(lvl, FVz(lvl), ncomp=1, nover=1, atface=[.false.,.true. ,.true. ])
         call this%amr%mfab_build(lvl, FWy(lvl), ncomp=1, nover=1, atface=[.false.,.true. ,.true. ])

         ! Build velocity fills with ghost cells
         call this%amr%mfab_build(lvl, Ufill, ncomp=1, nover=1, atface=[.true. ,.false.,.false.])
         call this%amr%mfab_build(lvl, Vfill, ncomp=1, nover=1, atface=[.false.,.true. ,.false.])
         call this%amr%mfab_build(lvl, Wfill, ncomp=1, nover=1, atface=[.false.,.false.,.true. ])
         call U%fill_mfab(Ufill, lvl, 0.0_WP)
         call V%fill_mfab(Vfill, lvl, 0.0_WP)
         call W%fill_mfab(Wfill, lvl, 0.0_WP)

         ! MFIter loop: compute fluxes
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            ! TODO: 2nd-order centered flux computation
         end do
         call this%amr%mfiter_destroy(mfi)

         ! Destroy velocity fills
         call amrex_multifab_destroy(Ufill)
         call amrex_multifab_destroy(Vfill)
         call amrex_multifab_destroy(Wfill)
      end do

      ! Phase 2: Average down fluxes
      ! Phase 3: Divergence

      ! Cleanup
      do lvl = 0, this%amr%clvl()
         call amrex_multifab_destroy(FUx(lvl))
         call amrex_multifab_destroy(FUy(lvl))
         call amrex_multifab_destroy(FUz(lvl))
         call amrex_multifab_destroy(FVx(lvl))
         call amrex_multifab_destroy(FVy(lvl))
         call amrex_multifab_destroy(FVz(lvl))
         call amrex_multifab_destroy(FWx(lvl))
         call amrex_multifab_destroy(FWy(lvl))
         call amrex_multifab_destroy(FWz(lvl))
      end do

   end subroutine get_dmomdt

   !> Velocity boundary condition callback (shared by U, V, W)
   !> Identifies component by its name and applies appropriate BCs
   !> TODO: Study incflo BC implementation for proper face-centered velocity BCs
   subroutine velocity_fillbc(this, mf, scomp, ncomp, time, geom)
      use amrex_amr_module, only: amrex_filcc, amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp, ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer, dimension(4) :: plo, phi

      ! Access solver and check for all-periodic (skip if so)
      select type (solver => this%parent)
       class is (amrincomp)
         if (solver%amr%xper .and. solver%amr%yper .and. solver%amr%zper) return
      end select

      ! Loop over boxes and apply BCs
      ! Component identified by this%name: 'U', 'V', or 'W'
      call amrex_mfiter_build(mfi, mf, tiling=.false.)
      do while (mfi%next())
         p => mf%dataptr(mfi)
         if (.not.geom%domain%contains(p)) then
            plo = lbound(p); phi = ubound(p)
            ! Apply BCs using this amrdata's lo_bc/hi_bc
            ! Note: amrex_filcc is designed for cell-centered data
            ! TODO: Implement proper face-centered BC logic (see incflo)
            call amrex_filcc(p, plo, phi, geom%domain%lo, geom%domain%hi, geom%dx, &
            &   geom%get_physical_location(plo), this%lo_bc, this%hi_bc)
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine velocity_fillbc

end module amrincomp_class
