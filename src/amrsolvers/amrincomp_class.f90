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
   &                           amrex_boxarray, amrex_distromap
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
      procedure :: fill_velocity              !< Fill velocity ghosts on all levels
      procedure :: average_down_velocity     !< Average down MAC velocity for C/F consistency
      procedure :: average_down_velocity_to  !< Average down MAC velocity for single level
      procedure :: average_down_pressure     !< Average down pressure for C/F consistency
      procedure :: average_down_pressure_to  !< Average down pressure for single level
      ! Deferred from amrsolver base class
      procedure :: get_info
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
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

      ! Initialize velocity fields (face-centered)
      call this%U%initialize(amr, name='U', ncomp=1, ng=1, nodal=[.true., .false., .false.])
      call this%V%initialize(amr, name='V', ncomp=1, ng=1, nodal=[.false., .true., .false.])
      call this%W%initialize(amr, name='W', ncomp=1, ng=1, nodal=[.false., .false., .true.])
      call this%Uold%initialize(amr, name='Uold', ncomp=1, ng=1, nodal=[.true., .false., .false.])
      call this%Vold%initialize(amr, name='Vold', ncomp=1, ng=1, nodal=[.false., .true., .false.])
      call this%Wold%initialize(amr, name='Wold', ncomp=1, ng=1, nodal=[.false., .false., .true.])

      ! Initialize pressure (cell-centered)
      call this%P%initialize(amr, name='P', ncomp=1, ng=1)

      ! Initialize divergence (cell-centered, no ghosts needed)
      call this%div%initialize(amr, name='div', ncomp=1, ng=0)

      ! Set parent pointers for BC context access
      this%U%parent => this; this%V%parent => this; this%W%parent => this
      this%Uold%parent => this; this%Vold%parent => this; this%Wold%parent => this
      this%P%parent => this
      this%div%parent => this

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
      call this%U%mf(lvl)%setval(0.0_WP)
      call this%V%mf(lvl)%setval(0.0_WP)
      call this%W%mf(lvl)%setval(0.0_WP)
      call this%Uold%mf(lvl)%setval(0.0_WP)
      call this%Vold%mf(lvl)%setval(0.0_WP)
      call this%Wold%mf(lvl)%setval(0.0_WP)
      call this%P%mf(lvl)%setval(0.0_WP)
      call this%div%mf(lvl)%setval(0.0_WP)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse
   subroutine on_coarse(this, lvl, time, ba, dm)
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Velocity gets interpolated from coarse
      call this%U%on_coarse(this%U, lvl, time, ba, dm)
      call this%V%on_coarse(this%V, lvl, time, ba, dm)
      call this%W%on_coarse(this%W, lvl, time, ba, dm)
      ! Old velocity just needs geometry
      call this%Uold%reset_level(lvl, ba, dm)
      call this%Vold%reset_level(lvl, ba, dm)
      call this%Wold%reset_level(lvl, ba, dm)
      ! Pressure
      call this%P%on_coarse(this%P, lvl, time, ba, dm)
      ! Divergence just needs geometry
      call this%div%reset_level(lvl, ba, dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid
   subroutine on_remake(this, lvl, time, ba, dm)
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Delegate to amrdata on_remake
      call this%U%on_remake(this%U, lvl, time, ba, dm)
      call this%V%on_remake(this%V, lvl, time, ba, dm)
      call this%W%on_remake(this%W, lvl, time, ba, dm)
      call this%Uold%reset_level(lvl, ba, dm)
      call this%Vold%reset_level(lvl, ba, dm)
      call this%Wold%reset_level(lvl, ba, dm)
      call this%P%on_remake(this%P, lvl, time, ba, dm)
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

   !> Fill velocity ghost cells on all levels
   !> @param time Time for interpolation (used by FillPatch)
   subroutine fill_velocity(this, time)
      class(amrincomp), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%U%fill(lvl, time)
         call this%V%fill(lvl, time)
         call this%W%fill(lvl, time)
      end do
   end subroutine fill_velocity

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
         this%divmax = max(this%divmax, this%div%mf(lvl)%norm0(1))
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

end module amrincomp_class
