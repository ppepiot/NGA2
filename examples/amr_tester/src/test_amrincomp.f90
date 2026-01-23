!> Test amrincomp solver - advection test with time loop
module mod_test_amrincomp
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrincomp_class,   only: amrincomp
   use amrdata_class,     only: amrdata, amrex_interp_none
   implicit none
   private
   public :: test_amrincomp

   ! Module-level objects
   type(amrgrid), allocatable, target :: amr
   type(amrincomp), allocatable, target :: fs
   type(amrviz), allocatable :: viz
   type(amrdata), allocatable :: dPdx, dPdy, dPdz
   type(amrdata), allocatable :: drhoUdt, drhoVdt, drhoWdt

contains

   !> User-provided initialization callback for velocity - uniform flow with perturbation
   subroutine velocity_init(solver, lvl, time, ba, dm)
      use amrex_amr_module,  only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_distromap, amrex_mfiter_build, amrex_mfiter_destroy
      use mathtools, only: Pi
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      real(WP) :: x, y, z, dx, dy, dz, Lx, Ly, Lz
      integer :: i, j, k

      dx = solver%amr%dx(lvl); dy = solver%amr%dy(lvl); dz = solver%amr%dz(lvl)
      Lx = solver%amr%xhi - solver%amr%xlo
      Ly = solver%amr%yhi - solver%amr%ylo
      Lz = solver%amr%zhi - solver%amr%zlo

      ! Build mfiter from cell-centered ba/dm (passed during regrid)
      call amrex_mfiter_build(mfi, ba, dm, tiling=solver%amr%default_tiling)
      do while (mfi%next())
         pU => solver%U%mf(lvl)%dataptr(mfi)
         pV => solver%V%mf(lvl)%dataptr(mfi)
         pW => solver%W%mf(lvl)%dataptr(mfi)

         ! U (x-faces): uniform + Gaussian perturbation
         fbx = mfi%nodaltilebox(1)
         do k = fbx%lo(3), fbx%hi(3); do j = fbx%lo(2), fbx%hi(2); do i = fbx%lo(1), fbx%hi(1)
            x = solver%amr%xlo + real(i,WP) * dx
            y = solver%amr%ylo + (real(j,WP) + 0.5_WP) * dy
            z = solver%amr%zlo + (real(k,WP) + 0.5_WP) * dz
            pU(i,j,k,1) = 1.0_WP + 0.5_WP * exp(-((x-0.5_WP)**2 + (y-0.5_WP)**2 + (z-0.5_WP)**2) / 0.1_WP**2)
         end do; end do; end do

         ! V (y-faces): zero
         fbx = mfi%nodaltilebox(2)
         do k = fbx%lo(3), fbx%hi(3); do j = fbx%lo(2), fbx%hi(2); do i = fbx%lo(1), fbx%hi(1)
            pV(i,j,k,1) = 0.0_WP
         end do; end do; end do

         ! W (z-faces): zero
         fbx = mfi%nodaltilebox(3)
         do k = fbx%lo(3), fbx%hi(3); do j = fbx%lo(2), fbx%hi(2); do i = fbx%lo(1), fbx%hi(1)
            pW(i,j,k,1) = 0.0_WP
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine velocity_init

   !> Geometric tagger: refine center box
   subroutine geometric_tagger(solver, lvl, tags_ptr, time)
      use iso_c_binding,    only: c_ptr, c_char
      use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags_ptr
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP) :: xc, yc, zc, x, y, z, dx, dy, dz, radius
      integer :: i, j, k

      ! Center and radius for refinement region
      xc = 0.5_WP * (solver%amr%xlo + solver%amr%xhi)
      yc = 0.5_WP * (solver%amr%ylo + solver%amr%yhi)
      zc = 0.5_WP * (solver%amr%zlo + solver%amr%zhi)
      radius = 0.25_WP * min(solver%amr%xhi - solver%amr%xlo, &
      &                      solver%amr%yhi - solver%amr%ylo, &
      &                      solver%amr%zhi - solver%amr%zlo)

      dx = solver%amr%dx(lvl); dy = solver%amr%dy(lvl); dz = solver%amr%dz(lvl)
      tags = tags_ptr

      call solver%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
            x = solver%amr%xlo + (real(i,WP) + 0.5_WP) * dx
            y = solver%amr%ylo + (real(j,WP) + 0.5_WP) * dy
            z = solver%amr%zlo + (real(k,WP) + 0.5_WP) * dz
            if (sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2) .lt. radius) tagarr(i,j,k,1) = SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine geometric_tagger

   !> Main test routine
   subroutine test_amrincomp()
      use messager, only: log
      use string,   only: itoa, rtoa
      real(WP) :: time, dt, time_end, factor, CFL
      real(WP) :: Umax, divmax
      integer :: step

      call log("========================================")
      call log("TEST: amrincomp advection time loop")
      call log("========================================")

      time = 0.0_WP
      time_end = 1.0_WP
      CFL = 0.2_WP

      ! Create amrgrid
      allocate(amr)
      amr%name = 'advect_test'
      amr%nx = 32; amr%ny = 32; amr%nz = 32
      !amr%nx = 128; amr%ny = 128; amr%nz = 128
      amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
      amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
      amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
      amr%xper = .true.; amr%yper = .true.; amr%zper = .true.
      amr%maxlvl = 2  ! 3 levels: 0, 1, 2
      call amr%initialize()
      call log("Grid: "//trim(itoa(amr%nx))//"^3, maxlvl="//trim(itoa(amr%maxlvl)))

      ! Create flow solver
      allocate(fs)
      call fs%initialize(amr, name='advect_fs')
      fs%rho = 1.0_WP
      fs%user_init => velocity_init
      fs%user_tagging => geometric_tagger
      call log("Flow solver initialized")

      ! Create workspace for pressure gradients and momentum RHS
      allocate(dPdx); call dPdx%initialize(amr, name='dPdx', ncomp=1, ng=0, nodal=[.true., .false., .false.], interp=amrex_interp_none); call dPdx%register()
      allocate(dPdy); call dPdy%initialize(amr, name='dPdy', ncomp=1, ng=0, nodal=[.false., .true., .false.], interp=amrex_interp_none); call dPdy%register()
      allocate(dPdz); call dPdz%initialize(amr, name='dPdz', ncomp=1, ng=0, nodal=[.false., .false., .true.], interp=amrex_interp_none); call dPdz%register()
      allocate(drhoUdt); call drhoUdt%initialize(amr, name='drhoUdt', ncomp=1, ng=0, nodal=[.true., .false., .false.], interp=amrex_interp_none); call drhoUdt%register()
      allocate(drhoVdt); call drhoVdt%initialize(amr, name='drhoVdt', ncomp=1, ng=0, nodal=[.false., .true., .false.], interp=amrex_interp_none); call drhoVdt%register()
      allocate(drhoWdt); call drhoWdt%initialize(amr, name='drhoWdt', ncomp=1, ng=0, nodal=[.false., .false., .true.], interp=amrex_interp_none); call drhoWdt%register()

      ! Initialize visualization
      allocate(viz); call viz%initialize(amr, 'test_advect')
      call viz%add_scalar(fs%P, 1, 'pressure')
      call viz%add_scalar(fs%div, 1, 'divergence')
      call viz%add_scalar(fs%U, 1, 'U')
      call viz%add_scalar(fs%V, 1, 'V')
      call viz%add_scalar(fs%W, 1, 'W')

      ! Build grid and initialize velocity
      call amr%init_from_scratch(time=time)
      call fs%average_down_velocity()
      call fs%fill_velocity(time)
      call log("Grid built: "//trim(itoa(amr%nlevels))//" levels")

      ! Setup pressure solver
      fs%psolver%verbose = 0
      call fs%psolver%setup()

      ! Initial projection to make velocity divergence-free
      call fs%get_div()
      call log("Initial divmax = "//trim(rtoa(fs%divmax)))
      call project_velocity(dt_proj=1.0_WP)
      call fs%get_div()
      call log("After projection divmax = "//trim(rtoa(fs%divmax)))

      ! Write initial state
      call viz%write(time=time)

      ! ========================================
      ! TIME LOOP
      ! ========================================
      step = 0
      do while (time .lt. time_end)
         step = step + 1

         ! Compute CFL-based timestep
         Umax = max(fs%U%norm0(lvl=0), fs%V%norm0(lvl=0), fs%W%norm0(lvl=0))
         dt = CFL * amr%dx(amr%clvl()) / max(Umax, 1.0e-10_WP)
         if (time + dt .gt. time_end) dt = time_end - time

         ! Compute advective momentum RHS: drhoUdt = -div(rho*u*u)
         call fs%get_dmomdt(fs%U, fs%V, fs%W, drhoUdt, drhoVdt, drhoWdt)

         ! Update velocity: U_star = U + dt/rho * drhoUdt
         factor = dt / fs%rho
         call fs%U%saxpy(a=factor, src=drhoUdt)
         call fs%V%saxpy(a=factor, src=drhoVdt)
         call fs%W%saxpy(a=factor, src=drhoWdt)

         ! Project to enforce divergence-free
         call fs%average_down_velocity()
         call fs%fill_velocity(time)
         call project_velocity(dt_proj=dt)

         ! Advance time
         time = time + dt

         ! Log progress and write output
         if (mod(step, 10) .eq. 0 .or. time .ge. time_end) then
            call fs%get_div()
            call log("Step "//trim(itoa(step))//": t="//trim(rtoa(time))// &
            &     ", dt="//trim(rtoa(dt))//", divmax="//trim(rtoa(fs%divmax)))
            call viz%write(time=time)
         end if
      end do

      ! Write final state
      call viz%write(time=time)
      call log("Final time: "//trim(rtoa(time))//", steps: "//trim(itoa(step)))

      ! Cleanup
      call viz%finalize()
      call fs%finalize()
      call amr%finalize()
      deallocate(viz, fs, amr, dPdx, dPdy, dPdz, drhoUdt, drhoVdt, drhoWdt)

      call log("Advection test complete")

   contains

      !> Helper: project velocity to be divergence-free
      subroutine project_velocity(dt_proj)
         real(WP), intent(in) :: dt_proj
         real(WP) :: factor_proj
         ! Compute divergence
         call fs%get_div()
         ! Scale RHS: rhs = rho/dt * div
         factor_proj = fs%rho / dt_proj
         call fs%div%mult(val=factor_proj)
         ! Solve pressure Poisson
         call fs%psolver%solve(rhs=fs%div, phi=fs%P)
         ! Get gradients and correct velocity
         call fs%psolver%get_fluxes(fs%P, dPdx, dPdy, dPdz)
         factor_proj = dt_proj / fs%rho
         call fs%U%saxpy(a=factor_proj, src=dPdx)
         call fs%V%saxpy(a=factor_proj, src=dPdy)
         call fs%W%saxpy(a=factor_proj, src=dPdz)
         ! Average down and fill ghosts
         call fs%average_down_velocity()
         call fs%fill_velocity(time)
      end subroutine project_velocity

   end subroutine test_amrincomp

end module mod_test_amrincomp
