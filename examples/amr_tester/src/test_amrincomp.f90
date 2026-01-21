!> Test amrincomp solver - 3-level projection with random velocity and geometric tagger
module mod_test_amrincomp
   use precision,         only: WP
   use random,            only: random_uniform
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrincomp_class,   only: amrincomp
   use amrdata_class,     only: amrdata, amrex_interp_none
   use amrex_amr_module,  only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_distromap, &
   &                            amrex_mfiter_build, amrex_mfiter_destroy
   implicit none
   private
   public :: test_amrincomp

   ! Module-level objects
   type(amrgrid), allocatable, target :: amr
   type(amrincomp), allocatable, target :: fs
   type(amrviz), allocatable :: viz
   type(amrdata), allocatable :: dPdx, dPdy, dPdz  ! Workspace for pressure gradients

contains

   !> User-provided initialization callback for velocity - sinusoidal values
   subroutine velocity_init(solver, lvl, time, ba, dm)
      use mathtools, only: Pi
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx, fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      integer :: i, j, k

      ! For fine levels: fill from coarse (no direct computation)
      !if (lvl > 0) then
      !   call solver%U%fill_from_coarse(lvl, time)
      !   call solver%V%fill_from_coarse(lvl, time)
      !   call solver%W%fill_from_coarse(lvl, time)
      !   return
      !end if

      ! Build mfiter from cell-centered ba/dm (passed during regrid)
      call amrex_mfiter_build(mfi, ba, dm, tiling=solver%amr%default_tiling)
      do while (mfi%next())
         bx = mfi%tilebox()
         pU => solver%U%mf(lvl)%dataptr(mfi)
         pV => solver%V%mf(lvl)%dataptr(mfi)
         pW => solver%W%mf(lvl)%dataptr(mfi)

         ! U (x-faces): random values in [-1, 1]
         fbx = mfi%nodaltilebox(1)
         do k = fbx%lo(3), fbx%hi(3)
            do j = fbx%lo(2), fbx%hi(2)
               do i = fbx%lo(1), fbx%hi(1)
                  pU(i,j,k,1) = random_uniform(lo=-1.0_WP, hi=1.0_WP)
                  !pU(i,j,k,1) = sin(2.0_WP*Pi*real(i,WP)*solver%amr%dx(lvl))
                  !pU(i,j,k,1) = sin(2.0_WP*Pi*real(i,WP)*solver%amr%dx(lvl)) * cos(2.0_WP*Pi*(real(j,WP)+0.5_WP)*solver%amr%dy(lvl))
               end do
            end do
         end do

         ! V (y-faces): random values in [-1, 1]
         fbx = mfi%nodaltilebox(2)
         do k = fbx%lo(3), fbx%hi(3)
            do j = fbx%lo(2), fbx%hi(2)
               do i = fbx%lo(1), fbx%hi(1)
                  pV(i,j,k,1) = random_uniform(lo=-1.0_WP, hi=1.0_WP)
                  !pV(i,j,k,1) = sin(2.0_WP*Pi*real(j,WP)*solver%amr%dy(lvl))
                  !pV(i,j,k,1) = -cos(2.0_WP*Pi*(real(i,WP)+0.5_WP)*solver%amr%dx(lvl)) * sin(2.0_WP*Pi*real(j,WP)*solver%amr%dy(lvl))
               end do
            end do
         end do

         ! W (z-faces): random values in [-1, 1]
         fbx = mfi%nodaltilebox(3)
         do k = fbx%lo(3), fbx%hi(3)
            do j = fbx%lo(2), fbx%hi(2)
               do i = fbx%lo(1), fbx%hi(1)
                  pW(i,j,k,1) = random_uniform(lo=-1.0_WP, hi=1.0_WP)
                  !pW(i,j,k,1) = sin(2.0_WP*Pi*real(k,WP)*solver%amr%dz(lvl))
                  !pW(i,j,k,1) = 0.0_WP
               end do
            end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine velocity_init

   !> Geometric tagger: refine center box
   subroutine geometric_tagger(solver, lvl, tags_ptr, time)
      use iso_c_binding,    only: c_ptr, c_char
      use amrex_amr_module, only: amrex_tagboxarray
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

      ! Convert c_ptr to tagboxarray
      tags = tags_ptr

      call solver%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3)
            z = solver%amr%zlo + (real(k,WP) + 0.5_WP) * dz
            do j = bx%lo(2), bx%hi(2)
               y = solver%amr%ylo + (real(j,WP) + 0.5_WP) * dy
               do i = bx%lo(1), bx%hi(1)
                  x = solver%amr%xlo + (real(i,WP) + 0.5_WP) * dx
                  ! Tag if inside sphere
                  if (sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2) .lt. radius) then
                     tagarr(i,j,k,1) = SETtag
                  end if
               end do
            end do
         end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine geometric_tagger

   !> Main test routine
   subroutine test_amrincomp()
      use messager, only: log
      use string,   only: itoa, rtoa
      use iso_c_binding, only: c_char
      real(WP) :: divmax_before, divmax_after, dt, time, factor
      integer :: lvl, i, j, k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx, fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW, pP
      real(WP) :: dxi, dyi, dzi

      call log("========================================")
      call log("TEST: amrincomp 3-level projection")
      call log("========================================")

      time = 0.0_WP
      dt = 1.0e-2_WP

      ! Create amrgrid
      allocate(amr)
      amr%name = 'incomp_test'
      amr%nx = 32; amr%ny = 32; amr%nz = 32
      amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
      amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
      amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
      amr%xper = .true.; amr%yper = .true.; amr%zper = .true.
      amr%maxlvl = 2  ! 3 levels: 0, 1, 2
      call amr%initialize()
      call log("Grid initialized: "//trim(itoa(amr%nx))//"x"//trim(itoa(amr%ny))//"x"//trim(itoa(amr%nz)))
      call log("Max levels: "//trim(itoa(amr%maxlvl+1)))

      ! Create flow solver
      allocate(fs)
      call fs%initialize(amr, name='test_fs')
      fs%rho = 1.0_WP   ! Set density directly
      fs%user_init => velocity_init       ! Set velocity init callback
      fs%user_tagging => geometric_tagger ! Set tagging callback
      call log("Flow solver initialized, rho="//trim(rtoa(fs%rho)))

      ! Create workspace for pressure gradients
      allocate(dPdx); call dPdx%initialize(amr, name='dPdx', ncomp=1, ng=0, nodal=[.true., .false., .false.], interp=amrex_interp_none); call dPdx%register()
      allocate(dPdy); call dPdy%initialize(amr, name='dPdy', ncomp=1, ng=0, nodal=[.false., .true., .false.], interp=amrex_interp_none); call dPdy%register()
      allocate(dPdz); call dPdz%initialize(amr, name='dPdz', ncomp=1, ng=0, nodal=[.false., .false., .true.], interp=amrex_interp_none); call dPdz%register()

      ! Initialize visualization
      allocate(viz); call viz%initialize(amr, 'test_incomp')
      call viz%add_scalar(fs%P, 1, 'pressure')
      call viz%add_scalar(fs%div, 1, 'divergence')
      call viz%add_scalar(fs%U, 1, 'U')
      call viz%add_scalar(fs%V, 1, 'V')
      call viz%add_scalar(fs%W, 1, 'W')

      ! Build grid - triggers on_init + user_init callbacks, creates all levels
      call amr%init_from_scratch(time=time)
      call log("Grid built: "//trim(itoa(amr%nlevels))//" levels")

      ! Average down velocity first to make coarse consistent with fine
      call fs%average_down_velocity()

      ! Now fill velocity ghosts (fine level will interpolate from corrected coarse)
      call fs%fill_velocity(time)
      call log("Velocity averaged down and ghosts filled")

      ! Setup pressure solver with verbose output for debugging
      fs%psolver%verbose = 2
      call fs%psolver%setup()
      call log("Pressure solver setup complete")

      ! ===== USER-DRIVEN PROJECTION SEQUENCE =====
      ! 1. Compute divergence (assumes velocity ghosts filled)
      call fs%get_div()
      divmax_before = fs%divmax
      call log("Before projection: divmax = "//trim(rtoa(divmax_before)))
      call log("RHS sum = "//trim(rtoa(fs%div%mf(0)%sum(1))))

      call viz%write(time=0.0_WP)

      ! 2. Scale RHS: rhs = rho/dt * div
      factor = fs%rho / dt
      do lvl = 0, amr%clvl()
         call fs%div%mf(lvl)%mult(factor, 1, 1, 0)
      end do

      ! 3. Solve pressure Poisson equation
      call log("RHS norm0 (after scaling) = "//trim(rtoa(fs%div%mf(0)%norm0(1))))
      call fs%psolver%solve(rhs=fs%div, phi=fs%P)
      call log("Pressure solved, P norm0 = "//trim(rtoa(fs%P%mf(0)%norm0(1))))
      call log("Pressure solver residual = "//trim(rtoa(fs%psolver%res)))

      ! 4. Fill pressure ghosts before gradient computation
      !do lvl = 0, amr%clvl()
      !   call fs%P%fill(lvl, time)
      !end do

      ! 5. Get C/F-consistent gradients and correct velocity
      ! get_fluxes returns -grad(P) for Poisson, so U = U + (dt/rho)*flux
      call fs%psolver%get_fluxes(fs%P, dPdx, dPdy, dPdz)
      call log("dPdx norm0 = "//trim(rtoa(dPdx%mf(0)%norm0(1)))//", dPdz norm0 = "//trim(rtoa(dPdz%mf(0)%norm0(1))))
      factor = dt / fs%rho

      ! Correct velocity: U = U - (dt/rho) * grad(P) = U + factor * dPdx
      do lvl = 0, amr%clvl()
         call amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            pU => fs%U%mf(lvl)%dataptr(mfi)
            pV => fs%V%mf(lvl)%dataptr(mfi)
            pW => fs%W%mf(lvl)%dataptr(mfi)

            ! Get pointers to flux workspace
            block
               real(WP), dimension(:,:,:,:), contiguous, pointer :: pdPdx, pdPdy, pdPdz
               pdPdx => dPdx%mf(lvl)%dataptr(mfi)
               pdPdy => dPdy%mf(lvl)%dataptr(mfi)
               pdPdz => dPdz%mf(lvl)%dataptr(mfi)

               ! U correction (x-faces)
               fbx = mfi%nodaltilebox(1)
               do k = fbx%lo(3), fbx%hi(3)
                  do j = fbx%lo(2), fbx%hi(2)
                     do i = fbx%lo(1), fbx%hi(1)
                        pU(i,j,k,1) = pU(i,j,k,1) + factor * pdPdx(i,j,k,1)
                     end do
                  end do
               end do

               ! V correction (y-faces)
               fbx = mfi%nodaltilebox(2)
               do k = fbx%lo(3), fbx%hi(3)
                  do j = fbx%lo(2), fbx%hi(2)
                     do i = fbx%lo(1), fbx%hi(1)
                        pV(i,j,k,1) = pV(i,j,k,1) + factor * pdPdy(i,j,k,1)
                     end do
                  end do
               end do

               ! W correction (z-faces)
               fbx = mfi%nodaltilebox(3)
               do k = fbx%lo(3), fbx%hi(3)
                  do j = fbx%lo(2), fbx%hi(2)
                     do i = fbx%lo(1), fbx%hi(1)
                        pW(i,j,k,1) = pW(i,j,k,1) + factor * pdPdz(i,j,k,1)
                     end do
                  end do
               end do
            end block
         end do
         call amr%mfiter_destroy(mfi)
      end do
      call log("Velocity corrected")
      call viz%write(time=0.5_WP)
      ! Multi-level sync: use proper MAC face averaging from fine to coarse
      call fs%average_down_velocity()

      ! Sync periodic ghosts AFTER average_down, then check divergence
      call fs%fill_velocity(time)
      call fs%get_div()
      divmax_after = fs%divmax
      call log("After projection: divmax = "//trim(rtoa(divmax_after)))

      ! Report result
      if (divmax_after .lt. 1.0e-9_WP) then
         call log("TEST PASSED: Velocity is divergence-free (divmax < 1e-9)")
      else if (divmax_after .lt. divmax_before * 1.0e-6_WP) then
         call log("TEST PASSED: Divergence reduced by > 6 orders of magnitude")
      else
         call log("TEST FAILED: Divergence not sufficiently reduced")
      end if

      ! Write visualization output
      call viz%write(time=1.0_WP)
      call log("Visualization written to amrviz/test_incomp/")

      ! Cleanup
      call viz%finalize()
      call fs%finalize()
      call amr%finalize()
      deallocate(viz, fs, amr)

      call log("Test complete")

   end subroutine test_amrincomp

end module mod_test_amrincomp
