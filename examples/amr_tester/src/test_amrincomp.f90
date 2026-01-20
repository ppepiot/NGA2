!> Test amrincomp solver - 3-level projection with random velocity and geometric tagger
module mod_test_amrincomp
   use precision,         only: WP
   use mathtools,         only: Pi
   use amrgrid_class,     only: amrgrid
   use amrincomp_class,   only: amrincomp
   use amrex_amr_module,  only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_distromap, &
   &                            amrex_mfiter_build, amrex_mfiter_destroy
   implicit none
   private
   public :: test_amrincomp

   ! Module-level objects
   type(amrgrid), allocatable, target :: amr
   type(amrincomp), allocatable, target :: fs

contains

   !> User-provided initialization callback for velocity - sinusoidal values
   subroutine velocity_init(solver, lvl, time, ba, dm)
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx, fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      integer :: i, j, k

      ! Build mfiter from cell-centered ba/dm (passed during regrid)
      call amrex_mfiter_build(mfi, ba, dm, tiling=solver%amr%default_tiling)
      do while (mfi%next())
         bx = mfi%tilebox()
         pU => solver%U%mf(lvl)%dataptr(mfi)
         pV => solver%V%mf(lvl)%dataptr(mfi)
         pW => solver%W%mf(lvl)%dataptr(mfi)

         ! U (x-faces): nodaltilebox in x
         fbx = mfi%nodaltilebox(1)
         do k = fbx%lo(3), fbx%hi(3)
            do j = fbx%lo(2), fbx%hi(2)
               do i = fbx%lo(1), fbx%hi(1)
                  pU(i,j,k,1) = sin(2.0_WP*Pi*real(i,WP)/32.0_WP)
               end do
            end do
         end do

         ! V (y-faces): nodaltilebox in y
         fbx = mfi%nodaltilebox(2)
         do k = fbx%lo(3), fbx%hi(3)
            do j = fbx%lo(2), fbx%hi(2)
               do i = fbx%lo(1), fbx%hi(1)
                  pV(i,j,k,1) = sin(2.0_WP*Pi*real(j,WP)/32.0_WP)
               end do
            end do
         end do

         ! W (z-faces): nodaltilebox in z
         fbx = mfi%nodaltilebox(3)
         do k = fbx%lo(3), fbx%hi(3)
            do j = fbx%lo(2), fbx%hi(2)
               do i = fbx%lo(1), fbx%hi(1)
                  pW(i,j,k,1) = sin(2.0_WP*Pi*real(k,WP)/32.0_WP)
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
                  if (sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2) < radius) then
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
      dt = 0.01_WP

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

      ! Build grid - triggers on_init + user_init callbacks, creates all levels
      call amr%init_from_scratch(time=time)
      call log("Grid built: "//trim(itoa(amr%nlevels))//" levels")

      ! Fill velocity ghosts after initialization (required before get_div)
      do lvl = 0, amr%clvl()
         call fs%U%fill(lvl, time)
         call fs%V%fill(lvl, time)
         call fs%W%fill(lvl, time)
      end do

      ! Average down velocity for C/F consistency before computing divergence
      call fs%average_down_velocity()
      call log("Velocity ghosts filled and averaged down")

      ! Setup pressure solver with verbose output for debugging
      fs%psolver%verbose = 2
      call fs%psolver%setup()
      call log("Pressure solver setup complete")

      ! ===== USER-DRIVEN PROJECTION SEQUENCE =====

      ! 1. Compute divergence (assumes velocity ghosts filled)
      call fs%get_div()
      divmax_before = fs%divmax
      call log("Before projection: divmax = "//trim(rtoa(divmax_before)))

      ! 2. Scale RHS: rhs = rho/dt * div
      factor = fs%rho / dt
      do lvl = 0, amr%clvl()
         call fs%div%mf(lvl)%mult(factor, 1, 1, 0)
      end do

      ! 3. Solve pressure Poisson equation
      call fs%psolver%solve(rhs=fs%div, phi=fs%P)
      call log("Pressure solved, P norm0 = "//trim(rtoa(fs%P%mf(0)%norm0(1))))
      call log("Pressure solver residual = "//trim(rtoa(fs%psolver%res)))

      ! 4. Fill pressure ghosts before gradient computation
      do lvl = 0, amr%clvl()
         call fs%P%fill(lvl, time)
      end do



      ! 5. Correct velocity: U = U - dt/rho * dP/dx (manual gradient)
      factor = dt / fs%rho
      do lvl = 0, amr%clvl()
         dxi = 1.0_WP / amr%dx(lvl)
         dyi = 1.0_WP / amr%dy(lvl)
         dzi = 1.0_WP / amr%dz(lvl)

         ! Correct velocity using cell-centered mfiter with nodalize
         call amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pU => fs%U%mf(lvl)%dataptr(mfi)
            pV => fs%V%mf(lvl)%dataptr(mfi)
            pW => fs%W%mf(lvl)%dataptr(mfi)
            pP => fs%P%mf(lvl)%dataptr(mfi)

            ! U correction (x-faces): nodaltilebox in x
            fbx = mfi%nodaltilebox(1)
            do k = fbx%lo(3), fbx%hi(3)
               do j = fbx%lo(2), fbx%hi(2)
                  do i = fbx%lo(1), fbx%hi(1)
                     pU(i,j,k,1) = pU(i,j,k,1) - factor * (pP(i,j,k,1) - pP(i-1,j,k,1)) * dxi
                  end do
               end do
            end do

            ! V correction (y-faces): nodaltilebox in y
            fbx = mfi%nodaltilebox(2)
            do k = fbx%lo(3), fbx%hi(3)
               do j = fbx%lo(2), fbx%hi(2)
                  do i = fbx%lo(1), fbx%hi(1)
                     pV(i,j,k,1) = pV(i,j,k,1) - factor * (pP(i,j,k,1) - pP(i,j-1,k,1)) * dyi
                  end do
               end do
            end do

            ! W correction (z-faces): nodaltilebox in z
            fbx = mfi%nodaltilebox(3)
            do k = fbx%lo(3), fbx%hi(3)
               do j = fbx%lo(2), fbx%hi(2)
                  do i = fbx%lo(1), fbx%hi(1)
                     pW(i,j,k,1) = pW(i,j,k,1) - factor * (pP(i,j,k,1) - pP(i,j,k-1,1)) * dzi
                  end do
               end do
            end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      call log("Velocity corrected")

      ! Multi-level sync: use proper MAC face averaging from fine to coarse
      call fs%average_down_velocity()

      ! Sync periodic ghosts AFTER average_down, then check divergence
      do lvl = 0, amr%clvl()
         call fs%U%fill(lvl, time)
         call fs%V%fill(lvl, time)
         call fs%W%fill(lvl, time)
      end do
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

      ! Cleanup
      call fs%finalize()
      call amr%finalize()
      deallocate(fs, amr)

      call log("Test complete")

   end subroutine test_amrincomp

end module mod_test_amrincomp
