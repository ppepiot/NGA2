!> Visualization test - create fields, output to ensight, verify in paraview
module mod_test_visualization
   use precision,        only: WP
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrensight_class, only: amrensight
   use messager,         only: log,die
   implicit none
   private
   public :: test_visualization

   ! Module-level amr object for tagger callback access
   type(amrgrid), pointer :: amr_ptr => null()

contains

   !> Simple geometric tagger - tag center of domain
   subroutine box_tagger(lvl, tags_ptr, time)
      use iso_c_binding,    only: c_ptr, c_char
      use amrex_amr_module, only: amrex_tagboxarray, amrex_mfiter, amrex_box
      implicit none
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags_ptr
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      character(kind=c_char), parameter :: SET = char(1)  ! Tag value
      real(WP) :: x, y, z, dx, dy, dz
      integer :: i, j, k

      ! Wrap the C pointer
      tags = tags_ptr

      ! Get mesh spacing at this level
      dx = amr_ptr%geom(lvl)%dx(1)
      dy = amr_ptr%geom(lvl)%dx(2)
      dz = amr_ptr%geom(lvl)%dx(3)

      ! Build mfiter for tagging
      call amr_ptr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)

         ! Tag cells in center of domain
         do k = bx%lo(3), bx%hi(3)
            z = amr_ptr%zlo + (real(k,WP) + 0.5_WP) * dz
            do j = bx%lo(2), bx%hi(2)
               y = amr_ptr%ylo + (real(j,WP) + 0.5_WP) * dy
               do i = bx%lo(1), bx%hi(1)
                  x = amr_ptr%xlo + (real(i,WP) + 0.5_WP) * dx
                  ! Tag center quarter of domain
                  if (x > 0.25_WP .and. x < 0.75_WP .and. &
                     y > 0.25_WP .and. y < 0.75_WP .and. &
                     z > 0.25_WP .and. z < 0.75_WP) then
                     tagarr(i,j,k,1) = SET
                  end if
               end do
            end do
         end do
      end do
      call amr_ptr%mfiter_destroy(mfi)

   end subroutine box_tagger


   subroutine test_visualization()
      use iso_c_binding,    only: c_associated
      use amrex_amr_module, only: amrex_mfiter, amrex_box
      implicit none
      type(amrgrid), target :: amr
      type(amrdata), target :: velocity, pressure
      type(amrensight) :: ens
      type(amrex_mfiter) :: mfi
      type(amrex_box)    :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: vel, p
      integer :: lvl, i, j, k
      real(WP) :: x, y, z, dx, dy, dz

      call log("---------------------------------------------------")
      call log("Running Test: Visualization (Multi-Level)")
      call log("---------------------------------------------------")

      ! Set up grid parameters
      amr%nx   = 32
      amr%ny   = 32
      amr%nz   = 32
      amr%xlo  = 0.0_WP
      amr%xhi  = 1.0_WP
      amr%ylo  = 0.0_WP
      amr%yhi  = 1.0_WP
      amr%zlo  = 0.0_WP
      amr%zhi  = 1.0_WP
      amr%xper = .true.
      amr%yper = .true.
      amr%zper = .true.
      amr%nlvl = 1      ! Allow 1 refinement level
      amr%nmax = 16

      ! Initialize AMR grid
      call amr%initialize("viz_amr")

      ! Set module pointer for tagger callback
      amr_ptr => amr

      ! Set tagging callback
      call amr%set_tagging(box_tagger)

      ! Register data fields
      call amr%register(data=velocity, ncomp=3, ng=1, name="velocity")
      call amr%register(data=pressure, ncomp=1, ng=1, name="pressure")

      ! Build the grid (will call tagger to create refined region)
      call amr%initialize_grid(0.0_WP)
      call amr%get_info()
      call log("Grid built with "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")

      ! Initialize fields with simple test patterns
      do lvl = 0, amr%clvl()
         ! Get mesh spacing at this level
         dx = amr%geom(lvl)%dx(1)
         dy = amr%geom(lvl)%dx(2)
         dz = amr%geom(lvl)%dx(3)

         ! Build iterator at this level
         call amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            vel => velocity%mf(lvl)%dataptr(mfi)
            p   => pressure%mf(lvl)%dataptr(mfi)

            ! Fill with test pattern
            do k = bx%lo(3), bx%hi(3)
               z = amr%zlo + (real(k,WP) + 0.5_WP) * dz
               do j = bx%lo(2), bx%hi(2)
                  y = amr%ylo + (real(j,WP) + 0.5_WP) * dy
                  do i = bx%lo(1), bx%hi(1)
                     x = amr%xlo + (real(i,WP) + 0.5_WP) * dx
                     ! Velocity: vortex-like pattern
                     vel(i,j,k,1) = -sin(2.0_WP * 3.14159_WP * y)
                     vel(i,j,k,2) =  sin(2.0_WP * 3.14159_WP * x)
                     vel(i,j,k,3) = 0.0_WP
                     ! Pressure: sinusoidal
                     p(i,j,k,1) = sin(2.0_WP * 3.14159_WP * x) * &
                        sin(2.0_WP * 3.14159_WP * y) * &
                        sin(2.0_WP * 3.14159_WP * z)
                  end do
               end do
            end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      call log("Fields initialized with test patterns")

      ! Initialize ensight output
      call ens%initialize(amr=amr, name="viz_test")
      call ens%add_scalar(data=pressure, comp=1, name="pressure")
      call ens%add_scalar(data=velocity, comp=1, name="U")
      call ens%add_scalar(data=velocity, comp=2, name="V")
      call ens%add_scalar(data=velocity, comp=3, name="W")
      call ens%add_vector(datax=velocity, compx=1, &
         datay=velocity, compy=2, &
         dataz=velocity, compz=3, name="velocity")
      call log("Ensight output initialized")

      ! Write output
      call ens%write_data(time=0.0_WP)
      call log("Ensight data written to ensight/viz_test/")

      ! Cleanup
      call amr%set_tagging()  ! Clear tagger
      amr_ptr => null()
      call amr%unregister(data=velocity)
      call amr%unregister(data=pressure)
      call amr%finalize()
      call log("PASS: Visualization test complete!")
      call log("Open ensight/viz_test/nga.lev0.case and nga.lev1.case in ParaView")

   contains
      function itoa(i) result(str)
         integer, intent(in) :: i
         character(len=32) :: str
         write(str,'(i0)') i
      end function itoa

   end subroutine test_visualization

end module mod_test_visualization
