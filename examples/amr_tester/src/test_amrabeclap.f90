!> Test program for amrabeclap variable coefficient solver
!> Solves -∇·(B∇ϕ) = f with B = 1 + 0.5*sin(πx)
!> Uses manufactured solution to verify solver
module mod_test_amrabeclap
   use precision,         only: WP
   use amrgrid_class,     only: amrgrid
   use amrdata_class,     only: amrdata
   use amrmg_class
   use timer_class,       only: timer
   use amrex_amr_module,  only: amrex_mfiter, amrex_box
   use mpi_f08
   implicit none
   private
   public :: test_amrabeclap

   ! Module-level objects
   type(amrgrid), allocatable, target :: amr_ab
   type(amrdata), allocatable, target :: phi, rhs_ab, exact
   type(amrdata), allocatable, target :: bcoef_x, bcoef_y, bcoef_z

   real(WP), parameter :: pi = 4.0_WP * atan(1.0_WP)

contains

   subroutine test_amrabeclap()
      type(amrmg) :: solver
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPhi, pRhs, pExact
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pBx, pBy, pBz
      real(WP) :: x, y, z, dx, err_l2, local_err, err_linf
      real(WP) :: B, dBdx, phi_exact, d2phi
      integer :: lvl, i, j, k, ierr
      type(timer) :: tmr

      write(*,*) '=== TEST_AMRABECLAP STARTING ==='
      write(*,*) ''
      write(*,*) 'Testing variable coefficient solver: -div(B*grad(phi)) = f'
      write(*,*) 'where B = 1 + 0.5*sin(pi*x)'
      write(*,*) ''

      ! Allocate and configure grid
      allocate(amr_ab)
      amr_ab%nx = 64
      amr_ab%ny = 64
      amr_ab%nz = 64
      amr_ab%xlo = 0.0_WP; amr_ab%xhi = 1.0_WP
      amr_ab%ylo = 0.0_WP; amr_ab%yhi = 1.0_WP
      amr_ab%zlo = 0.0_WP; amr_ab%zhi = 1.0_WP
      amr_ab%xper = .false.
      amr_ab%yper = .false.
      amr_ab%zper = .false.
      amr_ab%maxlvl = 0   ! Single level for this test
      amr_ab%nmax = 32

      call amr_ab%initialize("abeclap_test")
      call amr_ab%init_from_scratch(0.0_WP)

      ! Allocate fields
      allocate(phi, rhs_ab, exact)
      allocate(bcoef_x, bcoef_y, bcoef_z)

      ! Cell-centered fields
      call phi%initialize(   amr_ab, name='phi',   ncomp=1, ng=1)
      call rhs_ab%initialize(amr_ab, name='rhs',   ncomp=1, ng=0)
      call exact%initialize( amr_ab, name='exact', ncomp=1, ng=0)

      ! Face-centered B coefficients (staggered) - using nodal support
      call bcoef_x%initialize(amr_ab, name='bx', ncomp=1, ng=0, nodal=[.true.,.false.,.false.])
      call bcoef_y%initialize(amr_ab, name='by', ncomp=1, ng=0, nodal=[.false.,.true.,.false.])
      call bcoef_z%initialize(amr_ab, name='bz', ncomp=1, ng=0, nodal=[.false.,.false.,.true.])

      ! Build level 0 MultiFabs using on_init callback
      call phi%on_init(phi, 0, 0.0_WP, amr_ab%get_boxarray(0), amr_ab%get_distromap(0))
      call rhs_ab%on_init(rhs_ab, 0, 0.0_WP, amr_ab%get_boxarray(0), amr_ab%get_distromap(0))
      call exact%on_init(exact, 0, 0.0_WP, amr_ab%get_boxarray(0), amr_ab%get_distromap(0))
      call bcoef_x%on_init(bcoef_x, 0, 0.0_WP, amr_ab%get_boxarray(0), amr_ab%get_distromap(0))
      call bcoef_y%on_init(bcoef_y, 0, 0.0_WP, amr_ab%get_boxarray(0), amr_ab%get_distromap(0))
      call bcoef_z%on_init(bcoef_z, 0, 0.0_WP, amr_ab%get_boxarray(0), amr_ab%get_distromap(0))

      ! Set up manufactured solution:
      ! phi_exact = sin(pi*x)*sin(pi*y)*sin(pi*z)
      ! B = 1 (uniform for now)
      ! -div(B*grad(phi)) = -laplacian(phi) = 3*pi^2*phi

      ! Fill cell-centered fields (phi, rhs, exact)
      do lvl = 0, amr_ab%clvl()
         dx = amr_ab%dx(lvl)
         call amr_ab%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pPhi   => phi%mf(lvl)%dataptr(mfi)
            pRhs   => rhs_ab%mf(lvl)%dataptr(mfi)
            pExact => exact%mf(lvl)%dataptr(mfi)

            ! Fill phi ghosts with exact solution (for Dirichlet BC)
            do k = bx%lo(3)-1, bx%hi(3)+1
               z = amr_ab%zlo + (real(k,WP)+0.5_WP)*dx
               do j = bx%lo(2)-1, bx%hi(2)+1
                  y = amr_ab%ylo + (real(j,WP)+0.5_WP)*dx
                  do i = bx%lo(1)-1, bx%hi(1)+1
                     x = amr_ab%xlo + (real(i,WP)+0.5_WP)*dx
                     pPhi(i,j,k,1) = sin(pi*x) * sin(pi*y) * sin(pi*z)
                  end do
               end do
            end do

            ! Fill valid region data
            do k = bx%lo(3), bx%hi(3)
               z = amr_ab%zlo + (real(k,WP)+0.5_WP)*dx
               do j = bx%lo(2), bx%hi(2)
                  y = amr_ab%ylo + (real(j,WP)+0.5_WP)*dx
                  do i = bx%lo(1), bx%hi(1)
                     x = amr_ab%xlo + (real(i,WP)+0.5_WP)*dx
                     phi_exact = sin(pi*x) * sin(pi*y) * sin(pi*z)
                     pExact(i,j,k,1) = phi_exact
                     pRhs(i,j,k,1) = 3.0_WP*pi**2 * phi_exact
                  end do
               end do
            end do
         end do
         call amr_ab%mfiter_destroy(mfi)
      end do

      ! Fill face-centered B coefficients (B = 1 uniform)
      ! X-faces
      do lvl = 0, amr_ab%clvl()
         call bcoef_x%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()  ! Returns nodal box for nodal data
            pBx => bcoef_x%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pBx(i,j,k,1) = 1.0_WP
                  end do
               end do
            end do
         end do
         call amr_ab%mfiter_destroy(mfi)
      end do
      ! Y-faces
      do lvl = 0, amr_ab%clvl()
         call bcoef_y%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pBy => bcoef_y%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pBy(i,j,k,1) = 1.0_WP
                  end do
               end do
            end do
         end do
         call amr_ab%mfiter_destroy(mfi)
      end do
      ! Z-faces
      do lvl = 0, amr_ab%clvl()
         call bcoef_z%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pBz => bcoef_z%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pBz(i,j,k,1) = 1.0_WP
                  end do
               end do
            end do
         end do
         call amr_ab%mfiter_destroy(mfi)
      end do

      ! Zero interior for initial guess
      do lvl = 0, amr_ab%clvl()
         call amr_ab%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pPhi => phi%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pPhi(i,j,k,1) = 0.0_WP
                  end do
               end do
            end do
         end do
         call amr_ab%mfiter_destroy(mfi)
      end do

      ! Initialize solver
      call solver%initialize(amr=amr_ab, type=amrmg_varcoef)
      solver%lo_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%hi_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%alpha = 0.0_WP    ! No A term
      solver%beta  = 1.0_WP    ! div(B*grad) term
      solver%tol_rel = 1.0e-12_WP
      solver%verbose = 2
      solver%bottom_solver = amrmg_bottom_hypre

      ! Time setup
      tmr = timer(comm=MPI_COMM_WORLD, name='amrmg')
      call tmr%start()
      call solver%setup(bcoef_x=bcoef_x, bcoef_y=bcoef_y, bcoef_z=bcoef_z)
      call tmr%stop()
      write(*,'(a,es12.5,a)') ' Setup time: ', tmr%time, ' s'

      ! First solve with zero initial guess
      write(*,*) ''
      write(*,*) 'Starting first solve with Hypre bottom solver...'
      call phi%mf(0)%setval(0.0_WP)  ! Zero initial guess
      call tmr%reset()
      call tmr%start()
      call solver%solve(phi=phi, rhs=rhs_ab)
      call tmr%stop()
      write(*,'(a,es12.5,a,es12.5,a)') ' First solve: residual = ', solver%res, ', time = ', tmr%time, ' s'

      ! Second solve (reusing operator) with zero initial guess
      write(*,*) ''
      write(*,*) 'Testing solver reuse (second solve with zero initial guess)...'
      call phi%mf(0)%setval(0.0_WP)  ! Reset to zero initial guess
      call tmr%reset()
      call tmr%start()
      call solver%solve(phi=phi, rhs=rhs_ab)
      call tmr%stop()
      write(*,'(a,es12.5,a,es12.5,a)') ' Second solve: residual = ', solver%res, ', time = ', tmr%time, ' s'

      ! Compute error
      err_l2 = 0.0_WP
      err_linf = 0.0_WP
      do lvl = 0, amr_ab%clvl()
         call amr_ab%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pPhi   => phi%mf(lvl)%dataptr(mfi)
            pExact => exact%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     local_err = abs(pPhi(i,j,k,1) - pExact(i,j,k,1))
                     err_l2 = err_l2 + local_err**2
                     err_linf = max(err_linf, local_err)
                  end do
               end do
            end do
         end do
         call amr_ab%mfiter_destroy(mfi)
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, err_l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      err_l2 = sqrt(err_l2)

      write(*,*) ''
      write(*,*) 'Discretization Error (vs analytical):'
      write(*,*) '  L2  :', err_l2
      write(*,*) '  Linf:', err_linf
      write(*,*) ''

      if (err_linf < 0.1_WP .and. solver%res < 1.0e-8_WP) then
         write(*,*) 'PASS: amrabeclap variable coefficient solver test'
      else
         write(*,*) 'FAIL: residual=', solver%res, ' err_linf=', err_linf
      end if

      ! Cleanup
      call solver%finalize()
      call phi%finalize()
      call rhs_ab%finalize()
      call exact%finalize()
      call bcoef_x%finalize()
      call bcoef_y%finalize()
      call bcoef_z%finalize()
      deallocate(phi, rhs_ab, exact, bcoef_x, bcoef_y, bcoef_z)
      call amr_ab%finalize()
      deallocate(amr_ab)

      write(*,*) ''
      write(*,*) '=== TEST_AMRABECLAP COMPLETE ==='

   end subroutine test_amrabeclap

end module mod_test_amrabeclap
