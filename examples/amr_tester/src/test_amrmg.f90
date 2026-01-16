!> Test program for unified amrmg solver
!> Tests both constant-coefficient (Poisson) and variable-coefficient modes
!> Uses manufactured solution to verify solver
module mod_test_amrmg
   use precision,         only: WP
   use string,            only: str_long
   use messager,          only: log
   use mathtools,         only: Pi
   use amrgrid_class,     only: amrgrid
   use amrdata_class,     only: amrdata
   use amrmg_class
   use timer_class,       only: timer
   use amrex_amr_module,  only: amrex_mfiter, amrex_box
   use mpi_f08
   implicit none
   private
   public :: test_amrmg

contains

   subroutine test_amrmg()
      call log('=== TEST_AMRMG STARTING ===')
      call log('')

      ! Test constant-coefficient solver (Poisson)
      call test_cstcoef()

      ! Test variable-coefficient solver
      call test_varcoef()

      call log('')
      call log('=== TEST_AMRMG COMPLETE ===')
   end subroutine test_amrmg

   !> Test constant-coefficient (Poisson) solver
   subroutine test_cstcoef()
      type(amrgrid), allocatable, target :: amr
      type(amrdata), allocatable, target :: phi, rhs, exact
      type(amrmg) :: solver
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(timer) :: tmr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPhi, pRhs, pExact
      real(WP) :: x, y, z, dx, phi_exact
      real(WP) :: err_l2, err_linf, local_err
      integer :: lvl, i, j, k, ierr
      character(len=str_long) :: msg

      call log('--- Testing constant-coefficient (Poisson) solver ---')
      call log('Solving: -∇²ϕ = f with ϕ = sin(πx)sin(πy)sin(πz)')

      ! Allocate and configure grid
      allocate(amr)
      amr%nx = 64; amr%ny = 64; amr%nz = 64
      amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
      amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
      amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
      amr%xper = .false.; amr%yper = .false.; amr%zper = .false.
      amr%maxlvl = 0
      amr%nmax = 32
      call amr%initialize(name='poisson_test')
      call amr%init_from_scratch(time=0.0_WP)

      ! Allocate fields
      allocate(phi, rhs, exact)
      call phi%initialize(amr=amr, name='phi', ncomp=1, ng=1)
      call rhs%initialize(amr=amr, name='rhs', ncomp=1, ng=0)
      call exact%initialize(amr=amr, name='exact', ncomp=1, ng=0)
      call phi%on_init(phi, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call rhs%on_init(rhs, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call exact%on_init(exact, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))

      ! Fill fields with manufactured solution
      dx = amr%dx(0)
      call amr%mfiter_build(0, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pPhi => phi%mf(0)%dataptr(mfi)
         pRhs => rhs%mf(0)%dataptr(mfi)
         pExact => exact%mf(0)%dataptr(mfi)

         ! Fill phi ghosts with exact solution (Dirichlet BC)
         do k = bx%lo(3)-1, bx%hi(3)+1
            z = amr%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2)-1, bx%hi(2)+1
               y = amr%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1)-1, bx%hi(1)+1
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  pPhi(i,j,k,1) = sin(pi*x) * sin(pi*y) * sin(pi*z)
               end do
            end do
         end do

         ! Fill valid region
         do k = bx%lo(3), bx%hi(3)
            z = amr%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2), bx%hi(2)
               y = amr%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1), bx%hi(1)
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  phi_exact = sin(pi*x) * sin(pi*y) * sin(pi*z)
                  pExact(i,j,k,1) = phi_exact
                  pRhs(i,j,k,1) = -3.0_WP*pi**2 * phi_exact  ! Poisson: ∇²ϕ = f
               end do
            end do
         end do
      end do
      call amr%mfiter_destroy(mfi)

      ! Initialize and setup solver
      call solver%initialize(amr=amr, type=amrmg_cstcoef)
      solver%lo_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%hi_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%tol_rel = 1.0e-12_WP
      solver%verbose = 0
      solver%bottom_solver = amrmg_bottom_hypre

      tmr = timer(comm=MPI_COMM_WORLD, name='amrmg_cst')
      call tmr%start()
      call solver%setup()
      call tmr%stop()
      write(msg,'(a,es12.5,a)') '  Setup time: ', tmr%time, ' s'; call log(msg)

      ! Solve with zero initial guess
      call phi%mf(0)%setval(0.0_WP)
      call tmr%reset(); call tmr%start()
      call solver%solve(phi=phi, rhs=rhs)
      call tmr%stop()
      write(msg,'(a,es12.5,a,es12.5,a)') '  Solve: residual = ', solver%res, ', time = ', tmr%time, ' s'
      call log(msg)

      ! Compute error
      err_l2 = 0.0_WP; err_linf = 0.0_WP
      call amr%mfiter_build(0, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pPhi => phi%mf(0)%dataptr(mfi)
         pExact => exact%mf(0)%dataptr(mfi)
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
      call amr%mfiter_destroy(mfi)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      err_l2 = sqrt(err_l2)

      write(msg,'(a,es12.5,a,es12.5)') '  Error: L2 = ', err_l2, ', Linf = ', err_linf; call log(msg)

      if (err_linf .lt. 0.1_WP .and. solver%res .lt. 1.0e-8_WP) then
         call log('  PASS: constant-coefficient (Poisson) test')
      else
         call log('  FAIL: constant-coefficient (Poisson) test')
      end if

      ! Cleanup
      call solver%finalize()
      call phi%finalize(); call rhs%finalize(); call exact%finalize()
      deallocate(phi, rhs, exact)
      call amr%finalize(); deallocate(amr)
   end subroutine test_cstcoef

   !> Test variable-coefficient solver
   subroutine test_varcoef()
      type(amrgrid), allocatable, target :: amr
      type(amrdata), allocatable, target :: phi, rhs, exact
      type(amrdata), allocatable, target :: bcoef_x, bcoef_y, bcoef_z
      type(amrmg) :: solver
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(timer) :: tmr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPhi, pRhs, pExact, pBx, pBy, pBz
      real(WP) :: x, y, z, dx, phi_exact
      real(WP) :: err_l2, err_linf, local_err
      integer :: lvl, i, j, k, ierr
      character(len=str_long) :: msg

      call log('')
      call log('--- Testing variable-coefficient solver ---')
      call log('Solving: -∇·(B∇ϕ) = f with B = 1 (uniform)')

      ! Allocate and configure grid
      allocate(amr)
      amr%nx = 64; amr%ny = 64; amr%nz = 64
      amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
      amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
      amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
      amr%xper = .false.; amr%yper = .false.; amr%zper = .false.
      amr%maxlvl = 0
      amr%nmax = 32
      call amr%initialize(name='varcoef_test')
      call amr%init_from_scratch(time=0.0_WP)

      ! Allocate fields
      allocate(phi, rhs, exact, bcoef_x, bcoef_y, bcoef_z)
      call phi%initialize(amr=amr, name='phi', ncomp=1, ng=1)
      call rhs%initialize(amr=amr, name='rhs', ncomp=1, ng=0)
      call exact%initialize(amr=amr, name='exact', ncomp=1, ng=0)
      call bcoef_x%initialize(amr=amr, name='bx', ncomp=1, ng=0, nodal=[.true.,.false.,.false.])
      call bcoef_y%initialize(amr=amr, name='by', ncomp=1, ng=0, nodal=[.false.,.true.,.false.])
      call bcoef_z%initialize(amr=amr, name='bz', ncomp=1, ng=0, nodal=[.false.,.false.,.true.])

      call phi%on_init(phi, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call rhs%on_init(rhs, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call exact%on_init(exact, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call bcoef_x%on_init(bcoef_x, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call bcoef_y%on_init(bcoef_y, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call bcoef_z%on_init(bcoef_z, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))

      ! Fill fields with manufactured solution
      dx = amr%dx(0)
      call amr%mfiter_build(0, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pPhi => phi%mf(0)%dataptr(mfi)
         pRhs => rhs%mf(0)%dataptr(mfi)
         pExact => exact%mf(0)%dataptr(mfi)

         ! Fill phi ghosts with exact solution (Dirichlet BC)
         do k = bx%lo(3)-1, bx%hi(3)+1
            z = amr%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2)-1, bx%hi(2)+1
               y = amr%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1)-1, bx%hi(1)+1
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  pPhi(i,j,k,1) = sin(pi*x) * sin(pi*y) * sin(pi*z)
               end do
            end do
         end do

         ! Fill valid region
         do k = bx%lo(3), bx%hi(3)
            z = amr%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2), bx%hi(2)
               y = amr%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1), bx%hi(1)
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  phi_exact = sin(pi*x) * sin(pi*y) * sin(pi*z)
                  pExact(i,j,k,1) = phi_exact
                  pRhs(i,j,k,1) = 3.0_WP*pi**2 * phi_exact
               end do
            end do
         end do
      end do
      call amr%mfiter_destroy(mfi)

      ! Fill B coefficients (uniform = 1)
      call bcoef_x%mfiter_build(0, mfi)
      do while (mfi%next())
         pBx => bcoef_x%mf(0)%dataptr(mfi)
         bx = mfi%tilebox()
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  pBx(i,j,k,1) = 1.0_WP
               end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      call bcoef_y%mfiter_build(0, mfi)
      do while (mfi%next())
         pBy => bcoef_y%mf(0)%dataptr(mfi)
         bx = mfi%tilebox()
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  pBy(i,j,k,1) = 1.0_WP
               end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      call bcoef_z%mfiter_build(0, mfi)
      do while (mfi%next())
         pBz => bcoef_z%mf(0)%dataptr(mfi)
         bx = mfi%tilebox()
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  pBz(i,j,k,1) = 1.0_WP
               end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)

      ! Initialize and setup solver
      call solver%initialize(amr=amr, type=amrmg_varcoef)
      solver%lo_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%hi_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%alpha = 0.0_WP
      solver%beta  = 1.0_WP
      solver%tol_rel = 1.0e-12_WP
      solver%verbose = 0
      solver%bottom_solver = amrmg_bottom_hypre

      tmr = timer(comm=MPI_COMM_WORLD, name='amrmg_var')
      call tmr%start()
      call solver%setup(bcoef_x=bcoef_x, bcoef_y=bcoef_y, bcoef_z=bcoef_z)
      call tmr%stop()
      write(msg,'(a,es12.5,a)') '  Setup time: ', tmr%time, ' s'; call log(msg)

      ! First solve
      call phi%mf(0)%setval(0.0_WP)
      call tmr%reset(); call tmr%start()
      call solver%solve(phi=phi, rhs=rhs)
      call tmr%stop()
      write(msg,'(a,es12.5,a,es12.5,a)') '  First solve: residual = ', solver%res, ', time = ', tmr%time, ' s'
      call log(msg)

      ! Second solve (reuse test)
      call phi%mf(0)%setval(0.0_WP)
      call tmr%reset(); call tmr%start()
      call solver%solve(phi=phi, rhs=rhs)
      call tmr%stop()
      write(msg,'(a,es12.5,a,es12.5,a)') '  Second solve (reuse): residual = ', solver%res, ', time = ', tmr%time, ' s'
      call log(msg)

      ! Compute error
      err_l2 = 0.0_WP; err_linf = 0.0_WP
      call amr%mfiter_build(0, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pPhi => phi%mf(0)%dataptr(mfi)
         pExact => exact%mf(0)%dataptr(mfi)
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
      call amr%mfiter_destroy(mfi)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      err_l2 = sqrt(err_l2)

      write(msg,'(a,es12.5,a,es12.5)') '  Error: L2 = ', err_l2, ', Linf = ', err_linf; call log(msg)

      if (err_linf .lt. 0.1_WP .and. solver%res .lt. 1.0e-8_WP) then
         call log('  PASS: variable-coefficient test')
      else
         call log('  FAIL: variable-coefficient test')
      end if

      ! Cleanup
      call solver%finalize()
      call phi%finalize(); call rhs%finalize(); call exact%finalize()
      call bcoef_x%finalize(); call bcoef_y%finalize(); call bcoef_z%finalize()
      deallocate(phi, rhs, exact, bcoef_x, bcoef_y, bcoef_z)
      call amr%finalize(); deallocate(amr)
   end subroutine test_varcoef

end module mod_test_amrmg
