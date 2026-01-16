!> Unified AMR Multigrid solver for elliptic problems
!> Supports both constant-coefficient (Poisson) and variable-coefficient problems.
!>   - cstcoef type: Uses amrex_poisson for ∇²ϕ = f
!>   - varcoef type: Uses amrex_abeclaplacian for (α·A - β·∇·(B∇))ϕ = f
module amrmg_class
   use precision,     only: WP
   use amrgrid_class, only: amrgrid
   use amrdata_class, only: amrdata
   use amrex_amr_module, only: amrex_multifab, amrex_geometry, amrex_boxarray, amrex_distromap
   use amrex_lo_bctypes_module
   use amrex_linear_solver_module
   implicit none
   private

   ! Solver type constants (required at init)
   integer, parameter, public :: amrmg_cstcoef = 1  !< Constant coefficient (amrex_poisson)
   integer, parameter, public :: amrmg_varcoef = 2  !< Variable coefficient (amrex_abeclaplacian)

   ! Bottom solver types
   integer, parameter, public :: amrmg_bottom_smoother = 0
   integer, parameter, public :: amrmg_bottom_bicgstab = 1
   integer, parameter, public :: amrmg_bottom_cg       = 2
   integer, parameter, public :: amrmg_bottom_hypre    = 3
   integer, parameter, public :: amrmg_bottom_petsc    = 4

   ! Expose AMReX boundary condition types for user convenience
   integer, parameter, public :: amrmg_bc_interior    = amrex_lo_bogus      !< Interior (not at domain boundary)
   integer, parameter, public :: amrmg_bc_dirichlet   = amrex_lo_dirichlet  !< Dirichlet BC
   integer, parameter, public :: amrmg_bc_neumann     = amrex_lo_neumann    !< Neumann BC
   integer, parameter, public :: amrmg_bc_reflect_odd = amrex_lo_reflect_odd
   integer, parameter, public :: amrmg_bc_periodic    = amrex_lo_periodic   !< Periodic BC

   !> Unified AMR Multigrid solver
   type, public :: amrmg
      ! Solver type (required at initialize)
      integer :: type = -1  !< amrmg_cstcoef or amrmg_varcoef (no default)

      ! Equation coefficients (varcoef type only)
      real(WP) :: alpha = 0.0_WP  !< Scalar for A term
      real(WP) :: beta  = 1.0_WP  !< Scalar for ∇·(B∇) term

      ! Solver parameters
      real(WP) :: tol_rel = 1.0e-10_WP
      real(WP) :: tol_abs = 0.0_WP
      integer  :: max_iter = 200
      integer  :: verbose = 0
      integer  :: bottom_solver = amrmg_bottom_bicgstab
      integer  :: maxorder = 2  !< BC/interpolation order

      ! Boundary conditions (set from grid periodicity at init)
      integer, dimension(3) :: lo_bc = amrmg_bc_interior
      integer, dimension(3) :: hi_bc = amrmg_bc_interior

      ! Output from solve
      real(WP) :: res   = 0.0_WP  !< Final residual
      integer  :: niter = 0       !< Number of iterations

      ! Private internals
      class(amrgrid), pointer, private :: amr => null()
      logical, private :: setup_done = .false.
      ! Stored AMReX objects for reuse
      type(amrex_poisson), private :: poisson
      type(amrex_abeclaplacian), private :: abeclap
      type(amrex_multigrid), private :: multigrid
   contains
      procedure :: initialize
      procedure :: setup
      procedure :: solve
      procedure :: solve_level
      procedure :: destroy
      procedure :: print_short
      procedure :: print => print_long
      procedure :: log => log_solver
      procedure :: finalize
   end type amrmg

contains

   !> Initialize solver with grid and type
   subroutine initialize(this, amr, type)
      use messager, only: die, log
      class(amrmg), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      integer, intent(in) :: type  !< amrmg_cstcoef or amrmg_varcoef (required)

      ! Validate type
      if (type .ne. amrmg_cstcoef .and. type .ne. amrmg_varcoef) then
         call die('[amrmg initialize] Invalid type: must be amrmg_cstcoef or amrmg_varcoef')
      end if

      this%amr => amr
      this%type = type
      this%setup_done = .false.

      ! Set smart BC defaults based on grid periodicity
      ! Periodic directions: amrmg_bc_periodic
      ! Non-periodic directions: amrmg_bc_neumann
      if (amr%xper) then
         this%lo_bc(1) = amrmg_bc_periodic
         this%hi_bc(1) = amrmg_bc_periodic
      else
         this%lo_bc(1) = amrmg_bc_neumann
         this%hi_bc(1) = amrmg_bc_neumann
      end if
      if (amr%yper) then
         this%lo_bc(2) = amrmg_bc_periodic
         this%hi_bc(2) = amrmg_bc_periodic
      else
         this%lo_bc(2) = amrmg_bc_neumann
         this%hi_bc(2) = amrmg_bc_neumann
      end if
      if (amr%zper) then
         this%lo_bc(3) = amrmg_bc_periodic
         this%hi_bc(3) = amrmg_bc_periodic
      else
         this%lo_bc(3) = amrmg_bc_neumann
         this%hi_bc(3) = amrmg_bc_neumann
      end if

      if (type .eq. amrmg_cstcoef) then
         call log('[amrmg] Initialized constant-coefficient solver')
      else
         call log('[amrmg] Initialized variable-coefficient solver')
      end if
   end subroutine initialize

   !> Setup solver - call after coefficients change or after regrid
   !> For varcoef type, bcoef fields must be provided
   subroutine setup(this, acoef, bcoef_x, bcoef_y, bcoef_z)
      use messager, only: die
      class(amrmg), intent(inout) :: this
      type(amrdata), intent(in), optional :: acoef
      type(amrdata), intent(in), optional :: bcoef_x, bcoef_y, bcoef_z

      type(amrex_geometry), dimension(:), allocatable :: geom
      type(amrex_boxarray), dimension(:), allocatable :: ba
      type(amrex_distromap), dimension(:), allocatable :: dm
      type(amrex_multifab) :: bcoef(3)
      integer :: lev

      if (this%type .eq. -1) call die('[amrmg setup] Solver not initialized')

      ! If already setup, destroy first
      if (this%setup_done) call this%destroy()

      ! Build AMReX wrapper arrays from amrgrid (use clvl() for current finest level)
      allocate(geom(0:this%amr%clvl()))
      allocate(ba(0:this%amr%clvl()))
      allocate(dm(0:this%amr%clvl()))

      do lev = 0, this%amr%clvl()
         geom(lev) = this%amr%geom(lev)
         ba(lev) = this%amr%get_boxarray(lev)
         dm(lev) = this%amr%get_distromap(lev)
      end do

      ! Build operator and multigrid based on type
      select case (this%type)

       case (amrmg_cstcoef)
         ! Build Poisson operator
         call amrex_poisson_build(this%poisson, geom, ba, dm, &
            metric_term=.false., agglomeration=.true., consolidation=.true., &
            max_coarsening_level=30)
         call this%poisson%set_domain_bc(this%lo_bc, this%hi_bc)

         ! Build multigrid
         call amrex_multigrid_build(this%multigrid, this%poisson)
         call this%multigrid%set_verbose(this%verbose)
         call this%multigrid%set_max_iter(this%max_iter)
         call this%multigrid%set_bottom_solver(this%bottom_solver)

       case (amrmg_varcoef)
         ! Build ABecLaplacian operator
         call amrex_abeclaplacian_build(this%abeclap, geom, ba, dm, &
            metric_term=.false., agglomeration=.true., consolidation=.true., &
            max_coarsening_level=30)
         call this%abeclap%set_domain_bc(this%lo_bc, this%hi_bc)
         call this%abeclap%set_maxorder(this%maxorder)
         call this%abeclap%set_scalars(this%alpha, this%beta)

         ! Set spatially-varying A coefficient if provided
         if (present(acoef)) then
            do lev = 0, this%amr%clvl()
               call this%abeclap%set_acoeffs(lev, acoef%mf(lev))
            end do
         end if

         ! Set spatially-varying B coefficients if provided
         if (present(bcoef_x) .and. present(bcoef_y) .and. present(bcoef_z)) then
            do lev = 0, this%amr%clvl()
               bcoef(1) = bcoef_x%mf(lev)
               bcoef(2) = bcoef_y%mf(lev)
               bcoef(3) = bcoef_z%mf(lev)
               call this%abeclap%set_bcoeffs(lev, bcoef)
            end do
         end if

         ! Build multigrid
         call amrex_multigrid_build(this%multigrid, this%abeclap)
         call this%multigrid%set_verbose(this%verbose)
         call this%multigrid%set_max_iter(this%max_iter)
         call this%multigrid%set_bottom_solver(this%bottom_solver)

      end select

      this%setup_done = .true.
   end subroutine setup

   !> Solve - uses stored operator and multigrid
   !> @param phi Solution field (in: initial guess with BC in ghosts, out: solution)
   !> @param rhs Right-hand side field
   subroutine solve(this, phi, rhs)
      use messager, only: die
      class(amrmg), intent(inout) :: this
      type(amrdata), intent(inout) :: phi
      type(amrdata), intent(in) :: rhs

      type(amrex_multifab), dimension(:), allocatable :: sol, rhsmf
      integer :: lev

      if (this%type .eq. -1) call die('[amrmg solve] Solver not initialized')
      if (.not. this%setup_done) call die('[amrmg solve] Solver not setup')

      ! Build solution/rhs arrays
      allocate(sol(0:this%amr%clvl()))
      allocate(rhsmf(0:this%amr%clvl()))
      do lev = 0, this%amr%clvl()
         sol(lev) = phi%mf(lev)
         rhsmf(lev) = rhs%mf(lev)
      end do

      ! Set level BCs from solution ghost cells
      select case (this%type)
       case (amrmg_cstcoef)
         do lev = 0, this%amr%clvl()
            call this%poisson%set_level_bc(lev, sol(lev))
         end do
       case (amrmg_varcoef)
         do lev = 0, this%amr%clvl()
            call this%abeclap%set_level_bc(lev, sol(lev))
         end do
      end select

      ! Solve
      this%res = this%multigrid%solve(sol, rhsmf, this%tol_rel, this%tol_abs)
      ! Note: AMReX Fortran interface doesn't expose iteration count yet
      this%niter = 0
   end subroutine solve

   !> Level-by-level solve (for subcycling)
   subroutine solve_level(this, lev, phi, rhs, phi_crse)
      use messager, only: die
      class(amrmg), intent(inout) :: this
      integer, intent(in) :: lev
      type(amrdata), intent(inout) :: phi
      type(amrdata), intent(in) :: rhs
      type(amrdata), intent(in), optional :: phi_crse

      call die('[amrmg solve_level] Not yet implemented - use solve()')
   end subroutine solve_level

   !> Destroy solver internals - call before setup when operator changes
   subroutine destroy(this)
      class(amrmg), intent(inout) :: this
      if (this%setup_done) then
         call amrex_multigrid_destroy(this%multigrid)
         select case (this%type)
          case (amrmg_cstcoef)
            call amrex_poisson_destroy(this%poisson)
          case (amrmg_varcoef)
            call amrex_abeclaplacian_destroy(this%abeclap)
         end select
      end if
      this%setup_done = .false.
   end subroutine destroy

   !> Log solver info
   subroutine log_solver(this)
      use string,   only: str_long
      use messager, only: log
      class(amrmg), intent(in) :: this
      character(len=str_long) :: message
      if (this%amr%amRoot) then
         write(message,'("AMR Multigrid solver for grid [",a,"]")') trim(this%amr%name); call log(message)
         if (this%type .eq. amrmg_cstcoef) then
            write(message,'(" >     type = constant-coefficient (Poisson)")'); call log(message)
         else if (this%type .eq. amrmg_varcoef) then
            write(message,'(" >     type = variable-coefficient (ABecLaplacian)")'); call log(message)
            write(message,'(" >    alpha = ",es12.5)') this%alpha; call log(message)
            write(message,'(" >     beta = ",es12.5)') this%beta; call log(message)
         end if
         write(message,'(" >   it/max = ",i0,"/",i0)') this%niter,this%max_iter; call log(message)
         write(message,'(" >      res = ",es12.5)') this%res; call log(message)
         write(message,'(" >  tol_rel = ",es12.5)') this%tol_rel; call log(message)
      end if
   end subroutine log_solver

   !> Print solver info to screen
   subroutine print_long(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      class(amrmg), intent(in) :: this
      if (this%amr%amRoot) then
         write(output_unit,'("AMR Multigrid solver for grid [",a,"]")') trim(this%amr%name)
         if (this%type .eq. amrmg_cstcoef) then
            write(output_unit,'(" >     type = constant-coefficient (Poisson)")')
         else if (this%type .eq. amrmg_varcoef) then
            write(output_unit,'(" >     type = variable-coefficient (ABecLaplacian)")')
            write(output_unit,'(" >    alpha = ",es12.5)') this%alpha
            write(output_unit,'(" >     beta = ",es12.5)') this%beta
         end if
         write(output_unit,'(" >   it/max = ",i0,"/",i0)') this%niter,this%max_iter
         write(output_unit,'(" >      res = ",es12.5)') this%res
         write(output_unit,'(" >  tol_rel = ",es12.5)') this%tol_rel
      end if
   end subroutine print_long

   !> Short print of solver info
   subroutine print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      class(amrmg), intent(in) :: this
      if (this%amr%amRoot) then
         write(output_unit,'("AMR MG solver [",a16,"] -> it/max = ",i3,"/",i3," res = ",es12.5)') &
            trim(this%amr%name),this%niter,this%max_iter,this%res
      end if
   end subroutine print_short

   !> Finalize and release resources
   subroutine finalize(this)
      class(amrmg), intent(inout) :: this
      if (this%setup_done) call this%destroy()
      nullify(this%amr)
      this%type = -1
   end subroutine finalize

end module amrmg_class
