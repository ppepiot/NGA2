module simulation
   use precision,         only: WP
   use, intrinsic :: iso_c_binding
   use fsundials_context_mod
   use fcvode_mod                 ! Fortran interface to CVODE
   use fnvector_serial_mod        ! Fortran interface to serial N_Vector
   use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
   use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
   use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
   use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
   use fsundials_nvector_mod      ! Fortran interface to generic N_Vector
   use mod_pyrolysis
   use temperature_mod
   use monitor_class,     only: monitor
   use string, only: str_medium
   use param, only: param_read

   implicit none
   private
   ! Polymer infomation
   real(WP) :: rho_0           ! initial polymer density
   real(WP) :: diameter_0      ! diameter of the polymer
   real(WP) :: mMW_0          ! initial monomer molecular weight of the polymer
   real(WP) :: MW_n_0          ! initial number average molecular weight of the polymer
   real(WP) :: MW_w_0          ! initial weight average molecular weight of the polymer
   real(WP) :: DP_0         ! initial degree of polymerization


   ! Time infomation
   real(WP) :: t_start = 0.0_WP         ! start time
   real(WP) :: t_end = 20732.0_WP                    ! end time
   real(WP) :: dt = 1e-1_WP                      ! time step
   real(WP) :: tcur               ! current time
   real(WP) :: tret(1)           ! return time

   ! working variables
   real(WP) :: Tloc                                    ! local temperature
   real(WP), dimension(nspec) :: c                     ! Concentration array
   real(WP), dimension(nreac + nreac_reverse) :: k     ! reaction kinetics frequency array
   real(WP), dimension(nreac + nreac_reverse) :: w     ! reaction rate array
   real(WP), dimension(nspec) :: rhsp                   ! right hand side of the ODE
   real(WP), dimension(nspec,nspec) :: Jacp               ! Jacobian matrix
   character(len=str_medium),dimension(nspec) :: names                   ! species name


   ! CVode variables
   type(c_ptr)                    :: ctx          ! SUNDIALS context
   type(c_ptr)                    :: cvode_mem    ! CVODE memory
   type(N_Vector),        pointer :: sunvec_y     ! sundials vector
   type(SUNMatrix),       pointer :: sunmat_A     ! sundials matrix
   type(SUNLinearSolver), pointer :: sunlinsol_LS ! sundials linear solver
   real(c_double)                 :: rtol=1.0d-5, atol=1.0d-10   ! relative and absolute tolerance
   integer(c_int)                 :: ierr         ! error flag from C functions

   ! Needed but not used variables
   real(WP), dimension(nTB + nFO) :: M
   real(WP) :: Ploc
   real(WP), dimension(nqss) :: cqss

   public :: simulation_init,simulation_run,simulation_final

   type(monitor) :: mfile

contains

   ! Define the rhs function and Jac function
   ! ----------------------------------------------------------------
   ! RhsFn provides the right hand side function for the
   ! ODE: dy/dt = f(t,y)
   !
   ! Return values:
   !    0 = success,
   !    1 = recoverable error,
   !   -1 = non-recoverable error
   ! ----------------------------------------------------------------
   integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
      result(ierr) bind(C,name='RhsFn')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod

      !======= Declarations =========
      implicit none

      ! calling variables
      real(c_double), value :: tn        ! current time
      type(N_Vector)        :: sunvec_y  ! solution N_Vector
      type(N_Vector)        :: sunvec_f  ! rhs N_Vector
      type(c_ptr),    value :: user_data ! user-defined data

      ! pointers to data in SUNDIALS vectors
      real(c_double), pointer :: yvec(:)
      real(c_double), pointer :: fvec(:)

      integer :: i

      ! get data arrays from SUNDIALS vectors
      yvec => FN_VGetArrayPointer(sunvec_y)
      if (.not. associated(yvec)) then
         write(*,*) 'ERROR: yvec = NULL'
      end if
      fvec => FN_VGetArrayPointer(sunvec_f)
      if (.not. associated(fvec)) then
         write(*,*) 'ERROR: fvec = NULL'
      end if

      call get_reaction_rates(w,k,m,yvec,cqss)
      call fill_rhs_matrix(rhsp, w)

      do i = 1, nspec
         fvec(i) = rhsp(i)
      end do

      ! return success
      ierr = 0
      return

   end function RhsFn

   ! ----------------------------------------------------------------
   ! JacFn: The Jacobian of the ODE hand side function J = df/dy
   !
   ! Return values:
   !    0 = success,
   !    1 = recoverable error,
   !   -1 = non-recoverable error
   ! ----------------------------------------------------------------
   integer(c_int) function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, &
      user_data, tmp1, tmp2, tmp3) &
      result(ierr) bind(C,name='JacFn')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod
      use fsunmatrix_dense_mod
      use fsundials_matrix_mod

      !======= Declarations =========
      implicit none

      integer :: i, j

      ! calling variables
      real(c_double), value :: tn               ! current time
      type(N_Vector)        :: sunvec_y         ! current solution N_Vector
      type(N_Vector)        :: sunvec_f         ! current rhs N_Vector
      type(SUNMatrix)       :: sunmat_J         ! Jacobian SUNMatrix
      type(c_ptr), value    :: user_data        ! user-defined data
      type(N_Vector)        :: tmp1, tmp2, tmp3 ! workspace N_Vectors

      ! pointers to data in SUNDIALS vectors
      real(c_double), pointer :: yvec(:)

      ! pointer to data in SUNDIALS matrix
      real(c_double), pointer :: Jmat(:)

      ! get data array from SUNDIALS vector
      yvec => FN_VGetArrayPointer(sunvec_y)
      if (.not. associated(yvec)) then
         write(*,*) 'ERROR: yvec = NULL'
      end if
      ! get data arrays from SUNDIALS vectors
      Jmat => FSUNDenseMatrix_Data(sunmat_J)
      if (.not. associated(Jmat)) then
         write(*,*) 'ERROR: Jmat = NULL'
      end if

      call fill_jac_matrix(Jacp, yvec, k)

      ! Fill the Jacobian matrix
      do i=1, nspec
         do j=1, nspec
            Jmat(i+(j-1)*nspec) = Jacp(j,i)
         end do
      end do


      ! return success
      ierr = 0
      return

   end function JacFn


   function schultz_dis(l,DP) result(res)
      implicit none
      integer, intent(in) :: l
      real(WP), intent(in) :: DP
      real(WP) :: res

      res = 4.0_WP*l/(DP-1.0_WP)**2*((DP-1.0_WP)/(DP+1.0_WP))**l

   end function schultz_dis


   subroutine simulation_init
      implicit none
      integer :: i


      init_concentration: block
         real(WP) :: vol, chain_number, Nec, Nmc

         integer :: l

         ! read polymer information
         call param_read('Density', rho_0)
         call param_read('Particle diameter', diameter_0)
         call param_read('Number Average Molecular weight', MW_n_0)
         call param_read('Monomer Molecular weight', mMW_0)
         call param_read('Degree of polymerization', DP_0)

         ! Calcute volume
         vol = 4.0_WP/3.0_WP*3.1415926_WP*(diameter_0/2.0_WP)**3

         ! Calculate chain number
         chain_number = vol*rho_0/(MW_n_0/1000.0_WP)     ! in mol

         l = 4
         Nec = 0.0_WP
         do while (l < 5*DP_0)
            Nec = Nec + schultz_dis(l,DP_0)
            l = l + 1
         end do
         print *, 'Nec = ', Nec

         l = 4
         Nmc = 0.0_WP
         do while (l < 5*DP_0)
            Nmc = Nmc + schultz_dis(l,DP_0)*(mMW_0*l)
            l = l + 1
         end do
         print *, 'Nmc = ', Nmc

         Nmc = Nmc-Nec*mMW_0*4.0_WP       ! 2 end-groups in each chain, each end-group contains 2 monomers
         Nmc = Nmc/(mMW_0*2.0_WP)         ! each MC contains 2 monomers

         ! till now, Nec and Nmc are the number fraction of end-chains and MCs, respect to chain number
         ! True MC and EC number need to be multiplied by chain number
         Nec = Nec*chain_number
         Nmc = Nmc*chain_number

         ! Initialize concentration
         c(sPXC16H16XPGLG) = Nmc/vol
         c(sPXC16H17GLG) = Nec/vol
         c(sPXC16H15GLG) = Nec/vol


      end block init_concentration


      init_kinetics: block
         k = 0.0_WP
         w = 0.0_WP
         Tloc = 0.0_WP
      end block init_kinetics


      init_CVODE: block

         ierr = FSUNContext_Create(c_null_ptr, ctx)

         ! create SUNDIALS N_Vector with initial values
         sunvec_y => FN_VMake_Serial(nspec, c, ctx)
         if (.not. associated(sunvec_y)) then
            print *, 'ERROR: sunvec = NULL'
            stop 1
         end if

         sunmat_A => FSUNDenseMatrix(nspec, nspec, ctx)
         if (.not. associated(sunmat_A)) then
            print *, 'ERROR: sunmat = NULL'
            stop 1
         end if

         ! create a linear solver
         sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, ctx)
         if (.not. associated(sunlinsol_LS)) then
            print *, 'ERROR: sunlinsol = NULL'
            stop 1
         end if

         ! create CVode memory
         cvode_mem = FCVodeCreate(CV_BDF, ctx)
         if (.not. c_associated(cvode_mem)) then
            print *, 'ERROR: cvode_mem = NULL'
            stop 1
         end if

         ! initialize CVode
         ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), t_start, sunvec_y)
         if (ierr /= 0) then
            print *, 'ERROR: FCVodeInit failed'
            stop 1
         end if

         ! set relative and absolute tolerances
         ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
         if (ierr /= 0) then
            print *, 'ERROR: FCVodeInit failed'
            stop 1
         end if

         ! attach linear solver
         ierr = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)
         if (ierr /= 0) then
            print *, 'ERROR: FCVodeInit failed'
            stop 1
         end if

         ! set Jacobian routine
         ! ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn))
         ! if (ierr /= 0) then
         !    print *, 'ERROR: FCVodeInit failed'
         !    stop 1
         ! end if


      end block init_CVODE

      ! Get species names
      call get_species_names(names)

      ! Create Monitor file
      mfile=monitor(.true.,'concentration')
      do i=1,nspec
         call mfile%add_column(c(i), names(i))
      end do
      call mfile%write()

   end subroutine simulation_init



   subroutine simulation_run
      implicit none

      tcur = t_start

      ! Time loop
      do while (tcur < t_end)
         print *, '[Solving] Time: ', tcur, 's'
         ! Update Temperature
         call get_temperature(T = Tloc, time = tcur)
         print *, 'Temperature: ', Tloc

         ! Update kinetics coefficients
         call get_rate_coefficients(k,M,Tloc,Ploc)

         ! Update CVODE
         ierr = FCVode(cvode_mem, (tcur+dt), sunvec_y, tret, CV_NORMAL)
         if (ierr /= 0) then
            print *, 'ERROR: FCVode failed'
            stop 1
         end if

         ! Write to monitor file
         call mfile%write()

         ! Update time
         tcur = tcur + dt
      end do


   end subroutine simulation_run


   subroutine simulation_final

   end subroutine simulation_final



end module simulation
