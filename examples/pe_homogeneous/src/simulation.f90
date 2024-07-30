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
   real(WP) :: m_0           ! initial mass of the polymer
   real(WP) :: vol_0         ! initial volume of the polymer

   real(WP) :: mass ! dynamic mass of the polymer
   real(WP) :: vol ! dynamic volume of the polymer

   ! Time infomation
   real(WP) :: t_start         ! start time
   real(WP) :: t_end                    ! end time
   real(WP) :: dt                      ! time step
   real(WP) :: tcur               ! current time
   real(WP) :: tret(1)           ! return time

   ! working variables
   real(WP) :: Tloc                                    ! local temperature
   real(WP), dimension(nspec) :: c, c_old                     ! Concentration array
   real(WP), dimension(nreac + nreac_reverse) :: k     ! reaction kinetics frequency array
   real(WP), dimension(nreac + nreac_reverse) :: w     ! reaction rate array
   real(WP), dimension(nspec) :: rhsp                   ! right hand side of the ODE
   real(WP), dimension(nspec,nspec) :: Jacp               ! Jacobian matrix
   character(len=str_medium),dimension(nspec) :: names                   ! species name
   logical, dimension(nspec) :: is_real, is_liquid
   real(WP), dimension(nspec) :: species_mass            ! species mass


   logical :: use_fixed_temperature

   ! CVode variables
   type(c_ptr)                    :: ctx          ! SUNDIALS context
   type(c_ptr)                    :: cvode_mem    ! CVODE memory
   type(N_Vector),        pointer :: sunvec_y     ! sundials vector
   type(SUNMatrix),       pointer :: sunmat_A     ! sundials matrix
   type(SUNLinearSolver), pointer :: sunlinsol_LS ! sundials linear solver
   real(c_double)                 :: rtol=1.0d-7, atol=1.0d-14   ! relative and absolute tolerance
   integer(c_int)                 :: ierr         ! error flag from C functions
   integer(c_long)                :: mxstep = 1000000000      ! maximum number of steps

   ! Needed but not used variables
   real(WP), dimension(nTB + nFO) :: M
   real(WP) :: Ploc
   real(WP), dimension(nqss) :: cqss

   ! Post-processing variables
   real(WP), dimension(nspec) :: mass_fractions
   real(WP) :: total_mass_ratio, LO_ratio, liquid_species_mass_fraction
   real(WP) :: rhsp_mass_sum


   public :: simulation_init,simulation_run,simulation_final

   type(monitor) :: mfile,fracfile

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
            Jmat(i+(j-1)*nspec) = Jacp(i,j)
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

   subroutine check_rhs_mass_balance()
      implicit none
      integer :: i
      real(WP) :: sum

      call get_reaction_rates(w,k,m,c,cqss)
      call fill_rhs_matrix(rhsp, w)

      sum = 0.0_WP
      do i=1,nspec
         sum = sum + rhsp(i)*W_sp(i)*vol_0
      end do

      rhsp_mass_sum = sum

   end subroutine check_rhs_mass_balance

   subroutine update_liquid_concentration()
      implicit none
      integer :: i
      real(WP) :: g

      ! update mass of all species
      do i=1,nspec
         species_mass(i) = species_mass(i) + (c(i)-c_old(i))*W_sp(i)*vol
      end do
      total_mass_ratio = sum(species_mass)/m_0

      ! species_mass = species_mass/total_mass_ratio

      ! calculate the mass and volume of liquid species
      mass = 0.0_WP
      do i=1,nspec
         if (is_liquid(i)) then
            mass = mass + species_mass(i)
         end if
      end do
      vol = mass/rho_0

      ! remove gas species from the concentration and update the concentration of liquid species
      do i=1,nspec
         if (.not. is_liquid(i)) then
            c(i) = 0.0_WP
         else
            c(i) = (species_mass(i)/W_sp(i))/vol
         end if
      end do

      ! calculate the mass fractions respect to initial mass
      do i=1,nspec
         mass_fractions(i) = species_mass(i)/m_0
      end do

      ! Calculate LO ratio (Light Olefins)
      LO_ratio = 0.0_WP
      LO_ratio = mass_fractions(sC2H4)+mass_fractions(sC3H6)+mass_fractions(sC4H8X1)

      ! Calculate mass fraction of liquid species
      liquid_species_mass_fraction = 0.0_WP
      do i=1,nspec
         if (is_liquid(i)) then
            liquid_species_mass_fraction = liquid_species_mass_fraction + mass_fractions(i)
         end if
      end do

   end subroutine update_liquid_concentration

   subroutine simulation_init
      implicit none
      integer :: i

      call isRealSpecies(is_real)
      call isLiquidSpecies(is_liquid)

      call param_read("start time", t_start)
      call param_read("end time", t_end)
      call param_read("time step", dt)

      call param_read("use fixed temperature", use_fixed_temperature)
      if (use_fixed_temperature) then
         call param_read("temperature", Tloc)
      end if

      init_concentration: block
         real(WP) :: chain_number, Nec, Nmc

         integer :: l

         ! read polymer information
         call param_read('Density', rho_0)
         call param_read('Particle diameter', diameter_0)
         call param_read('Number Average Molecular weight', MW_n_0)
         call param_read('Monomer Molecular weight', mMW_0)
         call param_read('Degree of polymerization', DP_0)

         ! Calcute volume
         vol_0 = 4.0_WP/3.0_WP*3.1415926_WP*(diameter_0/2.0_WP)**3
         m_0 = rho_0*vol_0


         ! Calculate chain number
         chain_number = vol_0*rho_0/(MW_n_0/1000.0_WP)     ! in mol

         l = 4
         Nec = 0.0_WP
         do while (l < 5*DP_0)
            Nec = Nec + schultz_dis(l,DP_0)
            l = l + 1
         end do

         l = 4
         Nmc = 0.0_WP
         do while (l < 5*DP_0)
            Nmc = Nmc + schultz_dis(l,DP_0)*(mMW_0*l)
            l = l + 1
         end do

         Nmc = Nmc-Nec*mMW_0*4.0_WP       ! 2 end-groups in each chain, each end-group contains 2 monomers
         Nmc = Nmc/(mMW_0*2.0_WP)         ! each MC contains 2 monomers

         print *, 'Nec: ', Nec
         print *, 'Nmc: ', Nmc

         ! till now, Nec and Nmc are the number fraction of end-chains and MCs, respect to chain number
         ! True MC and EC number need to be multiplied by chain number
         Nec = Nec*chain_number        ! in mol
         Nmc = Nmc*chain_number        ! in mol

         ! Normalize the concentration to make mass balance
         ! Fix the mole
         Nec = 1.0_WP/((W_sp(sPXC16H15GLG)+W_sp(sPXC16H17GLG))/m_0+W_sp(sPXC16H16XPGLG)/m_0*(Nmc/Nec))
         Nmc = (1.0_WP-Nec*(W_sp(sPXC16H15GLG)+W_sp(sPXC16H17GLG))/m_0)/(W_sp(sPXC16H16XPGLG)/m_0)

         ! Initialize concentration
         c(sPXC16H16XPGLG) = Nmc/vol_0
         c(sPXC16H17GLG) = Nec/vol_0
         c(sPXC16H15GLG) = Nec/vol_0

         species_mass = 0.0_WP
         mass = m_0
         vol = vol_0

         call update_liquid_concentration()

      end block init_concentration

      init_kinetics: block
         k = 0.0_WP
         w = 0.0_WP
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
         ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn))
         if (ierr /= 0) then
            print *, 'ERROR: FCVodeInit failed'
            stop 1
         end if

         ! set max step size
         ierr = FCVodeSetMaxStep(cvode_mem, dt)

         ! set max number of iterations in one step
         ierr = FCVodeSetMaxNumSteps(cvode_mem, mxstep)

      end block init_CVODE

      ! Get species names
      call get_species_names(names)

      ! Create Monitor file
      mfile=monitor(.true.,'concentration')
      call mfile%add_column(tcur, 'Time')
      do i=1,nspec
         call mfile%add_column(c(i), names(i))
      end do
      call mfile%write()

      fracfile=monitor(.true.,'mass_fraction')
      call fracfile%add_column(tcur, 'Time')
      call fracfile%add_column(total_mass_ratio, 'sum of mass fraction')
      call fracfile%add_column(rhsp_mass_sum, 'sum of residure')
      call fracfile%add_column(Tloc, 'temperature')
      call fracfile%add_column(liquid_species_mass_fraction, 'liquid mass fraction')
      call fracfile%add_column(LO_ratio, 'LO ratio')
      call fracfile%add_column(mass_fractions(sC6H6), 'Benzene')
      call fracfile%add_column(mass_fractions(sC7H8), 'C7H8')
      call fracfile%add_column(mass_fractions(sC8H10), 'C8H10')
      call fracfile%add_column(mass_fractions(sC8H8), 'C8H8')
      call fracfile%add_column(mass_fractions(sC9H10), 'C9H10')
      call fracfile%write()

   end subroutine simulation_init



   subroutine simulation_run
      implicit none

      integer :: i,count

      count = -1
      tcur = t_start

      ! Time loop
      do while (tcur < t_end)
         print *, '[Solving] Time: ', tcur, 's'

         ! Update Temperature
         if (.not. use_fixed_temperature) then
            call get_temperature(T = Tloc, time = tcur)
         end if
         print *, 'Temperature: ', Tloc

         ! Update kinetics coefficients
         call get_rate_coefficients(k,M,Tloc,Ploc)

         ! store the old concentration
         c_old = c
         ! Update CVODE
         ierr = FCVode(cvode_mem, (tcur+dt), sunvec_y, tret, CV_NORMAL)
         if (ierr /= 0) then
            print *, 'ERROR: FCVode failed'
            stop 1
         end if
         print *, 'Returned Time: ', tret(1), 's'

         ! Update concentration
         call update_liquid_concentration()

         ! Update time
         tcur = tcur + dt

         ! reinit CVODE
         ! free old sunvecy and linear solver
         ierr = FSUNLinSolFree(sunlinsol_LS)
         call FN_VDestroy(sunvec_y)
         sunvec_y => FN_VMake_Serial(nspec, c, ctx)
         ! re-create linear solver
         sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, ctx)
         ! re-set linear solver
         ierr = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)
         ! re-init cvode
         ierr = FCVodeReInit(cvode_mem, tcur, sunvec_y)

         call check_rhs_mass_balance()

         if (mod(count,100) == 0) then
            ! Write to monitor file
            call mfile%write()
            call fracfile%write()
         end if
         count = count + 1

      end do


   end subroutine simulation_run


   subroutine simulation_final

   end subroutine simulation_final



end module simulation
