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
   use param, only: param_read

   implicit none
   private


contains

    subroutine simulation_init

    end subroutine simulation_init



    subroutine simulation_run

    end subroutine simulation_run


    subroutine simulation_final

    end subroutine simulation_final



end module simulation