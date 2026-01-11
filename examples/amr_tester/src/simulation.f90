!> Simulation hooks for AMR Tester
module simulation
   use mod_test_registry,      only: test_registry
   use mod_test_visualization, only: test_visualization
   use mod_test_amrscalar,     only: test_amrscalar
   use mod_test_multidata,     only: test_multidata
   use messager, only: log
   implicit none

contains

   !> Initialization hook
   subroutine simulation_init()
      call log("Simulation Init: AMR Tester")
   end subroutine simulation_init

   !> Clean run hook
   subroutine simulation_run()
      call log("Simulation Run: Executing Tests")
      call test_multidata()  ! Test current_amrdata robustness first
      call test_registry()
      call test_visualization()
      call test_amrscalar()
   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      call log("Simulation Final: Done")
   end subroutine simulation_final

end module simulation
