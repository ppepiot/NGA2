!> Simulation hooks for AMR Tester
module simulation
   use mod_test_amrdata,       only: test_amrdata
   use mod_test_amrscalar,     only: test_amrscalar
   use mod_test_amrabeclap,    only: test_amrabeclap
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
      call test_amrdata()
      call test_amrscalar()
      call test_amrabeclap()
   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      call log("Simulation Final: Done")
   end subroutine simulation_final

end module simulation
