!> Simulation hooks for AMR Tester
module simulation
   ! use mod_test_registry,      only: test_registry        ! DISABLED
   use mod_test_amrdata,       only: test_amrdata
   ! use mod_test_visualization, only: test_visualization   ! DISABLED
   use mod_test_amrscalar,     only: test_amrscalar
   ! use mod_test_amrmg,         only: test_amrmg           ! DISABLED
   ! use mod_test_amrabeclap,    only: test_amrabeclap      ! DISABLED
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
      ! call test_visualization()
      call test_amrscalar()
      ! call test_amrmg()
      ! call test_amrabeclap()
   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      call log("Simulation Final: Done")
   end subroutine simulation_final

end module simulation
