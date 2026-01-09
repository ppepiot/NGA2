module mod_test_registry
   use iso_c_binding
   use amrgrid_class
   use amrdata_class
   use messager, only: log,warn
   use precision, only: WP
   implicit none
   private

   public :: test_registry

contains

   subroutine test_registry()
      type(amrgrid) :: amr
      type(amrdata), target :: velocity
      logical :: success

      call log("---------------------------------------------------")
      call log("Running Test: Registry Pattern (Auto-Allocation)")
      call log("---------------------------------------------------")

      ! 1. Initialize AMR Config
      ! We need to set parameters manually as we are bypassing command line here
      amr%nx=32
      amr%ny=32
      amr%nz=32
      amr%xlo=0.0_WP
      amr%xhi=1.0_WP
      amr%ylo=0.0_WP
      amr%yhi=1.0_WP
      amr%zlo=0.0_WP
      amr%zhi=1.0_WP
      amr%nlvl=0 ! Start with single level
      amr%nmax=32

      call amr%initialize("TestGrid")
      call log("AMR Config Initialized.")

      ! 2. Register Data with name, ncomp, ng
      call amr%register(data=velocity,name="velocity",ncomp=3,ng=1)
      call log("Data 'velocity' registered (ncomp=3, ng=1).")

      ! 3. Initialize Grid (Triggers Allocation)
      ! This calls C++ MakeNewLevelFromScratch -> dispatch_mak_lvl_init -> field%define
      call amr%initialize_grid(0.0d0)
      call log("Grid Initialized.")

      ! 4. Verify Allocation
      ! Access level 0 of the data
      if (allocated(velocity%mf)) then
         if (c_associated(velocity%mf(0)%p)) then
            call log("PASS: Data allocated automatically!")
         else
            call warn("FAIL: Data MF pointer is null!")
         end if
      else
         call warn("FAIL: Data MF array not allocated!")
      end if

      ! 4b. Verify BC arrays allocated
      if (allocated(velocity%lo_bc).and.allocated(velocity%hi_bc)) then
         if (size(velocity%lo_bc,2)==3) then
            call log("PASS: BC arrays allocated (3 components)!")
         else
            call warn("FAIL: BC array has wrong size!")
         end if
      else
         call warn("FAIL: BC arrays not allocated!")
      end if


      ! 5. Test Unregister
      call amr%unregister(data=velocity)
      if (allocated(amr%data)) then
         call warn("FAIL: Registry should be empty after unregister!")
      else
         call log("PASS: Unregister successful (Registry empty).")
      end if


      ! Cleanup
      call amr%finalize()

   end subroutine test_registry

end module mod_test_registry
