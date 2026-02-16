!> Test new irlfab
module mod_test_irlfab
   use precision, only: WP
   use messager, only: log, die
   use string, only: rtoa
   implicit none
   private
   public :: test_irlfab

contains

   !> Main test routine for geometry functions
   subroutine test_irlfab()
      use amrvof_geometry, only: cut_tet_vol, get_plane_dist
      real(WP) :: tol, max_err
      
      call log("=== Testing amrvof_irlfab ===")
      
      tol = 1.0e-12_WP
      max_err = 0.0_WP
      
      ! Test 1: cut_tet_vol with plane completely above tet (full tet below plane)
      test_full_tet: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plane
         real(WP) :: vol_liq, vol_gas, vol_exact, err
         real(WP), dimension(3) :: bary_liq, bary_gas
         
         ! Unit tet at origin: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         ! Plane above: n=(0,0,1), d=2 (z=2, all vertices have z<2 so all liquid)
         plane = [0.0_WP, 0.0_WP, 1.0_WP, 2.0_WP]
         
         call cut_tet_vol(tet, plane, vol_liq, vol_gas, bary_liq, bary_gas)
         vol_exact = 1.0_WP / 6.0_WP  ! Volume of unit tet
         err = abs(vol_liq - vol_exact)
         max_err = max(max_err, err)
         
         call log("  Test 1 (full tet): vol_liq = "//rtoa(vol_liq)//", err = "//rtoa(err))
         if (err .gt. tol) call die("Test 1 failed: cut_tet_vol full tet")
      end block test_full_tet
      
      call log("=== All irlfab tests passed! ===")
      
   end subroutine test_irlfab

end module mod_test_irlfab
