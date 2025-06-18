!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision, only: WP
   use ceq_types, only: sys_type
   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final
   
   type (sys_type), pointer :: sys 
   integer :: ns, ne, ncs, ng, diag, lu_op, iret, info
   real(WP), dimension(2,2) :: Ein
   real(WP), dimension(2,15) :: thermo_in
   real(WP), dimension(2) :: neq
   integer, dimension(:), allocatable :: CS
   real(WP), dimension(:,:), allocatable :: Bg
   real(WP) :: pressure

   integer :: i, it, nt, n, imin 
   real(WP), dimension(2) :: ss, hs
   real(WP), dimension(:), allocatable :: NL,g
   real(WP) :: N0, Nend, T0, Tend, T, gmin, gL, gG

contains

   !> Nothing to initialize
   subroutine simulation_init
      implicit none
      
      ! Number of species: 1=H2Ol, 2=H2Ov, 3=N2
      ns = 3  
      ! Number of elements: 1=H, 2=O, 3=N
      ne = 3
      ! Number of constrained species
      ncs = 0 
      ! Number of general linear constraints
      ng = 0 
      ! Ein - element matrix (ns x ne).  A molecule of species k contains Ein(k,j) atoms of element j.
      Ein (1,1) = 2
      Ein (1,2) = 1
      Ein (1,3) = 0
      Ein (2,1) = 2
      Ein (2,2) = 1
      Ein (2,3) = 0
      Ein (3,1) = 0
      Ein (3,2) = 0
      Ein (3,3) = 2

      ! Array containing indexes of constrained species (n_cs x 1)
      allocate(CS(ncs))
      ! General linear constraint matrix (ns x ng) - Unconstrained calculations
      allocate(Bg(ns,ng))

      ! Thermodynamic data for species (ns x 15)
      ! Liquid water
      thermo_in(1,1) = 1000_WP ! Tmid
      thermo_in(1,2) = 72.5575005_WP ! a1 - low T range
      thermo_in(1,3) = -0.662445402_WP ! a2 - low T range
      thermo_in(1,4) = 0.00256198746_WP ! a3 - low T range
      thermo_in(1,5) = -4.36591923e-6_WP ! a4 - low T range
      thermo_in(1,6) = 2.78178981e-9_WP ! a5 - low T range
      thermo_in(1,7) = -41886.5499_WP ! a6 - low T range
      thermo_in(1,8) = -288.280137_WP ! a7 - low T range
      thermo_in(1,9) = 0_WP ! A1 - high T range
      thermo_in(1,10) = 0_WP ! A2 - high T range
      thermo_in(1,11) = 0_WP ! A3 - high T range
      thermo_in(1,12) = 0_WP ! A4 - high T range
      thermo_in(1,13) = 0_WP ! A5 - high T range
      thermo_in(1,14) = 0_WP ! A6 - high T range
      thermo_in(1,15) = 0_WP ! A7 - high T range

      ! Water - gas
      thermo_in(2,1) = 1000_WP ! Tmid
      thermo_in(2,2) = 4.1986352_WP ! a1 - low T range
      thermo_in(2,3) = -0.0020364017_WP ! a2 - low T range
      thermo_in(2,4) = 6.5203416e-6_WP ! a3 - low T range
      thermo_in(2,5) = -5.4879269e-9_WP ! a4 - low T range
      thermo_in(2,6) = 1.771968e-12_WP ! a5 - low T range
      thermo_in(2,7) = -30293.726_WP ! a6 - low T range
      thermo_in(2,8) = -0.84900901_WP ! a7 - low T range
      thermo_in(2,9) = 2.6770389_WP ! A1 - high T range
      thermo_in(2,10) = 0.0029731816_WP ! A2 - high T range
      thermo_in(2,11) = -7.7376889e-7_WP ! A3 - high T range
      thermo_in(2,12) = 9.4433514e-11_WP ! A4 - high T range
      thermo_in(2,13) = -4.2689991e-15_WP ! A5 - high T range
      thermo_in(2,14) = -29885.894_WP ! A6 - high T range
      thermo_in(2,15) = 6.88255_WP ! A7 - high T range

      ! N2 - gas
      thermo_in(3,1) = 1000_WP ! Tmid
      thermo_in(3,2) = 3.53100528_WP ! a1 - low T range
      thermo_in(3,3) = -1.23660988e-4_WP ! a2 - low T range
      thermo_in(3,4) = -5.02999433e-7_WP ! a3 - low T range
      thermo_in(3,5) = 2.43530612e-9_WP ! a4 - low T range
      thermo_in(3,6) = -1.40881235e-12_WP ! a5 - low T range
      thermo_in(3,7) = -1046.97628_WP ! a6 - low T range
      thermo_in(3,8) = 2.96747038_WP ! a7 - low T range
      thermo_in(3,9) = 2.95257637_WP ! A1 - high T range
      thermo_in(3,10) = 1.3969004e-3_WP ! A2 - high T range
      thermo_in(3,11) = -4.92631603e-7_WP ! A3 - high T range
      thermo_in(3,12) = 7.86010195e-11_WP ! A4 - high T range
      thermo_in(3,13) = -4.60755204e-15_WP ! A5 - high T range
      thermo_in(3,14) = -923.948688_WP ! A6 - high T range
      thermo_in(3,15) = 5.87188762_WP ! A7 - high T range


   end subroutine simulation_init

   !> CEQ testing
   subroutine simulation_run
      implicit none

      ! Ideal gas constant J/mol.K
      R = 8.3145
      ! Pressure [Pa]
      P0 = 1e5_WP
      pref = 1e5_WP
      ! Initial volume [m^3]
      V0 = 1e-6
      ! Initial temperatures [K]
      Tl0 = 300.0_WP
      Tg0 = 320.0_WP
      ! Initial liquid and gas fraction
      alphal = 0.5
      alphag = 0.5
      ! Initial mole fractions in the gas phase (single component in liquid phase)
      XH20 = 0.5
      XN2 = 0.5
      ! Liquid water density (assume constant) [kg/m3]
      rhoH20l = 1e3
      ! Molar masses
      MM(1) = 0.016
      MM(2) = 0.016
      MM(3) = 0.028
      ! Calculate number of moles
      N0(1) = alphal*V0*rhoH20l/MM(1) ! Liquid water
      N0(2) = XH20*P0*alphag*V0/(R*T0)! Water vapor
      N0(3) = ! N2
 




      ! Temperature array
      T0 = 300
      Tend = 400
      nt = 100
      ! Moles of liquid
      N0 = 0
      Nend = 1
      n = 100
      allocate(NL(n),g(n))
      do i=1,n
         NL(i) = N0 + (i-1)*(Nend-N0)/(n-1)
      end do

      ! Go through temperatures
      do it = 1,nt
         T = T0 + (it-1)*(Tend-N0)/(nt-1)
         ! Calculate for each liquid/gas mole numbers the normalized free gibbs energy
         g = 0.0_WP
         do i=1,n
            if (i.eq.1) NL(i) = tiny(1.0_WP)
            if (i.eq.n) NL(i) = 1.0_WP-tiny(1.0_WP)
            call ceq_s(ns,T,pressure,thermo_in,ss(:))
            call ceq_h(ns,T,thermo_in,hs(:))
            gL = NL(i)*(hs(1)-ss(1)) ! Assume constant volume to start with
            gG = (1.0_WP-NL(i))*(hs(2)-ss(2)+log(pressure)) ! Mole fraction in gas is 1
            g(i) = gL + gG 
         end do
         ! Find location of minimum g
         imin = 1
         gmin = g(1)
         do i=2,n
            if (g(i)<gmin) then
               gmin = g(i) 
               imin = i
            end if
         end do
         ! Print solution
         print*,imin,T,NL(imin),gmin
   
      end do

      ! initialize SYS
      !call ceq_sys_init(ns, ne, ncs, ng, Ein, CS, Bg, thermo_in, lu_op, diag, sys, iret)

      ! equilibrium calculation for SYS
      !myp = 0.01_WP!0.57868_WP
      !myT = 350.0_WP
      !call ceq_state(sys, N = [0.5_WP,0.5_WP], p_atm=myp, T=myT, N_eq=neq, info=info )

      !print*, info
      !print*, neq

      ! 300K H = -241 762 J/mol = -15 110 125 J/kg
      !Pelanti: 2 476 160
      
   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      
   end subroutine simulation_final


   
   
end module simulation
