!> Chemical state class
module chem_state_class
   use precision,      only: WP
   use chem_sys_class, only: chem_sys,Lphase,Gphase,ncof
   implicit none
   private
   
   !> Expose type/constructor/methods
   public :: chem_state

   !> List of available equilibrium conditions
   integer, parameter, public :: fixed_PT=1                     !< Fixed temperature
   integer, parameter, public :: fixed_PH=2                     !< Fixed enthalpy

   !> Temperature bounds
   real(WP), parameter :: T_high=5000.0_WP
   real(WP), parameter :: T_low =250.0_WP

   !> Reference pressure (Pa)
   real(WP), parameter :: p0=1.01325e5

   !> Chemical statete object definition
   type :: chem_state

      ! This is our chemical system
      class(chem_sys), pointer :: sys                          !< This is the chemical system the solver is build for

      ! Equilibrium condition
      integer :: cond

      ! Thermochemical quantities
      real(WP) :: p                                            !< Pressure
      real(WP) :: T                                            !< Temperature
      real(WP) :: HoR                                          !< Enthalpy
      real(WP), dimension(:),   allocatable :: N               !< Moles of species (ns)
      real(WP), dimension(:),   allocatable :: Nbar            !< Moles of phases (np)
      real(WP), dimension(:),   allocatable :: Nd              !< Moles of determined species (nsd)
      real(WP), dimension(:),   allocatable :: Nu              !< Moles of undetermined species (nsu)
      real(WP), dimension(:),   allocatable :: lam             !< Lagrange multipliers (nrc)
      real(WP), dimension(:),   allocatable :: cr              !< Reduced constraint vector
      real(WP), dimension(:),   allocatable :: gu              !< Gibbs functions of undetermined species (nsu)
      real(WP), dimension(:),   allocatable :: R,Rd            !< Residual arrays
      real(WP), dimension(:),   allocatable :: x               !< Chemical state unknowns
      real(WP), dimension(:,:), allocatable :: Btilde,Ptilde   !< Coefficient matrices
      real(WP), dimension(:,:), allocatable :: BtildeT,PtildeT !< Coefficient matrices transposed

      ! Pointer to the avaiable chemical equilibrium procedures
      procedure(get_ceq_interface), pointer :: get_ceq=>NULL() !< Map the input to the latent variables

      ! Numerical parameters
      real(WP) :: tol_N                                        !< Tolerance for the residual norm
      real(WP) :: tol_T                                        !< Tolerance for the temperature
      real(WP) :: dT                                           !< Residual error for the temperature
      integer  :: iter_N                                       !< Number of Newton-Raphson iterations
      integer  :: iter_T                                       !< Number of temperature iterations
      integer  :: iter_N_max                                   !< Maximum number of Newton-Raphson iterations
      integer  :: iter_T_max                                   !< Maximum number of temperature iterations

   contains
      procedure :: initialize                                  !< Object initializer
      procedure :: get_hort                                    !< Get the normalized enthalpy
      procedure :: get_gort                                    !< Get the normalized Gibbs free energy
      procedure :: hor2T                                       !< Convert enthalpy to temperature
      procedure :: get_dgdT                                    !< 
      procedure :: perturb                                     !< Perturb the chemical equilibrium problem
      procedure :: get_Nming                                   !< Get the composition that minimized G and satisfies the constraints
      procedure :: min_pert                                    !< Get the purturbed maxmin composition
      procedure :: get_cpor                                    !< Get the normalized Cp
      procedure :: maxmin_comp                                 !< Get the minmax composition
      procedure :: solve_linprog                               !< Solve the linear programming problem
      procedure :: get_y                                       !< Get the square root of the mole numbers
      procedure :: get_res                                     !< Get the residual vector
      procedure :: x_init                                      !< Initialize the chemical state solution vector
      procedure :: equilibrate                                 !< Obtain the chemical equilibrium state of the system
      procedure :: get_Cp_eff                                  !< Get the effective Cp
      procedure :: get_dxdT                                    !< Get the temperature derivative of the solution vector
      procedure :: get_BP                                      !< Get the coefficient matrices for constraints and phase summation
      procedure, private :: get_ceq_PT                         !< Get the chemical equilibrium state at constant pressure and temperature
      procedure, private :: get_ceq_PH                         !< Get the chemical equilibrium state at constant pressure and emthalpy
   end type chem_state

   !> Interface for get_ceq
   interface
      subroutine get_ceq_interface(this,Neq)
         use precision, only: WP
         import chem_state
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%ns), intent(inout) :: Neq
      end subroutine get_ceq_interface
   end interface


   contains


      !> Chemical state initializer
      subroutine initialize(this,sys,cond,p,T,c,N,HoR,N_h,T_h,N_g,T_g)

         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation. 

         ! Initializes the constrained equilibrium state of 
         ! an ideal liquid-gas mixture consisting of ns species,either at fixed pressure and
         ! temperature (p,T),or at fixed pressure and enthalpy (p,H).

         ! The nc equality constraints are written:  B'*N=c,where B is the
         ! ns x nc basic constraint matrix,N is the ns-vector of species moles,
         ! and c is the nc-vector of constraint values.

         ! Input:
         !  sys -object type chem_sys created by subroutine.
         !         Required first argument

         !  cond -Equilibrium condition: Either T_fixed or H_fixed

         !  (Either c or N must be specified.)
         !  c   -values of the nc basic constraints (real(nc))
         !  N   -moles of species used to calculate c as : c=B'*N (real(ns))

         !  p   -pressure

         ! (For a fixed (p,T) equilibrium calculation,specify T: do not specify HoR,N_h or T_h.)
         ! T-temperature (K) for fixed-temperature problem

         ! (For a fixed (p,H) equilibrium calculation,specify either HoR or N_h and T_h: 
         !  do not specify T.)
         ! HoR the fixed value of H/R [moles K],where H=enthalpy,R=universal gas constant.
         ! N_h species moles used to calculate H  (real(ns))
         ! T_h temperature used to calculate H as:  H/R=sum(N_h h(T_h)/R),where
         !     h(T_h)/R [which has dimensions K] is the molar specific species enthalpy.

         ! (Initial guesses are not needed, and should not be specified unless they
         !  are good guesses.)
         ! N_g  -initial guess for species moles
         ! T_g  -initial guess for temperature (for fixed (p,H) only)

         !  To diagnose an error condition,the case can be repeated with diagnostics turned on,
         !  by: call param_set(sys,diag=5).

         use, intrinsic :: iso_fortran_env, only: output_unit
         use messager,  only: die
         use mathtools, only: reorder_rows
         implicit none
         class(chem_state), intent(inout) :: this
         class(chem_sys), target, intent(in) :: sys
         integer,  intent(in) :: cond
         real(WP), intent(in) :: p
         real(WP), intent(in), optional :: c(sys%nc),N(sys%ns),T,HoR,N_h(sys%ns),T_h,N_g(sys%ns),T_g
         real(WP), parameter :: frac_Nm =0.1_WP! fraction of zm used in initial guess
         integer  :: np,nb,nc,ns,nsd,nsu,nrc,npert,iret,i
         real(WP) :: max_pert,Numin,cb(sys%nc),cmod(sys%nb),Nd(sys%nsd),cr_norm,cb_norm,res,N_low,res_tol=1e-9
         real(WP), dimension(sys%ns)  :: N0,N1,h,Neq
         real(WP), dimension(sys%nsu) :: Nu,Nu0,Nm,Nupper,Ng,gu
         real(WP), dimension(sys%nrc) :: cr
         logical :: fail,diag

         ! Point to chemical system
         this%sys=>sys

         ! Set the eqiuilibrium condition
         select case (cond)
            case (fixed_PT)
               if (present(T)) then
                  if (T.lt.T_low.or.T.gt.T_high) call die('[chem_sys initialize] Temperature out of range')
                  this%T=T
               else
                  call die('[chem_sys initialize] Temperature is required for the fixed temperature condition')
               end if
               this%get_ceq=>get_ceq_PT
            case (fixed_PH)
               if (present(HoR)) then
                  this%HoR=HoR
               elseif(present(N_h).and.present(T_h)) then
                  call reorder_rows(N_h,sys%sp_order,N0)
                  call this%get_hort(sys%ns,T_h,sys%thermo,h)
                  this%HoR=sum(N0*h)*T_h
               else
                  call die('[chem_sys initialize] Both N_h and T_h are required for the fixed enthalpy case')
               end if
               if (present(T_g)) then
                  ! Guess provided
                  this%T=T_g
                  if (T_g.lt.T_low.or.T_g.gt.T_high) call die('[chem_sys initialize] Guessed temperature out of range')
               else
                  this%T=sqrt(T_low*T_high)
                  this%T=max(this%T,0.1_WP*T_high)
               endif
               this%get_ceq=>get_ceq_PH
            case default
               call die('[chem_sys initialize] The chemical state must be at either constant temperature or constant enthalpy')
         end select
         this%cond=cond

         ! Determine pressure in standard atmospheres
         if (p.le.0.0_WP) call die('[chem_sys initialize] Pressure must be strictly positive')
         this%p=p/p0

         ! Obtain indexes
         np =sys%np
         nb =sys%nb
         nc =sys%nc
         nrc=sys%nrc
         ns =sys%ns
         nsd=sys%nsd
         nsu=sys%nsu

         ! Allocate arrays
         allocate(this%N   (ns));           this%N      =0.0_WP
         allocate(this%Nbar(np));           this%Nbar   =0.0_WP
         allocate(this%Nd  (nsd));          this%Nd     =0.0_WP
         allocate(this%Nu  (nsu));          this%Nu     =0.0_WP
         allocate(this%lam (nrc));          this%lam    =0.0_WP
         allocate(this%cr  (nrc));          this%cr     =0.0_WP
         allocate(this%gu  (nsu));          this%gu     =0.0_WP
         allocate(this%R   (nrc+np));       this%R      =0.0_WP
         allocate(this%Rd  (np));           this%Rd     =0.0_WP
         allocate(this%x   (nrc+np));       this%x      =0.0_WP
         allocate(this%Btilde(nsu,nrc));    this%Btilde =0.0_WP
         allocate(this%Ptilde(nsu,np));     this%Ptilde =0.0_WP
         allocate(this%BtildeT(nrc,nsu));   this%BtildeT=0.0_WP
         allocate(this%PtildeT(np,nsu));    this%PtildeT=0.0_WP

         ! Initialize the Gibbs function
         call this%get_gort(nsu,this%T,p,sys%thermo(nsd+1:ns,:),sys%P(nsd+1:ns,Gphase),gu)

         ! Form the basic constraint vector
         if (present(c)) then
            cb(1:nc)=c
         elseif(present(N)) then
            cb(1:nc)=matmul(N,sys%B)
         else
            call die('[chem_sys initialize] Neither c nor N specified')
         endif

         ! Form modified and reduced constraints
         cmod(1:nb)=matmul(sys%A(1:nb,1:nc),cb)
         Nd(1:nsd) =cmod(1:nsd)

         ! Treat the special case of no undetermined species
         if (nsu.eq.0) then
            Neq=Nd
            if (cond.eq.fixed_PH) call this%hor2T(ns,Neq,this%HoR,sys%thermo,this%T)
            go to 500
         endif
         
         ! Reduced constraint vector
         cr(1:nrc) =cmod(nsd+1:nb)
         cr_norm   =norm2(cr)

         if (cr_norm.le.0.0_WP) then
         ! SBP added 4/9/2009
            if (cr_norm.eq.0.0_WP.and.nsd.gt.0.and.sum(Nd(1:nsd)).gt.0.0_WP) then
               !  only determined species
               Neq=0.0_WP
               Neq(1:nsd)=Nd(1:nsd)
               if (cond.eq.fixed_PH) call this%hor2T(ns,Neq,this%HoR,sys%thermo,this%T)
               go to 500
            endif
            ! SBP end of added
            call die('[chem_sys initialize] All zero composition')
         endif

         ! Use initial guess N_g if provided
         if (present(N_g)) then
            call reorder_rows(N_g,sys%sp_order,N0)
            ! Guessed undetermined species
            Nu0(1:nsu)=N0(nsd+1:ns)
            ! Reduced c.v. based on N_g
            cb(1:nrc)=matmul(Nu0(1:nsu),sys%BR)
            cb_norm  =norm2(cb)
            if (cb_norm.eq.0.0_WP) go to 50
            res=norm2((cb(1:nrc)/cb_norm-cr/cr_norm))
            ! Adjust initial guess Nu0,store in Nm
            if (res.gt.res_tol) then
               N_low=sum(cr(1:sys%neu))*1e-15
               call this%min_pert(nsu,nrc,sys%BR,cr,Nu0,N_low,Nm,iret)
               ! min_pert failed
               if (iret.ne.0) then
                  call die('[chem_state initialize] min_pert failed')
                  go to 50
               endif
               Nu0=Nm
            endif

            ! Accept initial guess (Nu0); skip max-min and min_g
            this%Nd=Nd
            this%cr=cr
            go to 100

         endif

         50    continue  !  proceed without using initial guess N_g

         ! Perturb if necessary
         call this%perturb(ns,nsd,nsu,sys%ne,sys%ned,sys%neu,nrc,Nd,cr,sys%BR,sys%E,sys%diag,sys%eps_el,sys%eps_sp, &
         &                 sys%pert_tol,sys%pert_skip,this%Nd,Nm,Nupper,this%cr,npert,max_pert,iret)

         if (iret.eq.-1) then 
            write(output_unit,'(" >   chem_state perturb: non-realizable constraint = ")')
         elseif(iret.eq.-2) then 
            call die('[chem_state initialize] Perturb failed')
         endif

         ! Determine min_g composition
         call this%get_Nming(nsu,nrc,sys%BR,this%cr,gu,Ng,iret)
         if (iret.lt.0) call die('[chem_sys initialize] get_Nming failed')

         ! Form initial guess Nu0
         Nu0=Ng+frac_Nm*(Nm-Ng)

         ! Skip to here if Nu0 based on N_g accepted
         100 continue

         ! Re-estimate T0 and re-evaluate gu if required
         if ((cond.eq.fixed_PH).and.(.not.present(T_g))) then
            N1(1:nsd)   =this%Nd
            N1(nsd+1:ns)=Nu0
            call this%hor2T(ns,N1,this%HoR,sys%thermo,this%T)
            ! Set gu based on T0
            call this%get_gort(nsu,this%T,p,sys%thermo(nsd+1:ns,:),sys%P(nsd+1:ns,Gphase),gu)
         endif

         ! Set the Gibbs functin and the undetermined species moles
         this%gu=gu
         this%Nu=Nu0

         ! Determine required output
         Neq=[this%Nd,this%Nu]

         ! Jump to here if there are no undetermined species
         500	   continue

         ! Update the enthalpy
         if (cond.eq.fixed_PT) then
            call this%get_hort(ns,this%T,sys%thermo,h)
            this%HoR=sum(Neq*h)*this%T
         endif

         ! Re-order species
         do i=1,ns
            this%N(sys%sp_order(i))=Neq(i)
         end do
         if (this%cond.eq.fixed_PH) print*,'Initial guess for temperature = ',this%T
         
      end subroutine initialize


      !> Get normalized enthalpies at temperature T
      subroutine get_hort(this,ns,T,thermo,hort)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns
         real(WP), intent(in) :: T,thermo(ns,2*ncof+1)
         real(WP), intent(out) :: hort(ns)
         ! input:
         !	ns	  -number of species
         !   T     -temperature (K)
         !   thermo-thermo data for all species
         ! output:
         !   hort  -h_j/(RT) -normalized enthalpies
         ! S. B. Pope 9/26/02
         real(WP) :: th(6),Tpnm1
         integer :: k,n
         th(1)=1.0_WP  ! coefficient multipliers for enthalpy
         th(6)=1./T
         Tpnm1=1.0_WP
         do n=2,5
            Tpnm1=Tpnm1*T     ! =T.^(n-1)
            th(n)=Tpnm1/float(n)     ! =T.^(n-1) ./ n
         end do
         do k=1,ns
            if (T<thermo(k,1)) then
               hort(k)=dot_product(thermo(k,2:7),th)  ! coefficients in lower temperature range
            else
               hort(k)=dot_product(thermo(k,9:14),th) ! coefficients in upper temperature range
            endif
         end do
      end subroutine get_hort


      !> Get normalized Gibbs functions at temperature T
      subroutine get_gort(this,ns,T,p,thermo,isGas,gort)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns
         real(WP), intent(in) :: T,p,thermo(ns,2*ncof+1),isGas(ns)
         real(WP), intent(out) :: gort(ns)
         ! input:
         !   isGas -1 if gas,0 if liquid
         !   T     -temperature (K)
         !   p     -pressure (atm)
         !   thermo-thermo data for all species
         ! output:
         !   gort  -g_j/(RT) -normalized Gibbs functions
         ! S. B. Pope 9/26/02
         real(WP) :: tc(ncof),th(ncof),ts(ncof),tg(ncof)
         integer :: k,n
         if (ns.le.0) return
         tc=0.0_WP  ! coefficient multipliers for specific heats
         th=0.0_WP  ! coefficient multipliers for enthalpy
         ts=0.0_WP  ! coefficient multipliers for entropy
         tc(1)=1.0_WP
         th(1)=1.0_WP
         th(6)=1./T
         ts(1)=log(T)
         ts(7)=1.0_WP
         do n=2,5
            tc(n)=T*tc(n-1)   ! =T.^(n-1)
            th(n)=tc(n)/float(n)     ! =T.^(n-1) ./ n
            ts(n)=tc(n)/float((n-1)) ! =T.^(n-1) ./ (n-1)
         end do
         tg=th-ts
         do k=1,ns
            if (T<thermo(k,1)) then
               gort(k)=dot_product(thermo(k,2:8),tg)  ! coefficients in lower temperature range
            else
               gort(k)=dot_product(thermo(k,9:15),tg) ! coefficients in upper temperature range
            endif
         end do
         gort=gort+isGas*log(p)
      end subroutine get_gort


      !> Determine temperature given enthalpy
      subroutine hor2T(this,ns,z,hin,thermo,T)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         use messager, only: die
         implicit none
         class(chem_state),  intent(in)  :: this
         integer,            intent(in)  :: ns
         real(WP), intent(in)  :: z(ns),hin,thermo(ns,2*ncof+1)
         real(WP), intent(out) :: T
         ! input:
         !	ns		- number of species
         !   z      -moles of species
         !   hin    -enthalpy/R (K)=z'*h
         !   thermo -thermo data
         ! output:
         !   T   -temperature (K)

         ! Notes:  if the temperature is outside the range [T_low T_high]
         !   then T is returned as the closest of these bounds.
         !   If iteration fails,T is returned as T=-1.

         ! S. B. Pope 9/26/02

         integer :: itmax,it
         real(WP) :: T_tol,T0,hort(ns),hor,h_a,T_a,h_b,T_b,dT,&
            cpor(ns),hh,cpp

         itmax=100     ! maximum number of Newton iterations (usually only 3 required)
         T_tol=1e-6    ! error tolerance
         T0=1500.0_WP  ! initial guess

         !  determine if T>T0 and bracket T in [T_a T_b]
         call this%get_hort(ns,T0,thermo,hort)
         hor=dot_product(z,hort)*T0

         if (hin>hor) then	! T > T0=T_a
            h_a=hor
            T_a=T0
            call this%get_hort(ns,T_high,thermo,hort)
            h_b=dot_product(z,hort)*T_high
            if (hin.ge.h_b) then
               T=T_high   ! T > T_high (return T=T_high)
               return
            endif
            T_b=T_high
         else
            h_b=hor	! T < T0=T_b
            T_b=T0
            call this%get_hort(ns,T_low,thermo,hort)
            h_a=dot_product(z,hort)*T_low
            if (hin.le.h_a) then
               T=T_low    ! T < T_low (return T=T_low)
               return
            endif
            T_a=T_low
         endif

         !  estimate of T based on linear interpolation
         T=T_a+(hin-h_a)*(T_b-T_a)/(h_b-h_a)
            
         !  Newton iterations
         do it=1,itmax
            call this%get_cpor(ns,T,thermo,cpor)
            call this%get_hort(ns,T,thermo,hort)
            hh=dot_product(z,hort)*T
            cpp=dot_product(z,cpor)
            dT=(hin-hh)/cpp
            T=T+dT
            if (abs(dT).lt.T_tol) return  ! success
         end do

         ! Failure
         call die('[chem_state hor2T] Iterations failed')

      end subroutine hor2T


      !> Return d/dT of the normalized Gibbs functions at temperature T
      subroutine get_dgdT(this,ns,T,thermo,dgdT)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer,  intent(in)  :: ns
         real(WP), intent(in)  :: T,thermo(ns,2*ncof+1)
         real(WP), intent(out) :: dgdT(ns)
         ! input:
         !   T      - temperature (K)
         !   thermo - thermo data for all species
         ! output:
         !   dgdT   - d/dT (g_j/(RT))
         ! S. B. Pope 7/1/03
         real(WP) :: tc(ncof),th(ncof),ts(ncof),tg(ncof)
         integer  :: k,n
         tc=0.d0  ! coefficient multipliers for specific heats
         th=0.d0  ! coefficient multipliers for enthalpy
         ts=0.d0  ! coefficient multipliers for entropy
         tc(1)=1.d0/T
         th(1)=0.d0
         th(6)=-1.d0/T**2
         ts(1)=1.d0/T
         ts(7)=0.d0
         do n=2,5
            tc(n)=T*tc(n-1)                  ! =T.^(n-2)
            th(n)=tc(n)*float(n-1)/float(n)  ! =T.^(n-2) * (n-1) / n
            ts(n)=tc(n)                      ! =T.^(n-2) 
         end do
         tg=th-ts
         do k=1,ns
            if (T<thermo(k,1)) then
               dgdT(k)=dot_product(thermo(k,2:8),tg)  ! coefficients in lower temperature range
            else
               dgdT(k)=dot_product(thermo(k,9:15),tg) ! coefficients in upper temperature range
            endif
         end do
      end subroutine get_dgdT


      !> Generate (possibly) perturbed CE problem
      subroutine perturb(this,ns,nsd,nsu,ne,ned,neu,nrc,Nd,cr,BR,E,ifop,eps_el,eps_sp,pert_tol,pert_skip,zdp,zup,Nupper,crp,npert,max_pert,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         use, intrinsic :: iso_fortran_env, only: output_unit
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns,nsd,nsu,ne,ned,neu,nrc,ifop,pert_skip
         real(WP), intent(in) :: Nd(nsd),cr(nrc),BR(nsu,nrc),E(ns,ne),eps_el,eps_sp,pert_tol
         integer, intent(out) :: npert,iret
         real(WP), intent(out) :: zdp(nsd),zup(nsu),Nupper(nsu),crp(nrc),max_pert

         !  Input:
         !	ns		- number of species
         !	nsd		- number of determined species
         !	nsu		- number of undetermined species
         !	ne		- number of elements
         !	ned		- number of determined elements
         !	neu		- number of undetermined elements
         !	nrc		- number of reduced constraints
         !   Nd     -moles of determined species
         !   cr     -reduced constraint vector
         !   BR     -reduced constraint matrix
         !   E      -element matrix
         !   ifop    >0 for output
         !   eps_el -relative lower bound on element moles
         !   eps_sp -relative lower bound on species moles
         !   pert_tol- largest allowed perturbation (moles/moles of atoms)
         !  pert_skip>0 to skip perturbing undetermined species

         !  Output:
         !   zdp    -perturbed moles of determined species
         !   zup    -min-max solution for undetermined species
         !   Nupper -upper bound on undetermined species moles
         !   crp    -perturbed reduced constraint vector
         !   npert  -number of perturbations made
         !   max_pert- largest normalized perturbation made
         !   iret   = 0  successful operation
         !          =-1,non-realizable large perturbation made
         !          =-2,failed to determine max-min composition
         !          =-3,zero atoms

         integer :: j,k
         real(WP) :: zdlim,zatoms,&
         cre(neu),sumcre,cref,zed(ned),zeu(neu),ze(ne),zeu_in(neu),zau,zelow,&
         zemax,zumm(nsu),Numin,zulow,zeumax

         iret=0         ! anticipate success
         npert=0        ! number of perturbations
         max_pert=0.0_WP  ! largest normalized perturbation

         zdp=Nd       ! check that determined species are non-negative
         zatoms=0.0_WP  ! estimate of moles of atoms
         if (nsd.gt.0) then
         do j=1,ne
            zatoms=zatoms+abs(dot_product(Nd,E(1:nsd,j)))
         end do
         endif

         do j=1,neu
         zatoms=zatoms+abs(cr(j))
         end do

         if (zatoms.le.0.0_WP) then
            write(output_unit,'(" >   chem_state perturb: no atoms")')
            iret=-3
            return
         endif

         zdlim=zatoms*pert_tol
         do k=1,nsd
            if (zdp(k)<0.0_WP) then
               if (abs(zdp(k))>zdlim ) then ! significantly negative
                     max_pert=max(max_pert,abs(zdp(k))/zatoms)
                     npert=npert+1
                     if (ifop.ge.1) write(output_unit,'(" >   chem_state perturb: negative determined species")')
               endif
               zdp(k)=0.0_WP
            endif
         end do

         crp=cr
         cre=cr(1:neu) ! check that undetermined elements are positive
         sumcre=0.0_WP
         do j=1,neu
         sumcre=sumcre+abs(cre(j))
         end do
         cref=sumcre*pert_tol

         do k=1,neu
            if (cre(k)<0.0_WP) then
               if (abs(cre(k))>cref) then ! significantly negative
                     max_pert=max(max_pert,abs(cre(k))/zatoms)
                     npert=npert+1
                     if (ifop.ge.1) write(output_unit,'(" >   chem_state perturb: negative undetermined element")')
               endif
               cre(k)=0.0_WP
            endif
         end do
         crp(1:neu)=cre
                     
         zed=matmul(zdp,E(1:nsd,1:ned))                !  moles of determined elements
         zeu=matmul(zdp,E(1:nsd,ned+1:ne))+crp(1:neu)  ! moles of undetermined elements
         ze(1:ned)=zed
         ze(ned+1:ne)=zeu    ! moles of elements
         zeu_in=zeu
         zemax=maxval(ze)
         zeumax=maxval(zeu)
         zau=sum(zeu)       ! moles of atoms in undetermined species

         !  impose lower bound on moles of undetermined elements
         zelow=max(eps_el*zeumax,eps_el**2*zemax)
         do j=1,neu      
            if (zeu(j)<zelow) then
               npert=npert+1
               max_pert=max(max_pert,(zelow-zeu(j))/zatoms)
               zeu(j)=zelow
            endif
         end do

         Nupper=0.0_WP    ! determine upper bound on undetermined species
         do j=1,nsu
            Nupper(j)=1/maxval(E(nsd+j,ned+1:ne)/zeu)
         end do

         !  determine max-min moles of undetermined species

         call this%maxmin_comp(nsu,nrc,BR,crp,zumm,Numin,iret)

         if (iret<0) then
            if (ifop.ge.1) write(output_unit,'(" >   chem_state perturb: maxmin_comp failed")')
            iret=-2 
            return 
         endif

         zup=zumm
         do j=1,nsu     ! impose lower limit on undetermined species
            zulow=eps_sp*Nupper(j)
            if (zumm(j)<zulow) then
               zup(j)=zulow
               npert=npert+1
               max_pert=max(max_pert,(zup(j)-zumm(j))/zatoms)
            endif
         end do

         if (pert_skip.gt.0) return  ! do not modify constraints

         !  modify constraints according to the perturbation in undetermined species
         crp=crp+matmul(zup-zumm,BR)

         if (max_pert.gt.pert_tol) then
            if (ifop.ge.1) write(output_unit,'(" >   chem_state perturb: large perturbation made")')
            iret=-1 
         endif

      end subroutine perturb


      !> Get Ng,the value of N which minimizes g'N,subject to B'*N=c,N(i)>=0.
      subroutine get_Nming(this,nz,nc,B,c,g,Ng,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: nz,nc
         real(WP), intent(in) :: B(nz,nc),c(nc),g(nz)
         integer, intent(out) :: iret
         real(WP), intent(out) :: Ng(nz)
         !  Exit flag from linprog: iret<0 for failure.
         !  S.B. Pope 10/1/02
         integer :: i,iftest  
         real(WP) :: BT(nc,nz)
         BT=transpose(B)
         call this%solve_linprog(nz,nc,g,BT,c,Ng,iret)
         do i=1,nz  ! guard against small negative values due to round-off
         Ng(i)=max(Ng(i),0.0_WP)
         end do
      end subroutine get_Nming


      !> Determine N=N0+d which satisfies:
      !>    1) N(i).ge.eps
      !>    2) B'*N=c
      !>    3) t=max_i(|d(i)|) is minimized.
      subroutine min_pert(this,nz,nc,B,c,N0,eps,z,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         ! Input:
         !  nz -length of N
         !  nc -length of c
         !  B  -nz x nc equality constraint matrix
         !  c  -constraint vector
         !  N0 -nz-vector,N0
         !  eps-positive threshold
         ! Output:
         !  z   -solution
         !  iret=0 for success
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in)           :: nz,nc
         real(WP), intent(in)  :: B(nz,nc),c(nc),N0(nz),eps
         real(WP), intent(out) :: z(nz)
         integer, intent(out)          :: iret
         real(WP) :: A(nc+nz,3*nz),x(3*nz),f(3*nz),r(nc+nz),d(nz),tp,tm,res,Bsc,csc,zsc
         integer :: i,nx,nr
         logical :: linear=.false.
         ! x=[ u v w ]=[ z-eps (t+d)/2   (t-d)/2  ]
         ! f=[ 0 1 1 1   [ 0 0 0     1 1 1    1 1 1   ] 
         ! minimize sum(f*x)=nz*t subject to  A x=r
         nx=3*nz
         nr=nc+nz
         Bsc=max(maxval(B),-minval(B))  !  scale factors
         csc=max(maxval(c),-minval(c))
         zsc=csc/Bsc
         A(1:nz+nc,1:3*nz)=0.0_WP
         if (linear) then  !  min. sum of dz
            do i=1,nz
               A(i,i)     = 1.0_WP
               A(i,nz+i)  =-1.0_WP
               A(i,2*nz+i)= 1.0_WP
            end do
         else  !  min. sum of dz/N0
            do i=1,nz
               A(i,i)     = 1.0_WP
               A(i,nz+i)  =-(max(N0(i),eps)/zsc)**0.0_WP
               A(i,2*nz+i)= (max(N0(i),eps)/zsc)**0.0_WP
            end do
         endif
         A(nz+1:nz+nc,1:nz)  =transpose(B)/Bsc
         r(1:nz)      =(N0(1:nz) -eps) /zsc
         do i=1,nc
            r(nz+i)=(c(i)-eps*sum(B(:,i)))/csc
         end do
         f(1:nz)     =0.0_WP
         do i=1,nz
            f(nz+i)  =(eps/(eps+N0(i)))**0.0_WP
            f(2*nz+i)=f(nz+i)
         end do
         call this%solve_linprog(nx,nr,f,A,r,x,iret)
         if (iret.eq.0) then
            do i=1,nz
               z(i)=max(x(i),0.0_WP)*zsc+eps
            end do
         endif
         if (.false.) return
      end subroutine min_pert


      !> Get the normalized Cp's at temperature T
      subroutine get_cpor(this,ns,T,thermo,cpor)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns
         real(WP), intent(in) :: T,thermo(ns,2*ncof+1)
         real(WP), intent(out) :: cpor(ns)
         ! input:
         !	ns	  -number of species
         !   T     -temperature (K)
         !   thermo-thermo data for all species
         ! output:
         !   cpor  -Cp_j/R   -normalized constant-pressure specific heats
         ! S. B. Pope 9/26/02
         real(WP) :: tc(5)
         integer :: k,n
         cpor=0.0_WP
         tc(1)=1.0_WP  ! coefficient multipliers for specific heats
         do n=2,5
            tc(n)=T*tc(n-1)   ! =T.^(n-1) 
         end do
         do k=1,ns
            if (T.lt.thermo(k,1)) then
               cpor(k)=dot_product(thermo(k,2:6),tc)  ! coefficients in lower temperature range
            else
               cpor(k)=dot_product(thermo(k,9:13),tc) ! coefficients in upper temperature range
            endif
         end do
      end subroutine get_cpor


      !> Determine the max-min composition.
      subroutine maxmin_comp(this,nz,nc,B,c,Nm,zmin,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: nz,nc
         real(WP), intent(in) :: B(nz,nc),c(nc)
         integer, intent(out) :: iret
         real(WP), intent(out) :: Nm(nz),zmin
         !  Find Nm which maximizes zmin=min_i(z(i)) subject to B'*z=c.
         !  iret<0 indicates failure
         !  Method: 
         !   Initially assume zmin>=0.
         !   Define: x=[ (z-zmin)' zmin]'
         !   Maximize zmin (i.e.,minimize -x(n)) subject to x(i)>=0
         !       and B'*z=c,which is equivalent to A*x=c,
         !       where A=[ B' -sum(B)'].
         !   If feasible solution not found,zmin<0,and
         !   Define: x=[ (z-zmin)' -zmin]'
         !   Maximize zmin (i.e.,minimize x(n)) subject to x(i)>=0
         !       and B'*z=c,which is equivalent to A*x=c,
         !       where A=[ B' sum(B)'].
         !   S.B. Pope 10/1/02
         integer :: nx,tries, j,jj
         real(WP) :: A(nc,nz+1),bsum(nc),f(nz+1),x(nz+1)
         nx=nz+1
         bsum=sum(B,dim=1)
         f=0
         ! First assume zmin>0,x=[z'-zmin zmin]'
         f(nx)=-1.0_WP  ! minimize -zmin
         A(1:nc,1:nz)=transpose(B)
         A(1:nc,nx)=bsum
         call this%solve_linprog(nx,nc,f,A,c,x,iret)
         if (iret.eq.0) then  !  success,zmin>=0
            zmin=x(nx)
            Nm=x(1:nz)+zmin
            return
         elseif(iret.ne.-2) then  !  failure
            return
         endif
         ! zmin<0,re-define x=[z'-zmin -zmin]'
         f(nx)=1.0_WP  ! minimize -zmin
         A(1:nc,nx)=-bsum
         call this%solve_linprog(nx,nc,f,A,c,x,iret)
         if (iret.ne.0) return  ! failure
         zmin=-x(nx)
         Nm=x(1:nz)+zmin
      end subroutine maxmin_comp


      !> Determine x=xm which minimizes g=f'*x subject to x(i)>=0 and A*x=b,where A has full rank.
      subroutine solve_linprog(this,nx,nb,f,A,b,xm,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         use linprog, only: lp
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: nx,nb
         real(WP), intent(in)  :: f(nx),A(nb,nx),b(nb)
         real(WP), intent(out) :: xm(nx)
         integer, intent(out) :: iret
         ! Input:
         !	nx	- number of components of x
         !	nb	- number of components of b
         !	f	- nx-vector f
         !	A	- nx x nb matrix A
         ! Output:
         !	xm	- solution	
         !	iret= 0 if solution is found 
         !	iret=-1 if g is unbounded
         !	iret=-2 if there is no feasible solution
         !	iret=-3 if A is rank deficient
         !  S.B. Pope 10/1/06
         real(WP) :: eps=1e-9,ale(1,1),age(1,1),ble(1),bge(1)
         call lp(nx,0,0,nb,ale,age,A,ble,bge,b,f,xm,iret,toler=eps)
         if (iret<0) then
            !XXX write(0,*)'solve_linprog,iret=',iret  ! SBP XXX
         endif
      end subroutine solve_linprog


      !> Get the square root of moles
      function get_y(this) result(y)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%nsu) :: y
         y=exp(0.5_WP*(-this%gu+matmul(this%sys%BR,this%x(1:this%sys%nrc))+matmul(this%sys%P(this%sys%nsd+1:this%sys%ns,:),this%x(this%sys%nrc+1:this%sys%nrc+this%sys%np))))
      end function get_y


      !> Get the residual vector
      subroutine get_res(this,y)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%nsu), intent(in) :: y
         this%R(1:this%sys%nrc)=matmul(this%BtildeT,y)-this%cr
         this%R(this%sys%nrc+1:this%sys%nrc+this%sys%np)=matmul(this%PtildeT,y)-this%Nbar+this%Rd
      end subroutine get_res


      !> Initialize the solution unknowns
      subroutine x_init(this)
         use mathtools, only: lss
         use messager,  only: die
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP), dimension(:), allocatable :: rhs,lam
         integer :: info
         ! Allocate intermediate arrays
         allocate(rhs(this%sys%nrc))
         allocate(lam(this%sys%nrc))
         ! Get the species and phase moles
         this%Nbar=matmul(transpose(this%sys%P),[this%Nd,this%Nu])
         ! Calculate the contribution of Nd in the residual
         this%Rd=matmul(transpose(this%sys%P(1:this%sys%nsd,:)),this%Nd)
         ! Calculate the Lagrange multipliers
         call this%get_gort(this%sys%nsu,this%T,this%p,this%sys%thermo(this%sys%nsd+1:this%sys%ns,:),this%sys%P(this%sys%nsd+1:this%sys%ns,Gphase),this%gu)
         rhs=log(this%Nu)-matmul(this%sys%P(this%sys%nsd+1:this%sys%ns,:),log(this%Nbar))+this%gu
         call lss(this%sys%nsu,this%sys%nrc,this%sys%BR,rhs,lam,info)
         if (info.ne.0) call die('[chem_state x_init] Least squares solver for lambda initialization failed.')
         ! Set the initial solution vector
         this%x(1:this%sys%nrc)=lam
         this%x(this%sys%nrc+1:this%sys%nrc+this%sys%np)=log(this%Nbar)
         ! Deallocate intermediate arrays
         deallocate(rhs,lam)
      end subroutine x_init


      !> Find the chemical equilibium state
      subroutine equilibrate(this)
         use messager, only: die
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP), dimension(:), allocatable :: Neq
         integer :: isc
         ! Allocate arrays
         allocate(Neq(this%sys%ns))
         ! Get the chemical equilibrium state
         call this%get_ceq(Neq)
         ! Reorder the composition
         do isc=1,this%sys%ns
            this%N(this%sys%sp_order(isc))=Neq(isc)
         end do
         ! Deallocate arrays
         deallocate(Neq)
      end subroutine equilibrate


      !> Find the chemical equilibium state at constant pressure and temperature
      subroutine get_ceq_PT(this,Neq)
         use messager, only: die
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%ns), intent(inout) :: Neq
         integer :: i,j,iJ,jJ,info
         real(WP), dimension(:,:), allocatable :: Jac
         real(WP), dimension(:),   allocatable :: dx,y
         real(WP), dimension(:,:), allocatable :: BTB,PTP,BTP
         real(WP), dimension(:),   allocatable :: S,work
         real(WP) :: Rnorm,rcond
         integer  :: rank,lwork
         ! Allocate arrays
         allocate(Jac   (this%sys%nrc+this%sys%np,this%sys%nrc+this%sys%np)); Jac=0.0_WP
         allocate(dx    (this%sys%nrc+this%sys%np));  dx=0.0_WP
         allocate(BTB   (this%sys%nrc,this%sys%nrc)); BTB=0.0_WP
         allocate(BTP   (this%sys%nrc,this%sys%np));  BTP=0.0_WP
         allocate(PTP   (this%sys%np,this%sys%np));   PTP=0.0_WP
         allocate(y     (this%sys%nsu));              y=0.0_WP
         allocate(S     (this%sys%nrc+this%sys%np))
         lwork=10*(this%sys%nrc+this%sys%np)
         allocate(work(lwork))
         rcond=-1.0_WP
         call this%x_init()
         ! Newton-Raphson
         this%iter_N=0
         y=this%get_y()
         Rnorm=10.0_WP*this%tol_N
         do while(Rnorm.ge.this%tol_N)
            ! Increment iteration number
            this%iter_N=this%iter_N+1
            if (this%iter_N.gt.this%iter_N_max) then
               this%iter_N=this%iter_N-1
               call die('[chem_state get_ceq_PT] Newton solver reached maximum number of iterations')
            end if
            ! Build the Jacobian matrix
            call this%get_BP(y)
            BTB=matmul(this%BtildeT,this%Btilde)
            PTP=matmul(this%PtildeT,this%Ptilde)
            BTP=matmul(this%BtildeT,this%Ptilde)
            do j=1,this%sys%nrc
               jJ=j
               do i=1,j
                  iJ=i
                  Jac(iJ,jJ)=BTB(i,j)
               end do
            end do
            do j=1,this%sys%np
               jJ=this%sys%nrc+j
               do i=1,this%sys%nrc
                  iJ=i
                  Jac(iJ,jJ)=BTP(i,j)
               end do
            end do
            do j=1,this%sys%np
               jJ=this%sys%nrc+j
               do i=1,j
                  iJ=this%sys%nrc+i
                  Jac(iJ,jJ)=PTP(i,j)
               end do
            end do
            do j=1,this%sys%nrc+this%sys%np-1
               do i=j+1,this%sys%nrc+this%sys%np
                  Jac(i,j)=Jac(j,i)
               end do
            end do
            do i=1,this%sys%np
               iJ=this%sys%nrc+i
               jJ=this%sys%nrc+i
               Jac(iJ,jJ)=Jac(iJ,jJ)-this%Nbar(i)
            end do
            ! Get the residual error
            call this%get_res(y)
            Rnorm=norm2(this%R)
            ! Solve for dx
            dx=-this%R
            call dgelss(this%sys%nrc+this%sys%np,this%sys%nrc+this%sys%np,1,Jac,this%sys%nrc+this%sys%np,dx,this%sys%nrc+this%sys%np,S,rcond,rank,work,lwork,info)
            if (rank.ne.this%sys%nrc+this%sys%np) call die('[chem_state get_ceq_PT]: Jacobian is not full rank')
            if (info.ne.0) call die('[chem_state get_ceq_PT]: Least-squares solver failed')
            ! Update the solution
            this%x=this%x+dx
            ! Get the species and phase moles
            y=this%get_y()
            this%Nu=y*y
            this%Nbar=exp(this%x(this%sys%nrc+1:this%sys%nrc+this%sys%np))
         end do
         ! Assemble the composition
         Neq=[this%Nd,this%Nu]
         ! Deallocate arrays
         deallocate(Jac,dx,BTB,BTP,PTP,y,S,work)
      end subroutine get_ceq_PT


      !> Find the chemical equilibium state at constant pressure and enthalpy
      subroutine get_ceq_PH(this,Neq)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         use messager, only: die
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%ns), intent(inout) :: Neq
         real(WP), dimension(:), allocatable :: hort
         real(WP) :: Tn,Tlo,Thi
         real(WP) :: HoR0,hlo,hhi,Cp_eff
         real(WP), dimension(:), allocatable :: xold
         ! Allocate arrays
         allocate(hort(this%sys%ns))
         allocate(xold(this%sys%nrc+this%sys%np))
         ! Initialize
         HoR0=this%HoR
         this%dT=1e5*this%tol_T
         this%iter_T=0
         Tlo=-1e30 ! Lowest temperature at which h has been evaluated
         Thi= 1e30 ! Highest temperature at which h has been evaluated
         ! Iterate over temperature
         do while(abs(this%dT/this%T).ge.this%tol_T)
            ! Increment the iterations
            this%iter_T=this%iter_T+1
            if (this%iter_T.gt.this%iter_T_max) then
               this%iter_T=this%iter_T-1
               call die('[chem_state get_ceq_PH] Temperature solver reached maximum number of iterations')
               exit
            end if
            ! Update the phase moles
            this%Nbar=exp(this%x(this%sys%nrc+1:this%sys%nrc+this%sys%np))
            ! Determine equilibrium composition at current temperature
            call this%get_gort(this%sys%nsu,this%T,this%p,this%sys%thermo(this%sys%nsd+1:this%sys%ns,:),this%sys%P(this%sys%nsd+1:this%sys%ns,Gphase),this%gu)
            xold=this%x
            call this%get_ceq_PT(Neq)
            ! Get the effective Cp
            call this%get_Cp_eff(Cp_eff)
            ! Obtain species h/(RT)
            call this%get_hort(this%sys%ns,this%T,this%sys%thermo,hort)
            ! Mixture H/R
            this%HoR=this%T*sum(Neq*hort)
            ! Predict dT
            this%dT=(HoR0-this%HoR)/Cp_eff
            ! Check that T is within limits
            if (this%T.eq.T_high.and.this%dT.gt.0.0_WP) call die('[chem_state get_ceq_PH] T > T_high')
            if (this%T.eq.T_low .and.this%dT.lt.0.0_WP) call die('[chem_state get_ceq_PH] T < T_low')
            ! Ensure that Tn is within limits
            Tn=this%T+this%dT
            Tn=max(min(Tn,T_high),T_low)
            ! Use linear interpolation instead if Tn is closer to known bound
            if (this%dT.gt.0.0_WP) then
               Tlo=this%T
               hlo=this%HoR
               if (Tn.gt.0.5_WP*(Tlo+Thi)) then
                  Tn=Tlo+(Thi-Tlo)*(HoR0-hlo)/(hhi-hlo)
               endif
            else
               Thi=this%T
               hhi=this%HoR
               if (Tn.lt.0.5_WP*(Tlo+Thi)) then
                  Tn=Tlo+(Thi-Tlo)*(HoR0-hlo)/(hhi-hlo)
               endif
            endif
            ! Update temperature increment
            this%dT=Tn-this%T   
            ! Update temperature
            this%T=Tn
         end do
         ! Assemble the composition
         Neq=[this%Nd,this%Nu]
         ! Dellocate arrays
         deallocate(hort)
      end subroutine get_ceq_PH


      !> Evaluate the effective specific heat
      !> Cp_eff = (dh/dT)/R at constant constraints and p
      subroutine get_Cp_eff(this,Cp_eff)
         class(chem_state), intent(inout) :: this
         real(WP), intent(out) :: Cp_eff
         real(WP), dimension(:), allocatable :: dgudT,dxdT,cpor,hort,N,dNddT,dNudT,dNdT
         ! Allocate arrays
         allocate(dgudT(this%sys%nsu))
         allocate(dxdT(this%sys%nrc+this%sys%np))
         allocate(cpor(this%sys%ns))
         allocate(hort(this%sys%ns))
         allocate(N(this%sys%ns))
         allocate(dNddT(this%sys%nsd))
         allocate(dNudT(this%sys%nsu))
         allocate(dNdT(this%sys%ns))
         ! Get d(gu/(RT))/dT
         call this%get_dgdT(this%sys%nsu,this%T,this%sys%thermo(this%sys%nsd+1:this%sys%ns,:),dgudT)
         ! Get d(lambda)/dT and d(ln(Nbar))/dT
         call this%get_dxdT(dgudT,dxdT)
         ! Get the Gibbs and enthalpy
         call this%get_cpor(this%sys%ns,this%T,this%sys%thermo,cpor)
         call this%get_hort(this%sys%ns,this%T,this%sys%thermo,hort)
         ! Rates of change of moles
         N=[this%Nd,this%Nu]
         dNddT=0.0_WP
         dNudT=this%Nu*(-dgudT+matmul(this%sys%BR,dxdT(1:this%sys%nrc))+dxdT(this%sys%nrc+1:this%sys%nrc+this%sys%np))
         dNdT=[dNddT,dNudT]
         Cp_eff=sum(cpor*N)+this%T*sum(hort*dNdT)
         Cp_eff=sum(cpor*N)
         ! Deallocate arrays
         deallocate(dgudT,dxdT,cpor,hort,N,dNddT,dNudT,dNdT)
      end subroutine get_Cp_eff


      !> Get the temperature derivative of the solution vector
      subroutine get_dxdT(this,dgudT,dxdT)
         use messager,  only: die
         use mathtools, only: lss
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%nsu), intent(in) :: dgudT
         real(WP), dimension(this%sys%nrc+this%sys%np), intent(out) :: dxdT
         real(WP), dimension(:),   allocatable :: lamdotg,Sig,Sinv,work,y,Ygdot,dlnNbardT,rhs
         real(WP), dimension(:,:), allocatable :: Btildeinv,lamdoty,U,VT,M,Btilde_cp
         real(WP) :: srlim=1e-9
         integer :: info,lwork,i,n_small
         ! Allocate arrays
         allocate(Btildeinv(this%sys%nrc,this%sys%nsu))
         allocate(lamdotg  (this%sys%nrc))
         allocate(Sig      (this%sys%nrc))
         allocate(Sinv     (this%sys%nrc))
         allocate(work     (20*(this%sys%nsu+this%sys%nrc)))
         allocate(y        (this%sys%nsu))
         allocate(Ygdot    (this%sys%nsu))
         allocate(dlnNbardT(this%sys%np))
         allocate(rhs      (this%sys%np))
         allocate(lamdoty  (this%sys%nrc,this%sys%np))
         allocate(U        (this%sys%nsu,this%sys%nrc))
         allocate(VT       (this%sys%nrc,this%sys%nrc))
         allocate(M        (this%sys%np,this%sys%np))
         allocate(Btilde_cp(this%sys%nsu,this%sys%nrc))
         lwork=size(work)
         ! Get the y vector and update the coefficient matrices
         y=this%get_y()
         call this%get_BP(y)
         Btilde_cp=this%Btilde(1:this%sys%nsu,1:this%sys%nrc)
         ! Get the SVD of Btilde
         call dgesvd('S','A',this%sys%nsu,this%sys%nrc,Btilde_cp,this%sys%nsu,Sig(1:this%sys%nrc), &
         &           U(1:this%sys%nsu,1:this%sys%nrc),this%sys%nsu,VT(1:this%sys%nrc,1:this%sys%nrc),this%sys%nrc,work(1:lwork),    &
         &           lwork,info)
         ! call dgesvd('S','A',this%sys%nsu,this%sys%nrc,Btilde_cp,this%sys%nsu,Sig,U,this%sys%nsu,VT,this%sys%nrc,work,lwork,info)
         if (info.ne.0) call die('[chem_state get_dxdT] SVD of B tilde failed')
         ! Get the inverse of Sigma
         call get_Sinv(this%sys%nrc,Sig,Sinv,srlim,n_small)
         ! Store Sinv * V' in VT
         do i=1,this%sys%nrc
            VT(i,:)=Sinv(i)*VT(i,:)
         end do
         ! Btilde^-1 = V * Sinv * U'
         Btildeinv=transpose(matmul(U,VT))
         ! Solve for dxdT
         Ygdot=y*dgudT
         lamdotg=matmul(Btildeinv,Ygdot)
         lamdoty=matmul(Btildeinv,this%Ptilde)
         M=matmul(this%PtildeT,matmul(this%Btilde,lamdoty))
         rhs=matmul(this%PtildeT,matmul(this%Btilde,lamdotg)-Ygdot)
         call lss(this%sys%np,this%sys%np,M,rhs,dlnNbardT,info)
         if (info.ne.0) call die('[chem_state get_dxdT] Least squares solver failed')
         dxdT(1:this%sys%nrc)=lamdotg-matmul(lamdoty,dlnNbardT)
         dxdT(this%sys%nrc+1:this%sys%nrc+this%sys%np)=dlnNbardT
         ! Deallocate arrays
         deallocate(Btildeinv,lamdotg,lamdoty,Sig,Sinv,work,y,Ygdot,dlnNbardT,rhs,U,VT,M,Btilde_cp)
         contains
            ! Get the inverse of Sigma
            subroutine get_Sinv(n,S,Si,srat_lim,n_s)
               ! Given the n-vector of singular values, S, which are in decreasing order,
               ! return the n-vector of pseudo-inverses, Si, and the number n_s of
               ! small singular values.  The j-th singular value is deemed to be small
               ! if S(j)/S(1) < srat_lim.
               integer,  intent(in)  :: n
               real(WP), intent(in)  :: S(n),srat_lim
               integer,  intent(out) :: n_s
               real(WP), intent(out) :: Si(n)
               integer  :: j
               real(WP) :: slim
               slim=srat_lim*S(1)
               do j=1,n
                  if (S(j).gt.slim) then
                     Si(j)=1.0_WP/S(j)
                     n_s=n-j
                  else
                     Si(j)=0.0_WP
                  endif
               end do
            end subroutine get_Sinv
      end subroutine get_dxdT


      !> Get Btilde and Ptilde
      subroutine get_BP(this,y)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%nsu), intent(in) :: y
         integer :: j
         do j=1,this%sys%nrc
            this%Btilde(:,j)=y*this%sys%BR(:,j)
         end do
         do j=1,this%sys%np
            this%Ptilde(:,j)=y*this%sys%P(this%sys%nsd+1:this%sys%ns,j)
         end do
         this%BtildeT=transpose(this%Btilde)
         this%PtildeT=transpose(this%Ptilde)
      end subroutine get_BP


end module chem_state_class