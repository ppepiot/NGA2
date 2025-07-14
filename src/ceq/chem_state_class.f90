!> Chemical state
module chem_state_class
   use precision,      only: WP
   use chem_sys_class, only: chem_sys,Lphase,Gphase,ncof
   implicit none
   private
   
   !> Expose type/constructor/methods
   public :: chem_state

   !> Chemical statete object definition
   type :: chem_state

      ! This is our chemical system
      class(chem_sys), pointer :: sys            !< This is the chemical system the solver is build for

      ! Thermochemical quantities
      real(WP) :: p                              !< Pressure
      real(WP) :: T                              !< Temperature
      real(WP) :: h                              !< Enthalpy
      real(WP), dimension(:),allocatable :: Nd   !< Moles of determined species (nsd)
      real(WP), dimension(:),allocatable :: Nu   !< Moles of undetermined species (nsu)
      real(WP), dimension(:),allocatable :: lam  !< Lagrange multipliers (nrc)
      real(WP), dimension(:),allocatable :: cr   !< Reduced constraint vector
      real(WP), dimension(:),allocatable :: Q    !< Consistency vector Q'*Xu=1 (nsu)
      real(WP), dimension(:),allocatable :: gu   !< Gibbs functions of undetermined species (nsu)

      ! Numerical parameters
      integer :: temp_its                        !< Number of temperature iterations performed
      integer :: time_steps                      !< Number of pseudo time steps
      integer :: newt_calls                      !< Number of Newton solves attempted
      integer :: newt_its                        !< Number of Newton iterations performed
      integer :: nyeval                          !< Number of y evaluations

   contains
      procedure :: initialize
      procedure :: get_hort
      procedure :: get_gort
      procedure :: hor2T
      procedure :: perturb
      procedure :: get_Nming
      procedure :: min_pert
      procedure :: get_cpor
      procedure :: maxmin_comp
      procedure :: solve_linprog
      procedure :: eqiulibrate
   end type chem_state


   contains


      !> Chemical state initializer
      subroutine initialize(this,sys,c,N,p_atm,p_Pa,p_cgs,T,HoR,N_h,T_h,N_g,T_g,N_eq,T_eq,HoR_eq,stats,info)

         !  Initializes the constrained equilibrium state of 
         !  an ideal liquid-gas mixture consisting of ns species,either at fixed pressure and
         !  temperature (p,T),or at fixed pressure and enthalpy (p,H).

         !  The nc equality constraints are written:  B'*N=c,where B is the
         !  ns x nc basic constraint matrix,N is the ns-vector of species moles,
         !  and c is the nc-vector of constraint values.

         !  On return,info should be checked. (info.ge.0 for success)

         !  Input:
         !    sys -object type chem_sys created by subroutine.
         !           Required first argument

         !    (Either c or N must be specified.)
         !    c   -values of the nc basic constraints (real(nc))
         !    N   -moles of species used to calculate c as : c=B'*N (real(ns))

         !    (Either p_atm or p_Pa or p_cgs must be specified.)
         !    p_atm-pressure in standard atmospheres (real)
         !    p_Pa -pressure in Pascals (SI units) (real)
         !    p_cgs-pressure in cgs units (dynes/cm^2) (real)

         !    (For a fixed (p,T) equilibrium calculation,specify T: do not specify HoR,N_h or T_h.)
         !    T-temperature (K) for fixed-temperature problem

         !    (For a fixed (p,H) equilibrium calculation,specify either HoR or N_h and T_h: 
         !     do not specify T.)
         !    HoR the fixed value of H/R [moles K],where H=enthalpy,R=universal gas constant.
         !    N_h species moles used to calculate H  (real(ns))
         !    T_h temperature used to calculate H as:  H/R=sum(N_h h(T_h)/R),where
         !        h(T_h)/R [which has dimensions K] is the molar specific species enthalpy.

         !    (Initial guesses are not needed,and should not be specified unless they
         !     are good guesses.)
         !    N_g  -initial guess for species moles
         !    T_g  -initial guess for temperature (for fixed (p,H) only)

         !  Output:
         !    N_eq    -moles of species in equilibrium mixture (real(ns))
         !    T_eq    -temperature of the equilibrium mixture (real)
         !    HoR_eq  -value of H/R [moles K] for the equilibrium mixture (real)
         !    stats   -numerical information about the solution-details below (real(20))
         !    info    -information about the solution:
         !             > 0,info=the number of Newton iterations performed
         !             < 0,an error occurred-see details below

         !  Error conditions:  
         !    info=-1  sys is not associated (ceq_sys_init must be called first)
         !    info=-2  sys has not been initialized (ceq_sys_init must be called first)
         !    info=-3  the pressure has been specified more than once
         !    info=-4  the pressure has not been specified
         !    info=-5  the pressure is non-positive
         !    info=-6  the specified temperature is outside the range (T_low,T_high)
         !    info=-7  neither T nor enthalpy specified
         !    info=-8  the guess for temperature T_g is outside the range (T_low,T_high)
         !    info=-9  neither c nor N has not been specified
         !    info=-10 the constraints (and determined species) are zero
         !    info=-11 the constraints are not realizable
         !    info=-12 (maxmin_comp failed)
         !    info=-13 (get_Nming failed)
         !    info=-14 (get_hort2T failed)
         !    info=-15 (ceq_lamg_init failed)
         !    info=-16 (ceq_fixed_T failed)
         !    info=-17 (ceq_fixed_h failed)
         !    info=-18 T_eq is greater than T_high
         !    info=-19 T_eq is less than T_low

         !  Values of info between -12 and -17 are internal failures (which could be caused by erroneous
         !  input).  Other values are caused by erroneous input.  A temperature outside the range
         !  (T_low,T_high) may be caused by erroneous input (e.g.,of T or HoR).  If the correct
         !  temperature is outside this range,the values of T_low and T_high can be changed prior to
         !  calling ceq_state by calling ceq_param_set,e.g.,call ceq_param_set(sys,T_low=200.0_WP).
         !  To diagnose an error condition,the case can be repeated with diagnostics turned on,
         !  by: call ceq_param_set(sys,diag=5).

         ! stats:
         !    stats(1)=temp_its  -number of temperature iterations
         !    stats(2)=time_steps-number of psuedo-time steps
         !    stats(3)=newt_calls-number of Newton calls (iteration not performed if initial residual small)
         !    stats(4)=newt_its  -number of Newton iterations
         !    stats(5)=nyeval    -number of temperature iterations
         !    stats(6)=err       -last residual error
         !    stats(7)=npert     -number of quantities perturbed
         !    stats(8)=max_pert  -maximum value of perturbation
         !    stats(9)=zmm       -minimum maxmin composition

         use messager,  only: die
         use mathtools, only: reorder_rows
         implicit none
         class(chem_state), intent(inout) :: this
         class(chem_sys), target, intent(in) :: sys
         real(WP), intent(in),   optional :: c(sys%nc),N(sys%ns),p_atm,p_Pa,p_cgs,T,HoR,N_h(sys%ns),T_h,N_g(sys%ns),T_g
         real(WP), intent(out),  optional :: N_eq(sys%ns),T_eq,HoR_eq,stats(20)
         integer,   intent(out), optional :: info
         integer  :: nb,nc,ns,nsd,nsu,nrc,npert,iret,i,lu
         real(WP) :: p,T0,errc,max_pert,Ndbar,Nu_atoms,zumin,zumax,cb(sys%nc),cmod(sys%nb),Nd(sys%nsd),err,cr_norm,cb_norm,res,z_low,res_tol,err2
         real(WP), dimension(sys%ns)  :: N0,N1,h,Neq
         real(WP), dimension(sys%nsu) :: Nu,zu0,zm,zupper,atoms,zg,gu,xu,rhs,gu0,zu00
         real(WP), dimension(sys%nrc) :: cr,lam0
         real(WP), parameter :: p0_Pa=1.01325d5 ! standard atmosphere [Pa]
         logical :: fixed_T,fail,diag

         ! Point to chem_sys object
         this%sys=>sys

         if(.not.present(info)) then
            write(0,*)'chem_state initializer: argument INFO must be present; stopping'
            call die('chem_state initializer: argument INFO must be present')
         endif

         ! Anticipate success
         info =0
         if(present(stats)) stats=0.0_WP
         err=0.0_WP

         ! Check that sys has been initialized  =====================================

         if(.not.sys%initialized) then
            info=-2
            call die('')
         endif

         ! Obtain indexes and allocate state
         nb =sys%nb
         nc =sys%nc
         nrc=sys%nrc
         ns =sys%ns
         nsd=sys%nsd
         nsu=sys%nsu

         allocate(this%Nd(sys%nsd))
         allocate(this%Nu(sys%nsu))
         allocate(this%lam(sys%nrc))
         allocate(this%cr(sys%nrc))
         allocate(this%Q(sys%nsu))
         allocate(this%gu(sys%nsu))
         
         this%temp_its  =0 ! number of temperature iterations performed
         this%time_steps=0 ! number of pseudo time steps
         this%newt_calls=0 ! number of Newton solves attempted
         this%newt_its  =0 ! number of Newton iterations performed
         this%nyeval    =0 ! number of y evaluations

         ! Determine if errors are to be reported
         if(sys%diag.ge.1) then
            diag=.true.
         else
            diag=.false.
         endif
         lu=sys%lu

         ! Determine pressure  p  in standard atmospheres  !=================
         if(present(p_atm)) then
            p=p_atm
            if(present(p_Pa).or.present(p_cgs)) then
               if(diag) then
                  write(lu,*) 'ceq_state: multiple input pressures'
               end if
               info=-3
               call die('chem_sys initialize: Multiple input pressures')
            endif

         elseif(present(p_Pa)) then
            p=p_Pa/p0_Pa 
            if(present(p_atm).or.present(p_cgs)) then
               if(diag) then
                  write(lu,*) 'ceq_state: multiple input pressures'
               end if
               info=-3
               call die('chem_sys initialize: Multiple input pressures')
            endif

         elseif(present(p_cgs)) then
            p=p_cgs*0.1_WP/p0_Pa  
            if(present(p_atm).or.present(p_Pa)) then
               if(diag) then
                  write(lu,*) 'ceq_state: multiple input pressures'
               end if
               info=-3
               call die('chem_sys initialize: Multiple input pressures')
            endif
         else
            if(diag) then
               write(lu,*) 'ceq_state: no pressure specified'
            end if
            info=-4
            call die('chem_sys initialize: No pressure specified')
         endif

         if(p.le.0.0_WP) then
            if(diag) then
               write(lu,'(a,1p,9e13.4)') 'ceq_state: pressure must be strictly positive: p=',p
            end if
            info=-5
            call die('')
         endif

         this%p=p

         ! Determine T or HoR  =============================================

         fixed_T=.false.
         if(present(T)) then
            this%T=T
            fixed_T=.true.

            if(T.lt.sys%T_low.or.T.gt.sys%T_high) then
               if(diag) write(lu,'(a,1p,9e13.4)') &
                     'ceq_state: temperature out of range: T,T_low,T_high=',T,sys%T_low,sys%T_high
               info=-6
               call die('chem_sys initialize: ')
            endif

         elseif(present(HoR)) then
            this%h=HoR

         elseif(present(N_h).and.present(T_h)) then
            call reorder_rows(N_h,sys%sp_order,N0)
            call this%get_hort(sys%ns,T_h,sys%thermo,h)
            this%h=sum(N0*h)*T_h  ! HoR

         else
            if(diag) write(lu,*) 'ceq_state: neither T nor HoR nor (N_h and T_h) specified'
            info=-7
            call die('chem_sys initialize: ')
         endif

         ! Determine T0 for initial guess,and evaluate gu0=gu(T0) ===============

         if(fixed_T) then
            T0=this%T  ! specified fixed T

         elseif(present(T_g)) then
            T0=T_g  ! guess provided

            if(T_g.lt.sys%T_low.or.T_g.gt.sys%T_high) then
               if(diag) write(lu,'(a,1p,9e13.4)') 'ceq_state: T_g out of range: T_g,T_low,T_high=',&
                                          T_g,sys%T_low,sys%T_high
               info=-8
               call die('chem_sys initialize: ')
            endif

         else
            T0=sqrt(sys%T_low*sys%T_high)
            T0=max(T0,0.1_WP*sys%T_high)
         endif

         call this%get_gort(nsu,T0,p,sys%thermo(nsd+1:ns,:),sys%P(nsd+1:ns,Gphase),gu)

         ! Form the basic constraint vector  =========================

         if(present(c)) then
            cb(1:nc)=c

         elseif(present(N)) then
            cb(1:nc)=matmul(N,sys%B)

         else
            if(diag) write(lu,*) 'ceq_state: neither c nor N specified'
            info=-9
            call die('chem_sys initialize: ')
         endif

         ! Form modified and reduced constraints  !===================

         cmod(1:nb)=matmul(sys%A(1:nb,1:nc) ,cb)   ! form modified constraint vector
         Nd(1:nsd) =cmod(1:nsd)                    ! determined species

         ! Treat the special case of no undetermined species

         if(nsu.eq.0) then
            Neq=Nd
            if(.not.fixed_T) then
               call this%hor2T(ns,Neq,this%h,sys%T_low,sys%T_high,sys%thermo,this%T,iret)
               if(iret.lt.0) then
                  info=-17+iret
                  call die('chem_sys initialize: ')
               endif
            endif
            
            go to 500
         endif
         
         cr(1:nrc) =cmod(nsd+1:nb)   ! reduced constraint vector
         cr_norm   =norm2(cr)        ! |cr|

         if(cr_norm.le.0.0_WP) then
         ! SBP added 4/9/2009
            if(cr_norm.eq.0.0_WP.and.nsd.gt.0.and.sum(Nd(1:nsd)).gt.0.0_WP) then
               !  only determined species
               Neq=0.0_WP
               Neq(1:nsd)=Nd(1:nsd)
               if(.not.fixed_T) then
                  call this%hor2T(ns,Neq,this%h,sys%T_low,sys%T_high,sys%thermo,this%T,iret)
                  if(iret.lt.0) then
                     info=-17+iret
                     call die('chem_sys initialize: ')
                  endif
               endif
         
               go to 500
            endif
            ! SBP end of added
         
            if(diag) write(lu,*) 'ceq_state: constraints are zero'
            info=-10  ! SBP bug fix 4/8/2009
            call die('chem_sys initialize: ')
         endif

         ! Use initial guess N_g if provided  !====================

         if(present(N_g)) then
            call reorder_rows(N_g,sys%sp_order,N0)
            zu0(1:nsu)=N0(nsd+1:ns)  ! guessed undetermined species

            cb(1:nrc)=matmul(zu0(1:nsu),sys%BR) ! reduced c.v. based on N_g
            cb_norm  =norm2(cb)

            if(cb_norm.eq.0.0_WP) then
               if(sys%diag .ge.2) write(lu,*)'ceq_state: N_g yields zero constraint'
               go to 50
            endif  

            res=norm2((cb(1:nrc)/cb_norm-cr/cr_norm))  !  residual 

            if(sys%diag.ge.2) write(lu,'(a,1p,2e13.4)')'ceq_state: initial guess residual=',res

            if(res.gt.sys%res_tol) then  !  adjust initial guess zu0,store in zm

               z_low=sum(cr(1:sys%neu))*1e-15

               call this%min_pert(nsu,nrc,sys%BR,cr,zu0,z_low,zm,info)
               if(info.ne.0) then  !  min_pert failed
                     if(sys%diag.ge.2) write(lu,*)'ceq_state: min_pert failed,info=',info
                  go to 50
               endif

               zu0=zm

               if(sys%diag.ge.2) then
                  cb(1:nrc)=matmul(zu0(1:nsu),sys%BR) ! reduced c.v. based on N_g
                  cb_norm  =norm2(cb)
                  res      =norm2((cb(1:nrc)/cb_norm-cr/cr_norm))
                  write(0,'(a,1p,9e13.4)')  &
                     'ceq_state: after min_pert: res,min(zu0)=',res,minval(zu0)
               endif
            endif

            ! Accept initial guess (zu0); skip max-min and min_g

            this%Nd=Nd
            this%cr=cr

            go to 100
         endif

         50    continue  !  proceed without using initial guess N_g

         ! Perturb if necessary  =====================================

         call this%perturb(ns,nsd,nsu,sys%ne,sys%ned,sys%neu,nrc,&
                           Nd,cr,sys%BR,sys%E,sys%diag,sys%lu, &
                     sys%eps_el,sys%eps_sp,sys%pert_tol,sys%pert_skip,&
                           this%Nd,zm,zupper,this%cr,npert,max_pert,iret)

         if(iret.eq.-1) then 
            if(diag) write(lu,*) 'ceq_state: non-realizable constraint'
            info=-11
            call die('chem_sys initialize: ')

         elseif(iret.eq.-2) then 
            if(diag) write(lu,*) 'ceq_state: ceq_max-min failed'
            info=-12
            call die('chem_sys initialize: ')
         endif

         if(present(stats)) then
            stats(7)=npert
            stats(8)=max_pert
            stats(9)=minval(zm)
         endif


         ! Determine min_g composition based on T0  !===========

         call this%get_Nming(nsu,nrc,sys%BR,this%cr,gu,zg,iret)
         
         if(iret.lt.0) then
            if(diag) write(lu,*) 'ceq_state: get_Nming failed'
            info=-13
            call die('chem_sys initialize: ')
         endif

         ! Form initial guess zu0   !============================

         zu0=zg+sys%frac_zm*(zm-zg)

         100 continue  !  skip to here if zu0 based on N_g accepted

         ! Determine normalization-condition vector Q !=========
         !     q=Q'*Xu -1=0
         Ndbar      =sum(this%Nd)               ! total moles of determined species
         atoms   =sum(sys%E(nsd+1:ns,:),2)   ! atoms in undetermined species
         Nu_atoms=sum(this%cr(1:sys%neu))    ! total moles of atoms in undetermined species
         this%Q  =1.0_WP+atoms*Ndbar/Nu_atoms   ! consistency condition vector

         ! Re-estimate T0 and re-evaluate gu if required !======

         if(.not.fixed_T.and..not.present(T_g)) then
            N1(1:nsd)   =this%Nd
            N1(nsd+1:ns)=zu0
            call this%hor2T(ns,N1,this%h,sys%T_low,sys%T_high,sys%thermo,T0,iret)

            if(iret.eq.-3) then
                  if(diag) write(lu,*) 'ceq_state: hor2T failed'
               info=-14
               call die('chem_sys initialize: ')
            endif

            call this%get_gort(nsu,T0,p,sys%thermo(nsd+1:ns,:),sys%P(nsd+1:ns,Gphase),gu)  ! set gu based on T0
         endif

         this%gu=gu
         this%Nu=zu0
         

         ! Determine required output  =======================================

         Neq=(/ this%Nd ,this%Nu /)  !  equilibrium moles (re-ordered)

         500	   continue  !  jump to here if there are no undetermined species

         if(present(T_eq)) T_eq=this%T

         if(present(HoR_eq)) then
            if(fixed_T) then
               call this%get_hort(ns,this%T,sys%thermo,h)
               this%h=sum(Neq*h)*this%T
         endif

         HoR_eq=this%h
         endif

         if(present(N_eq)) then
            do i=1,ns          ! re-order species
               N_eq(sys%sp_order(i))=Neq(i)
            end do
         endif

         if(nsu.gt.0) then
            if(present(stats)) then
               stats(1)=this%temp_its
               stats(2)=this%time_steps
               stats(3)=this%newt_calls
               stats(4)=this%newt_its
               stats(5)=this%nyeval
            stats(6)=err
            endif

            info=this%nyeval
         endif
         
      end subroutine initialize


      !> Get normalized enthalpies at temperature T
      subroutine get_hort(this,ns,T,thermo,hort)
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns
         real(kind(1.0_WP)), intent(in) :: T,thermo(ns,2*ncof+1)
         real(kind(1.0_WP)), intent(out) :: hort(ns)
         ! input:
         !	ns	  -number of species
         !   T     -temperature (K)
         !   thermo-thermo data for all species
         ! output:
         !   hort  -h_j/(RT) -normalized enthalpies
         ! S. B. Pope 9/26/02
         real(kind(1.0_WP)) :: th(6),Tpnm1
         integer :: k,n
         th(1)=1.0_WP  ! coefficient multipliers for enthalpy
         th(6)=1./T
         Tpnm1=1.0_WP
         do n=2,5
            Tpnm1=Tpnm1*T     ! =T.^(n-1)
            th(n)=Tpnm1/float(n)     ! =T.^(n-1) ./ n
         end do
         do k=1,ns
            if(T<thermo(k,1)) then
               hort(k)=dot_product(thermo(k,2:7),th)  ! coefficients in lower temperature range
            else
               hort(k)=dot_product(thermo(k,9:14),th) ! coefficients in upper temperature range
            endif
         end do
      end subroutine get_hort


      !> Get normalized Gibbs functions at temperature T
      subroutine get_gort(this,ns,T,p,thermo,isGas,gort)
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns
         real(kind(1.0_WP)), intent(in) :: T,p,thermo(ns,2*ncof+1),isGas(ns)
         real(kind(1.0_WP)), intent(out) :: gort(ns)
         ! input:
         !   isGas -1 if gas,0 if liquid
         !   T     -temperature (K)
         !   p     -pressure (atm)
         !   thermo-thermo data for all species
         ! output:
         !   gort  -g_j/(RT) -normalized Gibbs functions
         ! S. B. Pope 9/26/02
         real(kind(1.0_WP)) :: tc(ncof),th(ncof),ts(ncof),tg(ncof)
         integer :: k,n
         if(ns.le.0) return
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
            if(T<thermo(k,1)) then
               gort(k)=dot_product(thermo(k,2:8),tg)  ! coefficients in lower temperature range
            else
               gort(k)=dot_product(thermo(k,9:15),tg) ! coefficients in upper temperature range
            endif
         end do
         gort=gort+isGas*log(p)
      end subroutine get_gort


      !> Determine temperature given enthalpy
      subroutine hor2T(this,ns,z,hin,T_low,T_high,thermo,T,iret)
         implicit none
         class(chem_state),  intent(in)  :: this
         integer,            intent(in)  :: ns
         real(kind(1.0_WP)), intent(in)  :: z(ns),hin,T_low,T_high,thermo(ns,2*ncof+1)
         real(kind(1.0_WP)), intent(out) :: T
         integer,            intent(out) :: iret
         ! input:
         !	ns		- number of species
         !   z      -moles of species
         !   hin    -enthalpy/R (K)=z'*h
         !   T_low  -lower bound on temperature range
         !   T_high  _ upper bound on temperature range
         !   thermo -thermo data
         ! output:
         !   T   -temperature (K)
         !   iret= 0  for success
         !       =-1  for T > T_high
         !       =-2  for T < T_low
         !       =-3  for failure

         ! Notes:  if the temperature is outside the range [T_low T_high]
         !   then T is returned as the closest of these bounds.
         !   If iteration fails,T is returned as T=-1.

         ! S. B. Pope 9/26/02

         integer :: itmax,it
         real(kind(1.0_WP)) :: T_tol,T0,hort(ns),hor,h_a,T_a,h_b,T_b,dT,&
            cpor(ns),hh,cpp

         itmax=100     ! maximum number of Newton iterations (usually only 3 required)
         T_tol=1e-6    ! error tolerance
         T0=1500.0_WP  ! initial guess
         iret=0        ! anticipate success

         !  determine if T>T0 and bracket T in [T_a T_b]
         call this%get_hort(ns,T0,thermo,hort)
         hor=dot_product(z,hort)*T0

         if(hin>hor) then	! T > T0=T_a
            h_a=hor
            T_a=T0
            call this%get_hort(ns,T_high,thermo,hort)
            h_b=dot_product(z,hort)*T_high
            if(hin.ge.h_b) then
               T=T_high   ! T > T_high (return T=T_high)
               iret=-1
               return
            endif
            T_b=T_high
         else
            h_b=hor	! T < T0=T_b
            T_b=T0
            call this%get_hort(ns,T_low,thermo,hort)
            h_a=dot_product(z,hort)*T_low
            if(hin.le.h_a) then
               T=T_low    ! T < T_low (return T=T_low)
               iret=-2
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
            if(abs(dT).lt.T_tol) return  ! success
         end do

         ! failure
         iret=-3
         T=-1.0_WP

      end subroutine hor2T


      !> Generate (possibly) perturbed CE problem
      subroutine perturb(this,ns,nsd,nsu,ne,ned,neu,nrc,Nd,cr,BR,E,ifop,lud,eps_el,eps_sp,pert_tol,pert_skip,zdp,zup,zupper,crp,npert,max_pert,iret)
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns,nsd,nsu,ne,ned,neu,nrc,ifop,lud,pert_skip
         real(kind(1.0_WP)), intent(in) :: Nd(nsd),cr(nrc),BR(nsu,nrc),E(ns,ne),eps_el,eps_sp,pert_tol
         integer, intent(out) :: npert,iret
         real(kind(1.0_WP)), intent(out) :: zdp(nsd),zup(nsu),zupper(nsu),crp(nrc),max_pert

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
         !   lud    -logical unit for diagnostic output
         !   eps_el -relative lower bound on element moles
         !   eps_sp -relative lower bound on species moles
         !   pert_tol- largest allowed perturbation (moles/moles of atoms)
         !  pert_skip>0 to skip perturbing undetermined species

         !  Output:
         !   zdp    -perturbed moles of determined species
         !   zup    -min-max solution for undetermined species
         !   zupper -upper bound on undetermined species moles
         !   crp    -perturbed reduced constraint vector
         !   npert  -number of perturbations made
         !   max_pert- largest normalized perturbation made
         !   iret   = 0  successful operation
         !          =-1,non-realizable large perturbation made
         !          =-2,failed to determine max-min composition
         !          =-3,zero atoms

         integer :: j,k
         real(kind(1.0_WP)) :: zdlim,zatoms,&
         cre(neu),sumcre,cref,zed(ned),zeu(neu),ze(ne),zeu_in(neu),zau,zelow,&
         zemax,zumm(nsu),zumin,zulow,zeumax

         iret=0         ! anticipate success
         npert=0        ! number of perturbations
         max_pert=0.0_WP  ! largest normalized perturbation

         if(ifop.ge.5) then
            write(lud,*)'perturb'
            write(lud,'(a,9i4)')'ne ned neu ns nsd nsu nrc= ',ne,ned,neu,ns,nsd,nsu,nrc
         endif

         zdp=Nd       ! check that determined species are non-negative
         zatoms=0.0_WP  ! estimate of moles of atoms
         if(nsd.gt.0) then
         do j=1,ne
            zatoms=zatoms+abs(dot_product(Nd,E(1:nsd,j)))
         end do
         endif

         do j=1,neu
         zatoms=zatoms+abs(cr(j))
         end do

         if(zatoms.le.0.0_WP) then
            if(ifop .ge.1) write(lud,*) 'perturb: no atoms'
            iret=-3
            return
         endif

         zdlim=zatoms*pert_tol
         do k=1,nsd
            if(zdp(k)<0.0_WP) then
               if(abs(zdp(k))>zdlim  ) then ! significantly negative
                     max_pert=max(max_pert,abs(zdp(k))/zatoms)
                     npert=npert+1
                     if(ifop.ge.1) write(lud,*)'perturb: negative determined species '
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
            if(cre(k)<0.0_WP) then
               if(abs(cre(k))>cref ) then ! significantly negative
                     max_pert=max(max_pert,abs(cre(k))/zatoms)
                     npert=npert+1
                     if(ifop.ge.1) write(lud,*)'perturb: negative undetermined element '
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
            if(zeu(j)<zelow) then
               npert=npert+1
               max_pert=max(max_pert,(zelow-zeu(j))/zatoms)
               zeu(j)=zelow
            endif
         end do

         if(ifop.ge.5.and.npert>0) then
            write(lud,*) 'zeu_in'
            write(lud,'(1p,5e13.4)') zeu_in
            write(lud,*) 'zeu'
            write(lud,'(1p,5e13.4)') zeu
            write(lud,*) 'zeu-zeu_in'
            write(lud,'(1p,5e13.4)') zeu-zeu_in
         endif

         zupper=0.0_WP    ! determine upper bound on undetermined species
         do j=1,nsu
            zupper(j)=1/maxval(E(nsd+j,ned+1:ne)/zeu)
         end do

         if(ifop.ge.5) then
            write(lud,*) 'log10(zupper)= '
            write(lud,'(1p,5e13.4)') log10(zupper)
         endif

         !  determine max-min moles of undetermined species

         call this%maxmin_comp(nsu,nrc,BR,crp,zumm,zumin,iret)

         if(iret<0) then
            if(ifop.ge.1) write(lud,*)'perturb: maxmin_comp failed,iret=',iret
            iret=-2 
            return 
         endif

         zup=zumm
         do j=1,nsu     ! impose lower limit on undetermined species
            zulow=eps_sp*zupper(j)
            if(zumm(j)<zulow) then
               zup(j)=zulow
               npert=npert+1
               max_pert=max(max_pert,(zup(j)-zumm(j))/zatoms)
            endif
         end do

         if(pert_skip.gt.0) return  ! do not modify constraints

         !  modify constraints according to the perturbation in undetermined species
         crp=crp+matmul(zup-zumm,BR)

         if(max_pert.gt.pert_tol) then
            if(ifop.ge.1) then
               write(lud,*)'perturb: large perturbation made'
               write(lud,'(a,i5,1p,e13.4)') 'npert,maxpert=',npert,max_pert
            endif
            iret=-1 
         endif

      end subroutine perturb


      !> Get zg,the value of z which minimizes g'z,subject to B'*z=c,z(i)>=0.
      subroutine get_Nming(this,nz,nc,B,c,g,zg,iret)
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: nz,nc
         real(kind(1.0_WP)), intent(in) :: B(nz,nc),c(nc),g(nz)
         integer, intent(out) :: iret
         real(kind(1.0_WP)), intent(out) :: zg(nz)
         !  Exit flag from linprog: iret<=0 for failure.
         !  S.B. Pope 10/1/02
         integer :: i,iftest  
         real(kind(1.0_WP)) :: BT(nc,nz),errc(nc),errmax,gmin
         BT=transpose(B)
         call this%solve_linprog(nz,nc,g,BT,c,zg,iret)
         do i=1,nz  ! guard against small negative values due to round-off
         zg(i)=max(zg(i),0.0_WP)
         end do
         iftest=0	! for testing only
         if(iftest.eq.0 ) return
         ! Test satisfaction of constraints and value of g'z
         errc=matmul(zg,B)-c
         errmax=max(maxval(errc),-minval(errc))
         gmin=dot_product(zg,g)
         write(0,*)'get_Nming: iret,errmax,gmin=',iret,errmax,gmin
      end subroutine get_Nming


      !> Determine z=N0+d which satisfies:
      !>    1) z(i).ge.eps
      !>    2) B'*z=c
      !>    3) t=max_i(|d(i)|) is minimized.
      subroutine min_pert(this,nz,nc,B,c,N0,eps,z,iret)
         ! Input:
         !  nz -length of z
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
         real(kind(1.0_WP)), intent(in)  :: B(nz,nc),c(nc),N0(nz),eps
         real(kind(1.0_WP)), intent(out) :: z(nz)
         integer, intent(out)          :: iret
         real(kind(1.0_WP)) :: A(nc+nz,3*nz),x(3*nz),f(3*nz),r(nc+nz),d(nz),tp,tm,res,Bsc,csc,zsc
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
         if(linear) then  !  min. sum of dz
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
         if(iret.eq.0) then
            do i=1,nz
               z(i)=max(x(i),0.0_WP)*zsc+eps
            end do
         endif
         if(.false.) return
         !  diagnostics
         write(0,*)' '
         write(0,*)'min_pert: iret=',iret
         write(0,*)'min_pert: N0,z,N0-z'
         do i=1,nz
            write(0,'(1p,10e11.2)') N0(i),z(i),N0(i)-z(i),x(nz+i),x(2*nz+i),x(nz+i)+x(2*nz+i),x(nz+i)-x(2*nz+i)
         end do
         write(0,*) 'min_pert: min(x)=',minval(x)
         r(1:nc)=matmul(z(1:nz),B(1:nz,1:nc))-c(1:nc)
         write(0,*) 'min_pert: residual=',norm2(r(1:nc))
         write(0,*) 'min_pert: sum of perts.=',sum(x(nz+1:3*nz))
      end subroutine min_pert


      !> Get the normalized Cp's at temperature T
      subroutine get_cpor(this,ns,T,thermo,cpor)
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns
         real(kind(1.0_WP)), intent(in) :: T,thermo(ns,2*ncof+1)
         real(kind(1.0_WP)), intent(out) :: cpor(ns)
         ! input:
         !	ns	  -number of species
         !   T     -temperature (K)
         !   thermo-thermo data for all species
         ! output:
         !   cpor  -Cp_j/R   -normalized constant-pressure specific heats
         ! S. B. Pope 9/26/02
         real(kind(1.0_WP)) :: tc(5)
         integer :: k,n
         cpor=0.0_WP
         tc(1)=1.0_WP  ! coefficient multipliers for specific heats
         do n=2,5
            tc(n)=T*tc(n-1)   ! =T.^(n-1) 
         end do
         do k=1,ns
            if(T.lt.thermo(k,1)) then
               cpor(k)=dot_product(thermo(k,2:6),tc)  ! coefficients in lower temperature range
            else
               cpor(k)=dot_product(thermo(k,9:13),tc) ! coefficients in upper temperature range
            endif
         end do
      end subroutine get_cpor


      !> Determine the max-min composition.
      subroutine maxmin_comp(this,nz,nc,B,c,zm,zmin,iret)
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: nz,nc
         real(kind(1.0_WP)), intent(in) :: B(nz,nc),c(nc)
         integer, intent(out) :: iret
         real(kind(1.0_WP)), intent(out) :: zm(nz),zmin
         !  Find zm which maximizes zmin=min_i(z(i)) subject to B'*z=c.
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
         real(kind(1.0_WP)) :: A(nc,nz+1),bsum(nc),f(nz+1),x(nz+1)
         nx=nz+1
         bsum=sum(B,dim=1)
         f=0
         ! First assume zmin>0,x=[z'-zmin zmin]'
         f(nx)=-1.0_WP  ! minimize -zmin
         A(1:nc,1:nz)=transpose(B)
         A(1:nc,nx)=bsum
         call this%solve_linprog(nx,nc,f,A,c,x,iret)
         if(iret.eq.0) then  !  success,zmin>=0
            zmin=x(nx)
            zm=x(1:nz)+zmin
            return
         elseif(iret.ne.-2) then  !  failure
            return
         endif
         ! zmin<0,re-define x=[z'-zmin -zmin]'
         f(nx)=1.0_WP  ! minimize -zmin
         A(1:nc,nx)=-bsum
         call this%solve_linprog(nx,nc,f,A,c,x,iret)
         if(iret.ne.0) return  ! failure
         zmin=-x(nx)
         zm=x(1:nz)+zmin
      end subroutine maxmin_comp


      !> Determine x=xm which minimizes g=f'*x subject to x(i)>=0 and A*x=b,where A has full rank.
      subroutine solve_linprog(this,nx,nb,f,A,b,xm,iret)
         use linprog, only: lp
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: nx,nb
         real(kind(1.0_WP)), intent(in)  :: f(nx),A(nb,nx),b(nb)
         real(kind(1.0_WP)), intent(out) :: xm(nx)
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
         real(kind(1.0_WP)) :: eps=1e-9,ale(1,1),age(1,1),ble(1),bge(1)
         call lp(nx,0,0,nb,ale,age,A,ble,bge,b,f,xm,iret,toler=eps)
         if(iret<0) then
            !XXX write(0,*)'solve_linprog,iret=',iret  ! SBP XXX
         endif
      end subroutine solve_linprog


      !> Find the chemical equilibium state
      subroutine eqiulibrate(this)
         implicit none
         class(chem_state), intent(inout) :: this
      end subroutine eqiulibrate


end module chem_state_class