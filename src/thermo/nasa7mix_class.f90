!> Thermally perfect ideal-gas mixture with NASA7 polynomial thermodynamics: T-dependent cp, and enthalpies that
!> include the formation enthalpy, so the internal energy carried by a compressible solver is the total (sensible +
!> chemical) energy and heat release appears through T(rho,e,Y) without any explicit source term.
!> Generic: initialized from per-species arrays (molar mass, NASA7 switch temperature, 14 coefficients), hence no
!> dependency on a particular kinetics module. T(rho,e) is a bracketed Newton solve; the combined (rho,e)->(p,T,c)
!> accessor performs a single solve, warm-started from an optional temperature guess.
module nasa7mix_class
   use precision,      only: WP
   use string,         only: str_medium
   use messager,       only: die
   use material_class, only: material
   implicit none
   private

   public :: nasa7mix,nasa7mix_nclamp,nasa7mix_last_iter

   ! Solver diagnostics: module-level because thermodynamic accessors are intent(in) on the material
   integer, save :: nasa7mix_nclamp=0      !< Number of T(rho,e) solves clamped to the temperature bracket
   integer, save :: nasa7mix_last_iter=0   !< Newton iterations used by the last T(rho,e) solve

   type, extends(material) :: nasa7mix
      real(WP), dimension(:),   allocatable :: W      !< Molar mass [kg/mol]
      real(WP), dimension(:),   allocatable :: Rsp    !< Species gas constant Ru/W [J/(kg K)]
      real(WP), dimension(:),   allocatable :: Tmid   !< NASA7 switch temperature [K]
      real(WP), dimension(:,:), allocatable :: alo    !< (7,ns) low-temperature coefficients a1..a7
      real(WP), dimension(:,:), allocatable :: ahi    !< (7,ns) high-temperature coefficients a1..a7
      real(WP) :: Ru   =8.314462618_WP   !< Universal gas constant [J/(mol K)]
      real(WP) :: p_ref=101325.0_WP      !< NASA7 standard-state pressure [Pa]
      real(WP) :: T_lo =200.0_WP         !< Lower end of the T(rho,e) bracket [K]
      real(WP) :: T_hi =5000.0_WP        !< Upper end of the T(rho,e) bracket [K]
      real(WP) :: T_ref=300.0_WP         !< Linearization point for the cold-start guess [K]
      real(WP) :: tol_T=1.0e-10_WP       !< Relative Newton tolerance on T
      integer  :: max_iter=50            !< Newton iteration limit
      logical  :: clamp_T=.true.         !< Internal energy outside the bracket: clamp T (true) or die (false)
   contains
      procedure, private :: nasa7mix_initialize
      generic   :: initialize              => nasa7mix_initialize
      procedure :: species_index
      procedure :: get_p_from_rho_e        => nasa7mix_get_p_from_rho_e
      procedure :: get_T_from_p_rho        => nasa7mix_get_T_from_p_rho
      procedure :: get_c_from_p_rho        => nasa7mix_get_c_from_p_rho
      procedure :: get_e_from_p_rho        => nasa7mix_get_e_from_p_rho
      procedure :: get_e_from_p_T          => nasa7mix_get_e_from_p_T
      procedure :: get_p_from_rho_T        => nasa7mix_get_p_from_rho_T
      procedure :: get_rho_from_p_T        => nasa7mix_get_rho_from_p_T
      procedure :: get_cv_from_rho_T       => nasa7mix_get_cv_from_rho_T
      procedure :: get_h_from_p_T          => nasa7mix_get_h_from_p_T
      procedure :: get_hk_from_p_T         => nasa7mix_get_hk_from_p_T
      procedure :: get_s_from_p_T          => nasa7mix_get_s_from_p_T
      procedure :: get_g_from_p_T          => nasa7mix_get_g_from_p_T
      procedure :: get_gruneisen_from_rho_e=> nasa7mix_get_gruneisen_from_rho_e
      procedure :: get_T_from_rho_e        => nasa7mix_get_T_from_rho_e
      procedure :: get_c_from_rho_e        => nasa7mix_get_c_from_rho_e
      procedure :: get_cv_from_rho_e       => nasa7mix_get_cv_from_rho_e
      procedure :: get_rhoe_from_p_rho     => nasa7mix_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T       => nasa7mix_get_rhoe_from_p_T
      procedure :: get_pTc_from_rho_e      => nasa7mix_get_pTc_from_rho_e
      procedure :: print                   => nasa7mix_print
      procedure :: finalize                => nasa7mix_finalize
      procedure, private :: poly           !< Dimensionless NASA7 cp/R, h/(RT), s/R of one species
      procedure, private :: get_mix        !< Mixture gas constant, cp and h at T
      procedure, private :: solve_T        !< Bracketed Newton for T(e,y)
   end type nasa7mix

contains

   !> Initialize from per-species arrays (ns inferred from size(W)); coeffs(ns,14) holds the low-T a1..a7 in
   !> columns 1-7 and the high-T a1..a7 in columns 8-14 (the layout emitted by yaml2nga.py)
   subroutine nasa7mix_initialize(this,W,T_mid,coeffs,species_names,name,Ru,p_ref,T_lo,T_hi)
      class(nasa7mix), intent(inout) :: this
      real(WP), dimension(:),   intent(in) :: W,T_mid
      real(WP), dimension(:,:), intent(in) :: coeffs
      character(len=*), dimension(:), intent(in) :: species_names
      character(len=*), intent(in), optional :: name
      real(WP), intent(in), optional :: Ru,p_ref,T_lo,T_hi
      integer :: ns,k
      ns=size(W)
      if (size(T_mid).ne.ns) call die('[nasa7mix initialize] T_mid size mismatch')
      if (size(coeffs,1).ne.ns.or.size(coeffs,2).ne.14) call die('[nasa7mix initialize] coeffs must be dimensioned (ns,14)')
      if (size(species_names).ne.ns) call die('[nasa7mix initialize] species_names size mismatch')
      if (any(W.le.0.0_WP)) call die('[nasa7mix initialize] non-positive molar mass')
      if (present(name)) this%name=name
      if (present(Ru)) this%Ru=Ru
      if (present(p_ref)) this%p_ref=p_ref
      if (present(T_lo)) this%T_lo=T_lo
      if (present(T_hi)) this%T_hi=T_hi
      if (this%T_lo.le.0.0_WP.or.this%T_hi.le.this%T_lo) call die('[nasa7mix initialize] invalid temperature bracket')
      if (this%T_ref.lt.this%T_lo.or.this%T_ref.gt.this%T_hi) this%T_ref=0.5_WP*(this%T_lo+this%T_hi)
      call this%finalize()
      this%ns=ns
      allocate(this%W(ns),this%Rsp(ns),this%Tmid(ns),this%alo(7,ns),this%ahi(7,ns),this%species_names(ns))
      this%W=W
      this%Rsp=this%Ru/W
      this%Tmid=T_mid
      do k=1,ns
         this%alo(:,k)=coeffs(k,1:7)
         this%ahi(:,k)=coeffs(k,8:14)
         this%species_names(k)=trim(adjustl(species_names(k)))
      end do
   end subroutine nasa7mix_initialize

   !> Index of a species by name (0 if absent)
   integer function species_index(this,name) result(k)
      class(nasa7mix), intent(in) :: this
      character(len=*), intent(in) :: name
      integer :: n
      k=0
      do n=1,this%ns
         if (trim(this%species_names(n)).eq.trim(adjustl(name))) then; k=n; return; end if
      end do
   end function species_index

   !> NASA7 polynomials of species k at T (lnT passed in to share the log across species)
   subroutine poly(this,T,lnT,k,cp_R,h_RT,s_R)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: T,lnT
      integer, intent(in) :: k
      real(WP), intent(out) :: cp_R,h_RT,s_R
      real(WP), dimension(7) :: a
      real(WP) :: T2,T3,T4,Tinv
      if (T.le.this%Tmid(k)) then; a=this%alo(:,k); else; a=this%ahi(:,k); end if
      T2=T*T; T3=T2*T; T4=T3*T; Tinv=1.0_WP/T
      cp_R=a(1)+a(2)*T+a(3)*T2+a(4)*T3+a(5)*T4
      h_RT=a(1)+0.5_WP*a(2)*T+a(3)*T2/3.0_WP+0.25_WP*a(4)*T3+0.2_WP*a(5)*T4+a(6)*Tinv
      s_R =a(1)*lnT+a(2)*T+0.5_WP*a(3)*T2+a(4)*T3/3.0_WP+0.25_WP*a(5)*T4+a(7)
   end subroutine poly

   !> Mixture gas constant Rm, cp and h (mass units) at T; one polynomial pass over the species
   subroutine get_mix(this,T,y,Rm,cp,h)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: T
      real(WP), dimension(:), intent(in) :: y
      real(WP), intent(out) :: Rm,cp,h
      real(WP) :: cp_R,h_RT,s_R,yR
      integer :: k
      Rm=0.0_WP; cp=0.0_WP; h=0.0_WP
      do k=1,this%ns
         call this%poly(T,0.0_WP,k,cp_R,h_RT,s_R)
         yR=y(k)*this%Rsp(k)
         Rm=Rm+yR
         cp=cp+yR*cp_R
         h =h +yR*T*h_RT
      end do
   end subroutine get_mix

   !> Bracketed Newton for the temperature such that e_mix(T,y)=e; returns cp and Rm at the solution.
   !> Convergence is checked on the temperature increment (|e| is dominated by formation enthalpies).
   !> Out-of-range e: the bracket end is tested at most once, then T is clamped there (clamp_T) or the run dies.
   subroutine solve_T(this,e,y,T,cp,Rm,Tguess)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: e
      real(WP), dimension(:), intent(in) :: y
      real(WP), intent(out) :: T,cp,Rm
      real(WP), intent(in), optional :: Tguess
      real(WP) :: Ta,Tb,h,f,cv,dT,Tn
      integer :: it
      logical :: have_guess,lo_ok,hi_ok,clamped
      ! Initial guess: warm start when a guess inside the bracket is provided, else linearization at T_ref
      have_guess=.false.
      if (present(Tguess)) have_guess=(Tguess.gt.this%T_lo.and.Tguess.lt.this%T_hi)
      if (have_guess) then
         T=Tguess
      else
         call this%get_mix(this%T_ref,y,Rm,cp,h)
         T=this%T_ref+(e-(h-Rm*this%T_ref))/max(cp-Rm,1.0_WP)
         T=min(max(T,this%T_lo),this%T_hi)
      end if
      ! Newton iteration with bisection safeguard
      Ta=this%T_lo; Tb=this%T_hi; lo_ok=.false.; hi_ok=.false.; clamped=.false.
      do it=1,this%max_iter
         call this%get_mix(T,y,Rm,cp,h)
         f=h-Rm*T-e; cv=cp-Rm
         if (f.gt.0.0_WP) then; Tb=T; hi_ok=.true.; else; Ta=T; lo_ok=.true.; end if
         dT=-f/max(cv,1.0_WP)
         if (abs(dT).le.this%tol_T*T) exit
         Tn=T+dT
         if (Tn.le.Ta) then
            if (.not.lo_ok) then
               ! Newton wants to go below the bracket: check whether e is reachable at all
               call this%get_mix(Ta,y,Rm,cp,h)
               if (h-Rm*Ta.ge.e) then; T=Ta; clamped=.true.; exit; end if
               lo_ok=.true.
            end if
            Tn=0.5_WP*(Ta+Tb)
         else if (Tn.ge.Tb) then
            if (.not.hi_ok) then
               call this%get_mix(Tb,y,Rm,cp,h)
               if (h-Rm*Tb.le.e) then; T=Tb; clamped=.true.; exit; end if
               hi_ok=.true.
            end if
            Tn=0.5_WP*(Ta+Tb)
         end if
         T=Tn
      end do
      nasa7mix_last_iter=min(it,this%max_iter)
      if (it.gt.this%max_iter) call fail('Newton iteration did not converge (non-finite input?)')
      if (clamped) then
         if (.not.this%clamp_T) call fail('internal energy outside the temperature bracket')
         nasa7mix_nclamp=nasa7mix_nclamp+1
      end if
   contains
      subroutine fail(why)
         use string, only: str_long
         character(len=*), intent(in) :: why
         character(len=str_long) :: msg
         write(msg,'(a,a,a,es12.5,a,es12.5,a,es12.5,a,es12.5,a,i0)') '[nasa7mix solve_T] ',why,': e=',e,' T=',T, &
         &  ' sum(y)=',sum(y(1:this%ns)),' min(y)=',minval(y(1:this%ns)),' iterations=',it
         call die(trim(msg))
      end subroutine fail
   end subroutine solve_T

   ! ---------------------------------------------------------------------------------------------------------------
   ! (rho,e) family: one Newton solve each; use get_pTc_from_rho_e when several primitives are needed at once
   ! ---------------------------------------------------------------------------------------------------------------

   real(WP) function nasa7mix_get_p_from_rho_e(this,rho,e,y) result(p)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: T,cp,Rm
      call this%solve_T(e,y,T,cp,Rm)
      p=rho*Rm*T
   end function nasa7mix_get_p_from_rho_e

   real(WP) function nasa7mix_get_T_from_rho_e(this,rho,e,y) result(T)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cp,Rm
      call this%solve_T(e,y,T,cp,Rm)
   end function nasa7mix_get_T_from_rho_e

   real(WP) function nasa7mix_get_c_from_rho_e(this,rho,e,y) result(c)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: T,cp,Rm
      call this%solve_T(e,y,T,cp,Rm)
      c=sqrt(cp/(cp-Rm)*Rm*T)
   end function nasa7mix_get_c_from_rho_e

   real(WP) function nasa7mix_get_cv_from_rho_e(this,rho,e,y) result(cv)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: T,cp,Rm
      call this%solve_T(e,y,T,cp,Rm)
      cv=cp-Rm
   end function nasa7mix_get_cv_from_rho_e

   !> Gruneisen parameter (dp/de)_rho/rho = Rm/cv for an ideal gas
   real(WP) function nasa7mix_get_gruneisen_from_rho_e(this,rho,e,y) result(gruneisen)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: T,cp,Rm
      call this%solve_T(e,y,T,cp,Rm)
      gruneisen=Rm/(cp-Rm)
   end function nasa7mix_get_gruneisen_from_rho_e

   !> Combined flash: a single Newton solve warm-started from Tguess provides p, T and c
   subroutine nasa7mix_get_pTc_from_rho_e(this,rho,e,y,p,T,c,Tguess)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP), intent(out) :: p,T,c
      real(WP), intent(in), optional :: Tguess
      real(WP) :: cp,Rm
      call this%solve_T(e,y,T,cp,Rm,Tguess)
      p=rho*Rm*T
      c=sqrt(cp/(cp-Rm)*Rm*T)
   end subroutine nasa7mix_get_pTc_from_rho_e

   ! ---------------------------------------------------------------------------------------------------------------
   ! (p,rho), (p,T), (rho,T) families: closed form or one polynomial pass
   ! ---------------------------------------------------------------------------------------------------------------

   real(WP) function nasa7mix_get_T_from_p_rho(this,p,rho,y) result(T)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      T=p/(rho*sum(y(1:this%ns)*this%Rsp(1:this%ns)))
   end function nasa7mix_get_T_from_p_rho

   real(WP) function nasa7mix_get_c_from_p_rho(this,p,rho,y) result(c)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: T,Rm,cp,h
      T=this%get_T_from_p_rho(p,rho,y)
      call this%get_mix(T,y,Rm,cp,h)
      c=sqrt(cp/(cp-Rm)*p/rho)
   end function nasa7mix_get_c_from_p_rho

   real(WP) function nasa7mix_get_e_from_p_rho(this,p,rho,y) result(e)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: T,Rm,cp,h
      T=this%get_T_from_p_rho(p,rho,y)
      call this%get_mix(T,y,Rm,cp,h)
      e=h-Rm*T
   end function nasa7mix_get_e_from_p_rho

   real(WP) function nasa7mix_get_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      rhoe=rho*this%get_e_from_p_rho(p,rho,y)
   end function nasa7mix_get_rhoe_from_p_rho

   real(WP) function nasa7mix_get_e_from_p_T(this,p,T,y) result(e)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm,cp,h
      call this%get_mix(T,y,Rm,cp,h)
      e=h-Rm*T
   end function nasa7mix_get_e_from_p_T

   real(WP) function nasa7mix_get_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm,cp,h
      call this%get_mix(T,y,Rm,cp,h)
      rhoe=p/(Rm*T)*(h-Rm*T)
   end function nasa7mix_get_rhoe_from_p_T

   real(WP) function nasa7mix_get_p_from_rho_T(this,rho,T,y) result(p)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      p=rho*sum(y(1:this%ns)*this%Rsp(1:this%ns))*T
   end function nasa7mix_get_p_from_rho_T

   real(WP) function nasa7mix_get_rho_from_p_T(this,p,T,y) result(rho)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      rho=p/(sum(y(1:this%ns)*this%Rsp(1:this%ns))*T)
   end function nasa7mix_get_rho_from_p_T

   real(WP) function nasa7mix_get_cv_from_rho_T(this,rho,T,y) result(cv)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm,cp,h
      call this%get_mix(T,y,Rm,cp,h)
      cv=cp-Rm
   end function nasa7mix_get_cv_from_rho_T

   real(WP) function nasa7mix_get_h_from_p_T(this,p,T,y) result(h)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm,cp
      call this%get_mix(T,y,Rm,cp,h)
   end function nasa7mix_get_h_from_p_T

   !> Species specific enthalpies h_k(T) [J/kg], formation enthalpy included
   subroutine nasa7mix_get_hk_from_p_T(this,p,T,y,hk)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(:), intent(out) :: hk
      real(WP) :: cp_R,h_RT,s_R
      integer :: k
      do k=1,this%ns
         call this%poly(T,0.0_WP,k,cp_R,h_RT,s_R)
         hk(k)=this%Rsp(k)*T*h_RT
      end do
   end subroutine nasa7mix_get_hk_from_p_T

   !> Mixture entropy: sum_k y_k [s0_k(T) - R_k ln(X_k p/p_ref)] with mole-fraction partial pressures
   real(WP) function nasa7mix_get_s_from_p_T(this,p,T,y) result(s)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cp_R,h_RT,s_R,lnT,xsum,Xk
      integer :: k
      xsum=sum(y(1:this%ns)/this%W(1:this%ns))
      if (xsum.le.tiny(1.0_WP)) call die('[nasa7mix get_s_from_p_T] non-positive mole sum')
      lnT=log(T)
      s=0.0_WP
      do k=1,this%ns
         call this%poly(T,lnT,k,cp_R,h_RT,s_R)
         Xk=y(k)/this%W(k)/xsum
         s=s+y(k)*this%Rsp(k)*(s_R-log(max(Xk*p/this%p_ref,tiny(1.0_WP))))
      end do
   end function nasa7mix_get_s_from_p_T

   real(WP) function nasa7mix_get_g_from_p_T(this,p,T,y) result(g)
      class(nasa7mix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      g=this%get_h_from_p_T(p,T,y)-T*this%get_s_from_p_T(p,T,y)
   end function nasa7mix_get_g_from_p_T

   subroutine nasa7mix_print(this)
      use messager, only: log
      use string,   only: str_long
      class(nasa7mix), intent(in) :: this
      character(len=str_long) :: msg
      integer :: k
      write(msg,'(a,a)') '[material:nasa7mix] ',trim(this%name); call log(msg)
      write(msg,'(2x,a,i0)') 'ns    = ',this%ns; call log(msg)
      write(msg,'(2x,a,es12.5,a,es12.5)') 'Ru    = ',this%Ru,'  p_ref = ',this%p_ref; call log(msg)
      write(msg,'(2x,a,es12.5,a,es12.5)') 'T_lo  = ',this%T_lo,'  T_hi  = ',this%T_hi; call log(msg)
      do k=1,this%ns
         write(msg,'(2x,a,i0,a,a,a,es12.5,a,es12.5)') 'species(',k,') = ',trim(this%species_names(k)), &
         &  '  W=',this%W(k),'  T_mid=',this%Tmid(k)
         call log(msg)
      end do
   end subroutine nasa7mix_print

   subroutine nasa7mix_finalize(this)
      class(nasa7mix), intent(inout) :: this
      if (allocated(this%W))             deallocate(this%W)
      if (allocated(this%Rsp))           deallocate(this%Rsp)
      if (allocated(this%Tmid))          deallocate(this%Tmid)
      if (allocated(this%alo))           deallocate(this%alo)
      if (allocated(this%ahi))           deallocate(this%ahi)
      if (allocated(this%species_names)) deallocate(this%species_names)
      this%ns=0
   end subroutine nasa7mix_finalize

end module nasa7mix_class
