!> Simplified transport model for finite-rate chemistry on top of the generated fcmech module.
!> Mixture viscosity from a power law or Sutherland's law (or, for validation against Cantera, Wilke's rule on
!> the fcmech pure-species values), thermal conductivity from a constant Prandtl number (or the Mathur-Saxena
!> mixture rule), and a single species mass diffusivity rho*D from a constant Lewis or Schmidt number.
!> Turbulent Prandtl and Schmidt numbers are carried for the LES closure applied by the solver glue.
module fctransport_class
   use precision, only: WP
   use string,    only: str_medium
   use messager,  only: die
   use fcmech,    only: nS,W_sp,Rcst,fcmech_nasa7,fcmech_get_viscosity,fcmech_get_conductivity
   implicit none
   private

   public :: fctransport,tr_powerlaw,tr_sutherland,tr_mixavg

   integer, parameter :: tr_powerlaw  =1   !< mu=mu_ref*(T/T_ref)**alpha, lambda=mu*cp/Pr
   integer, parameter :: tr_sutherland=2   !< mu=mu_ref*(T/T_ref)**1.5*(T_ref+S)/(T+S), lambda=mu*cp/Pr
   integer, parameter :: tr_mixavg    =3   !< Wilke viscosity and Mathur-Saxena conductivity from pure-species data

   type :: fctransport
      character(len=str_medium) :: name='UNNAMED_FCTRANSPORT'
      integer  :: model=tr_powerlaw
      real(WP) :: mu_ref=1.716e-5_WP   !< Reference viscosity [Pa s] (air at 273.15 K)
      real(WP) :: T_ref =273.15_WP     !< Reference temperature [K]
      real(WP) :: alpha =0.7_WP        !< Power-law exponent
      real(WP) :: S_suth=110.4_WP      !< Sutherland temperature [K] (air)
      real(WP) :: Pr    =0.7_WP        !< Prandtl number
      real(WP) :: Le    =1.0_WP        !< Lewis number (used when Sc<=0): rho*D=lambda/(cp*Le)
      real(WP) :: Sc    =-1.0_WP       !< Schmidt number (>0 overrides Le): rho*D=mu/Sc
      real(WP) :: Pr_t  =0.7_WP        !< Turbulent Prandtl number
      real(WP) :: Sc_t  =0.7_WP        !< Turbulent Schmidt number
      real(WP), dimension(:,:), allocatable :: wilke_A   !< (W_j/W_i)**0.25
      real(WP), dimension(:,:), allocatable :: wilke_B   !< 1/sqrt(8*(1+W_i/W_j))
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_powerlaw
      procedure :: set_sutherland
      procedure :: set_mixavg
      procedure :: set_numbers
      procedure :: get_cp
      procedure :: get_hk
      procedure :: get_props
      procedure :: print
   end type fctransport

contains

   subroutine initialize(this,name)
      class(fctransport), intent(inout) :: this
      character(len=*), intent(in), optional :: name
      if (present(name)) this%name=trim(name)
   end subroutine initialize

   subroutine finalize(this)
      class(fctransport), intent(inout) :: this
      if (allocated(this%wilke_A)) deallocate(this%wilke_A)
      if (allocated(this%wilke_B)) deallocate(this%wilke_B)
   end subroutine finalize

   !> Power-law viscosity mu=mu_ref*(T/T_ref)**alpha
   subroutine set_powerlaw(this,mu_ref,T_ref,alpha)
      class(fctransport), intent(inout) :: this
      real(WP), intent(in), optional :: mu_ref,T_ref,alpha
      this%model=tr_powerlaw
      if (present(mu_ref)) this%mu_ref=mu_ref
      if (present(T_ref))  this%T_ref =T_ref
      if (present(alpha))  this%alpha =alpha
      if (this%mu_ref.le.0.0_WP.or.this%T_ref.le.0.0_WP) call die('[fctransport set_powerlaw] mu_ref and T_ref must be positive')
   end subroutine set_powerlaw

   !> Sutherland viscosity mu=mu_ref*(T/T_ref)**1.5*(T_ref+S)/(T+S)
   subroutine set_sutherland(this,mu_ref,T_ref,S)
      class(fctransport), intent(inout) :: this
      real(WP), intent(in), optional :: mu_ref,T_ref,S
      this%model=tr_sutherland
      if (present(mu_ref)) this%mu_ref=mu_ref
      if (present(T_ref))  this%T_ref =T_ref
      if (present(S))      this%S_suth=S
      if (this%mu_ref.le.0.0_WP.or.this%T_ref.le.0.0_WP) call die('[fctransport set_sutherland] mu_ref and T_ref must be positive')
   end subroutine set_sutherland

   !> Mixture-averaged mode: Wilke viscosity, Mathur-Saxena conductivity; composition-independent factors precomputed
   subroutine set_mixavg(this)
      class(fctransport), intent(inout) :: this
      integer :: i,j
      this%model=tr_mixavg
      if (.not.allocated(this%wilke_A)) allocate(this%wilke_A(nS,nS),this%wilke_B(nS,nS))
      do j=1,nS
         do i=1,nS
            this%wilke_A(i,j)=(W_sp(j)/W_sp(i))**0.25_WP
            this%wilke_B(i,j)=1.0_WP/sqrt(8.0_WP*(1.0_WP+W_sp(i)/W_sp(j)))
         end do
      end do
   end subroutine set_mixavg

   !> Dimensionless numbers (all optional)
   subroutine set_numbers(this,Pr,Le,Sc,Pr_t,Sc_t)
      class(fctransport), intent(inout) :: this
      real(WP), intent(in), optional :: Pr,Le,Sc,Pr_t,Sc_t
      if (present(Pr))   this%Pr  =Pr
      if (present(Le))   this%Le  =Le
      if (present(Sc))   this%Sc  =Sc
      if (present(Pr_t)) this%Pr_t=Pr_t
      if (present(Sc_t)) this%Sc_t=Sc_t
      if (this%Pr.le.0.0_WP.or.this%Le.le.0.0_WP.or.this%Pr_t.le.0.0_WP.or.this%Sc_t.le.0.0_WP) then
         call die('[fctransport set_numbers] Pr, Le, Pr_t and Sc_t must be positive')
      end if
   end subroutine set_numbers

   !> Mixture heat capacity at constant pressure [J/(kg K)]
   real(WP) function get_cp(this,T,Y) result(cp)
      class(fctransport), intent(in) :: this
      real(WP), intent(in) :: T
      real(WP), dimension(nS), intent(in) :: Y
      real(WP), dimension(nS) :: cp_R,h_RT,s_R
      call fcmech_nasa7(T,cp_R,h_RT,s_R)
      cp=Rcst*sum(Y(:)*cp_R(:)/W_sp(:))
   end function get_cp

   !> Species enthalpies h_k(T) [J/kg], formation enthalpy included
   subroutine get_hk(this,T,hk)
      class(fctransport), intent(in) :: this
      real(WP), intent(in) :: T
      real(WP), dimension(nS), intent(out) :: hk
      real(WP), dimension(nS) :: cp_R,h_RT,s_R
      call fcmech_nasa7(T,cp_R,h_RT,s_R)
      hk(:)=Rcst*T*h_RT(:)/W_sp(:)
   end subroutine get_hk

   !> Mixture cp [J/(kg K)], viscosity [Pa s], thermal conductivity [W/(m K)], and species mass diffusivity
   !> rho*D [kg/(m s)] at (T,Y)
   subroutine get_props(this,T,Y,cp,mu,lambda,rhoD)
      class(fctransport), intent(in) :: this
      real(WP), intent(in) :: T
      real(WP), dimension(nS), intent(in) :: Y
      real(WP), intent(out) :: cp,mu,lambda,rhoD
      real(WP), dimension(nS) :: mu_k,lam_k,X,sqmu
      real(WP) :: xsum,phisum
      integer :: i,j
      cp=this%get_cp(T,Y)
      select case (this%model)
      case (tr_powerlaw)
         mu=this%mu_ref*(T/this%T_ref)**this%alpha
         lambda=mu*cp/this%Pr
      case (tr_sutherland)
         mu=this%mu_ref*(T/this%T_ref)**1.5_WP*(this%T_ref+this%S_suth)/(T+this%S_suth)
         lambda=mu*cp/this%Pr
      case (tr_mixavg)
         if (.not.allocated(this%wilke_A)) call die('[fctransport get_props] call set_mixavg before using the mixavg model')
         xsum=sum(Y(:)/W_sp(:))
         X(:)=Y(:)/W_sp(:)/max(xsum,tiny(1.0_WP))
         call fcmech_get_viscosity(mu_k,T)
         call fcmech_get_conductivity(lam_k,T,mu_k)
         sqmu(:)=sqrt(mu_k(:))
         ! Wilke: mu=sum_i X_i mu_i/sum_j X_j Phi_ij, Phi_ij=(1+sqrt(mu_i/mu_j)*(W_j/W_i)^0.25)^2/sqrt(8(1+W_i/W_j))
         mu=0.0_WP
         do i=1,nS
            if (X(i).le.0.0_WP) cycle
            phisum=0.0_WP
            do j=1,nS
               phisum=phisum+X(j)*(1.0_WP+sqmu(i)/sqmu(j)*this%wilke_A(i,j))**2*this%wilke_B(i,j)
            end do
            mu=mu+X(i)*mu_k(i)/phisum
         end do
         ! Mathur-Saxena: lambda=0.5*(sum_k X_k lambda_k + 1/sum_k X_k/lambda_k)
         lambda=0.5_WP*(sum(X(:)*lam_k(:))+1.0_WP/sum(X(:)/lam_k(:)))
      case default
         call die('[fctransport get_props] unknown transport model')
      end select
      if (this%Sc.gt.0.0_WP) then
         rhoD=mu/this%Sc
      else
         rhoD=lambda/(cp*this%Le)
      end if
   end subroutine get_props

   subroutine print(this)
      use messager, only: log
      use string,   only: str_long
      class(fctransport), intent(in) :: this
      character(len=str_long) :: msg
      write(msg,'(a,a)') '[fctransport] ',trim(this%name); call log(msg)
      select case (this%model)
      case (tr_powerlaw)
         write(msg,'(2x,a,es12.5,a,es12.5,a,es12.5)') 'power law: mu_ref=',this%mu_ref,' T_ref=',this%T_ref,' alpha=',this%alpha
      case (tr_sutherland)
         write(msg,'(2x,a,es12.5,a,es12.5,a,es12.5)') 'Sutherland: mu_ref=',this%mu_ref,' T_ref=',this%T_ref,' S=',this%S_suth
      case (tr_mixavg)
         write(msg,'(2x,a)') 'mixture-averaged: Wilke viscosity, Mathur-Saxena conductivity'
      end select
      call log(msg)
      if (this%Sc.gt.0.0_WP) then
         write(msg,'(2x,a,es12.5,a,es12.5)') 'Pr=',this%Pr,' Sc=',this%Sc
      else
         write(msg,'(2x,a,es12.5,a,es12.5)') 'Pr=',this%Pr,' Le=',this%Le
      end if
      call log(msg)
      write(msg,'(2x,a,es12.5,a,es12.5)') 'Pr_t=',this%Pr_t,' Sc_t=',this%Sc_t; call log(msg)
   end subroutine print

end module fctransport_class
