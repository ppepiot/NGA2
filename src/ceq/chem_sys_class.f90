!> Chemical system class
module chem_sys_class
   use precision, only: WP
   implicit none
   private
   
   !> Expose type/constructor/methods
   public :: chem_sys

   !> Take these from the two_phase classes *** HAS TO START FROM 1 ***
   integer, parameter, public :: Lphase=1,Gphase=2

   !> NASA-7 polynomial coefficient number
   integer, parameter, public :: ncof=7

   !> Chemical system object definition
   type :: chem_sys

      ! Indices
      integer :: np                                   !< Number of phases
      integer :: ne                                   !< Number of elements
      integer :: ned                                  !< Number of determined elements
      integer :: neu                                  !< Number of undetermined elements
      integer :: ns                                   !< Number of species
      integer :: nsd                                  !< Number of determined species    
      integer :: nsu                                  !< Number of undetermined species
      integer :: ncs                                  !< Number of constrained species
      integer :: ng                                   !< Number of general linear constraints
      integer :: nc                                   !< Number of basic constraints =ne+ncs+ng
      integer :: nrc                                  !< Number of reduced constraints
      integer :: nb                                   !< Number of independent constraints =nsd+nrc

      ! Constrained species
      integer, dimension(:), allocatable :: CS        !< Indexes of constrained species (ncs x 1)

      ! Ordering arrays
      integer, dimension(:), allocatable :: sp_order  !< The k-th ordered species is sp_order(k) (ns x 1)
      integer, dimension(:), allocatable :: el_order  !< The j-th ordered element is el_order(j) (ne x 1)

      ! Phase matrix
      real(WP), dimension(:,:), allocatable :: P      !< Phase summation matrix (ns x np)

      ! Element matrices
      real(WP), dimension(:,:), allocatable :: Ein    !< Element matrix (ns x ne).  A molecule of species k contains Ein(k,j) atoms of element j.
      real(WP), dimension(:,:), allocatable :: E      !< The re-ordered element matrix (ns x ne)

      ! Constraint matrices
      real(WP), dimension(:,:), allocatable :: Bg     !< General linear constraint matrix (ns x ng)
      real(WP), dimension(:,:), allocatable :: B      !< Basic constraint matrix (ns x nc)
      real(WP), dimension(:,:), allocatable :: BR     !< The reduced constraint matrix (nsu x nrc)
      real(WP), dimension(:,:), allocatable :: A      !< Modified constraint transformation matrix (nb x nc)
      
      ! Thermodynamic coefficients
      real(WP), dimension(:,:), allocatable :: thermo !< Thermodynamic coefficients (re-ordered) (ns x 2*ncof+1)

      ! Numerical parameters (when parameters are altered,chem_sys_param_set and chem_param_def must be updated.)
      integer  :: diag                                !< Greater than 0 for diagnostics
      real(WP) :: eps_el                              !< Relative lower bound on element moles (ceq_perturb)
      real(WP) :: eps_sp                              !< Relative lower bound on species moles (ceq_perturb)
      real(WP) :: pert_tol                            !< Largest allowed normalized perturbation
      integer  :: pert_skip                           !< Set =1 to skip perturbing undetermined species

   contains
      procedure :: initialize                         !< Object initializer
      procedure :: param_def=>chem_sys_param_def      !< Define parameters
      procedure :: param_set=>chem_sys_param_set      !< Reset parameters
      procedure :: red_con                            !< Reduce the constraints
   end type chem_sys


   contains
      

      !> Chemical system initializer
      subroutine initialize(this,np,ns,ne,ncs,ng,P,Ein,CS,Bg,thermo_in,diag)

         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.

         !  Specify a constrained equilibrium system.
      
         ! Input:
         !  np    - number of phases
         !  ns    - number of species
         !  ne    - number of elements
         !  ncs   - number of constrained species
         !  ng    - number of general linear constraints
         !  Ein   - element matrix (ns x ne).  A molecule of species k contains
         !          Ein(k,j) atoms of element j.
         !  CS    - array containing indexes of constrained species (ncs x 1)
         !  Bg    - general linear constraint matrix (ns x ng)
         !  thermo_in  - thermodynamic data for species (ns x 2*ncof+1)
         !             - see below for details
         !  P     - Phase summation matrix
      
         ! Unconstrained equilibrium:
         !  For unconstrained equilibrium calculations (i.e.,in which the only 
         !  constraint is on the elements),set ncs=0,ng=0,CS=[] and Bg=[].
      
         ! Details of thermo_in:
         !  As in Chemkin,the thermodynamic properties of each species are given
         !  in terms of 7 non-dimensional coefficients,a(1:7).  Different values 
         !  of a(1:7) are given for two temperature ranges.  T* [K] denotes the
         !  upper limit of the lower temperature range,and the lower limit of the 
         !  upper temperature range.  The contents of the array thermo_in are: 
         !  thermo_in(k,1)   =T* for species k
         !  thermo_in(k,2:8) =a(1:7) for species k in the lower temperature range
         !  thermo_in(k,9:15)=a(1:7) for species k in the upper temperature range
         !  The relevant thermodynamic variables are:
         !  T   - temperature [K]
         !  H   - molar specific enthalpy [J/kmol]
         !  S   - molar specific entropy [J/(kmol K)]
         !  R   - universal gas constant [J/(kmol K)]
         !  Then for each species,H and S are given by:
         !  H/(RT)=a(6)/T+sum_{n=1}^{n=5} a(n) T^{n-1}/n
         !  S/R=a(1)log(T)+a(7)+sum_{n=2}^{n=5} a(n) T^{n-1}/(n-1)

         use messager,  only: die
         use mathtools, only: reorder_rows
         implicit none
         class(chem_sys), intent(inout) :: this
         integer,  intent(in) :: np,ns,ne,ncs,ng,CS(ncs),diag
         real(WP), intent(in) :: Ein(ns,ne),Bg(ns,ng),thermo_in(ns,2*ncof+1),P(ns,Lphase:Gphase)
         real(WP) :: BR(ns,ne+ncs+ng),A(ne+ncs+ng,ne+ncs+ng)
         integer  :: nc,j,k
         real(WP) :: cpu_redcon0,cpu_redcon
      
         ! Define parameters
         call this%param_def()

         ! Check input
         if (np.ge.1) then
            this%np=np
         else
            call die('[chem_sys initialize] np needs to be greater than zero')
         endif

         if (ns.ge.1) then
            this%ns=ns
         else
            call die('[chem_sys initialize] ns needs to be greater than zero')
         endif

         if (ne.ge.1) then
            this%ne=ne
         else
            call die('[chem_sys initialize] ne needs to be greater than zero')
         endif
      
         if (ncs.ge.0) then
            this%ncs=ncs
            if (ncs.ge.1) then
               allocate(this%CS(ncs))
               this%CS=CS(1:ncs)
            endif
         else
            call die('[chem_sys initialize] ncs cannot be negative')
         endif
      
         if (ng.ge.0) then
            this%ng=ng
            if (ng.ge.1) then
               allocate(this%Bg(ns,ng))
               this%Bg=Bg(1:ns,1:ng)
            endif
         else
            call die('[chem_sys initialize] ng cannot be negative') 
         endif
      
         ! Allocate arrays
         nc=ne+ncs+ng
         this%nc=nc
         allocate(this%el_order(ne))
         allocate(this%sp_order(ns))
         allocate(this%Ein(ns,ne))
         allocate(this%E(ns,ne))
         allocate(this%B(ns,nc))
         allocate(this%thermo(ns,2*ncof+1))
         allocate(this%P(ns,Lphase:Gphase))
      
         this%Ein=Ein(1:ns,1:ne)
         
         ! Form basic and reduced constraint equations
         call this%red_con(BR,A,diag)
         
         allocate(this%BR(this%nsu,this%nrc))
         allocate(this%A(this%nb,nc))
      
         this%BR=BR(1:this%nsu,1:this%nrc)
         this%A =A (1:this%nb,1:nc)
      
         ! Re-order thermo
         call reorder_rows(thermo_in,this%sp_order,this%thermo)
      
         ! Re-order P
         call reorder_rows(P,this%sp_order,this%P)
      
         this%diag=diag

      end subroutine initialize


      !> Reset values of parameters in chem_sys: apart from chem_sys,all arguments are optional.
      subroutine chem_sys_param_set(this,diag,eps_el,eps_sp,pert_tol,pert_skip)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_sys), intent(inout) :: this
         integer,  intent(in), optional :: diag,pert_skip
         real(WP), intent(in), optional :: eps_el,eps_sp,pert_tol
         if (present(diag))     this%diag    =diag
         if (present(eps_el))   this%eps_el  =eps_el
         if (present(eps_sp))   this%eps_sp  =eps_sp
         if (present(pert_tol)) this%pert_tol=pert_tol
         if (present(pert_skip)) this%pert_skip=pert_skip
      end subroutine chem_sys_param_set


      !> Reset parameters in chem_sys to their default settings.
      subroutine chem_sys_param_def(this)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_sys), intent(inout) :: this
         this%diag    =1          ! >0 for diagnostics
         this%eps_el  =1e-9       ! Selative lower bound on element moles (ceq_perturb)
         this%eps_sp  =1e-9       ! Selative lower bound on species moles (ceq_perturb)
         this%pert_tol=1e-4       ! Largest allowed normalized perturbation
         this%pert_skip= 0        ! Set =1 to skip perturbing undetermined species
      end subroutine chem_sys_param_def


      !> Set the reduced constraint matrix BR,and determine the ordering of the elements and species
      subroutine red_con(this,BR,A,ifop)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         use messager,  only: die
         use mathtools, only: ind_col,reorder_rows
         implicit none
         class(chem_sys), intent(inout) :: this
         integer,  intent(in)  :: ifop
         real(WP), intent(out) :: BR(this%ns,this%nc),A(this%nc,this%nc)
         integer  :: k,info,lwork,rank_def,sp_det(this%ns),kk,ndebg,el_det(this%ne),ngi,is,js,ie,je,j,indcol(this%ns+this%ne+this%ng)
         real(WP) :: Sc(this%ns,this%ncs),S(this%ns+this%nc),U(this%ns,this%ns),VT(this%nc,this%nc),work(4*(this%ns+this%nc)*(this%ns+this%nc))
         real(WP) :: DEBg(this%ns,this%ns+this%ne+this%ng),EU(this%ns,this%ne),BS(this%ns,this%nc),PU(this%ns,this%nc),BTTPU(this%nc,this%nc),BB(this%ns,this%nc)
         real(WP), parameter :: thresh=1e-10 ! threshold for singular values and vectors
         
         ! Check input
         call check_input()

         if (.true.) then !  new test for constrained species being distinct SBP 2/6/09 

            do j=1,this%ncs-1
               do k=j+1,this%ncs
                  if (this%CS(j).eq.this%CS(k)) call die('[chem_sys red_con] CS not distinct')
               end do
            end do

            else  ! original test which inappropriately assumes ordeering

               if (this%ncs.gt.1) then
                  if (minval(this%CS(2:this%ncs)-this%CS(1:this%ncs-1)).lt.1) call die('[chem_sys red_con] CS not distinct')
               endif

         endif  !  end of distinct test

         ! Set basic constraint matrix
         this%B=0.0_WP

         ! Set constrained species matrix
         if (this%ncs.gt.0) then
            Sc=0.0_WP
            do k=1,this%ncs 
               Sc(this%CS(k),k)=1.0_WP
            end do 
            this%B(:,1:this%ncs)=Sc              ! first ncs columns - constrained species
         endif
         
         this%B(1:this%ns,this%ncs+1:this%ncs+this%ne)=this%Ein(1:this%ns,1:this%ne)      ! next ne columns - element vectors

         if (this%ng.gt.0) &
         this%B(1:this%ns,this%ne+this%ncs+1:this%nc)=this%Bg(1:this%ns,1:this%ng)        ! last ng columns - general constraints

         BB(1:this%ns,1:this%nc)=this%B(1:this%ns,1:this%nc)
         
         lwork=size(work)
         call dgesvd('A','A',this%ns,this%nc,BB(1:this%ns,1:this%nc),this%ns,S(1:this%ns+this%nc),U(1:this%ns,1:this%ns),this%ns,&
                     VT(1:this%nc,1:this%nc),this%nc,work(1:lwork),lwork,info)
         
         if (info.ne.0) call die('[chem_sys red_con] SVD of the constraint matrix failed')

         ! Determine rank deficiency of B
         rank_def=0
         do j=1,this%nc
            if (S(j)<S(1)*thresh) then
               rank_def=this%nc-j+1
               exit
            endif
         end do

         this%nb=this%nc-rank_def

         ! Identify all determined species
         ! Species k is determined if U(k,nb+1:ns)=0
         sp_det=0
         this%sp_order=0
         kk=0
         do k=1,this%ns
            if (norm2(U(k,this%nb+1:this%ns)).lt.thresh) then
               sp_det(k)=1	! species k is determined
               kk=kk+1
               this%sp_order(kk)=k	! determined species are first in ordering
            endif
         end do
         this%nsd=sum(sp_det)
         this%nsu=this%ns-this%nsd

         do k=1,this%ns
            if (sp_det(k).eq.0) then
               kk=kk+1
               this%sp_order(kk)=k	! undetermined species are last in ordering
            endif
         end do
         
         ! Form DEBg=[D E Bg]
         ndebg=this%nsd+this%ne+this%ng
         DEBg=0.0_WP
         DEBg(1:this%ns,this%nsd+1:ndebg)=this%B(1:this%ns,this%ncs+1:this%nc)  ! set [E Bg]
         do kk=1,this%nsd
            k=this%sp_order(kk)
            DEBg(k,kk)=1.0_WP   ! set [D]
            DEBg(k,this%nsd+1:ndebg)=0.0_WP  ! set rows of [E Bg] to zero for determined species
         end do

         ! Determine independent columns of DEBg
         call ind_col(this%ns,ndebg,this%nsd,DEBg(1:this%ns,1:ndebg),thresh,indcol(1:ndebg),info)

         if (info.lt.0) call die('[chem_sys red_con] ind_col failed')

         ! Identify all determined elements
         el_det=0
         this%el_order=0
         this%neu=0
         this%ned=0

         do k=1,this%ne
            if (indcol(this%nsd+k).eq.0) then	!	element k is determined
               this%ned=this%ned+1
               el_det(k)=1
               this%el_order(this%ned)=k	! determined elements are first in ordering
            else
               this%neu=this%neu+1			! element is undetermined
            endif
         end do

         kk=0
         EU=0.0_WP
         do k=1,this%ne
            if (el_det(k).eq.0) then
               kk=kk+1
               this%el_order(kk+this%ned)=k		! undetermined elements are last in ordering
               EU(:,kk)=this%Ein(:,k)
            endif
         end do

         ! Identify any linearly dependent columns of Bg
         if (this%ng.gt.0) then
            ngi=sum(indcol(this%nsd+this%ne+1:ndebg))
         else
            ngi=0
         endif

         ! Assemble BS (BR prior to species re-ordering)
         this%nrc=this%neu+ngi		! number of reduced constraints
         BS=0.0_WP
         kk=0
         do k=1,this%ne
            if (el_det(k).eq.0) then
               kk=kk+1
               BS(:,kk)=this%Ein(:,k)
            endif
         end do

         do k=1,this%ng
            if (indcol(this%nsd+this%ne+k).eq.1) then
               kk=kk+1
               BS(:,kk)=this%Bg(:,k)
            endif
         end do
         
         ! Assemble BR
         BR=0.0_WP
         do kk=1,this%nsu
            k=this%sp_order(this%nsd+kk)
            BR(kk,1:this%nrc)=BS(k,1:this%nrc)
         end do
         
         ! Assemble E
         this%E=0.0_WP
         do is=1,this%ns
            js=this%sp_order(is)
            do ie=1,this%ne
               je=this%el_order(ie)
               this%E(is,ie)=this%Ein(js,je)
            end do
         end do

         ! Determine modified constraint transformation matrix
         call reorder_rows(U(1:this%ns,1:this%nb),this%sp_order,PU(1:this%ns,1:this%nb))

         BTTPU=0.0_WP
         if (this%nsd.gt.0) BTTPU(1:this%nsd,1:this%nb)=PU(1:this%nsd,1:this%nb)
         if (this%nsu.gt.0) BTTPU(this%nsd+1:this%nb,1:this%nb)=matmul(transpose(BR(1:this%nsu,1:this%nrc)),PU(this%nsd+1:this%ns,1:this%nb))
         do j=1,this%nb
            BTTPU(1:this%nb,j)=BTTPU(1:this%nb,j)/S(j)
         end do
         A=0.0_WP
         A(1:this%nb,1:this%nc)=matmul(BTTPU(1:this%nb,1:this%nc),VT)

         ! Set to zero near-zero components of A
         do j=1,this%nb
            do k=1,this%nc
               if (abs(A(j,k))<thresh) A(j,k)=0.0_WP
            end do
         end do

         if (ifop<3) return

         contains
            subroutine check_input()
               use, intrinsic :: iso_fortran_env, only: output_unit
               integer :: myerr(5),i
               myerr=0
               if (this%ns.lt.1) myerr(1)=1
               if (this%ne.lt.1) myerr(2)=1
               if (this%ncs.lt.0.or.this%ncs.gt.this%ns) myerr(3)=1
               if (this%ng.lt.0) myerr(4)=1
               if (this%nc.ne. this%ne+this%ncs+this%ng) myerr(5)=1
               if (maxval(myerr).eq.0) return
               if (ifop.ge.1) then
                  write(output_unit,'(" >   chem_sys red_con: error in input")')
                  do i=1,5
                     if (myerr(i).gt.0) write(output_unit,'(" >   chem_sys red_con: error in argument,i5")') i
                  end do
               endif
            end subroutine check_input
      end subroutine red_con


end module chem_sys_class