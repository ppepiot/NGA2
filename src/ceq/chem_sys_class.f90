!> Chemical system
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

      ! Is the system initialized
      logical :: initialized                          !<=.false. !=.true. when chem_sys has been initialized

      ! Numerical parameters (when parameters are altered,chem_sys_param_set and ceq_param_def must be updated.)
      integer  :: diag                                !< Greater than 0 for diagnostics
      integer  :: lu                                  !< Logical unit number for diagnostics
      real(WP) :: T_low                               !< Lowest allowed temperature
      real(WP) :: T_high                              !< Highest allowed temperature
      real(WP) :: frac_zm                             !< Fraction of zm used in initial guess
      real(WP) :: T_tol                               !< Convergence tolerance on log(T)
      real(WP) :: ds_inc                              !< Factor by which ds is increased after success (ceq_fixed_T)
      real(WP) :: res_tol                             !< Convergence tolerance for residual (ceq_fixed_T)
      real(WP) :: ires_fac                            !< Factor by which the irreducible residual can exceed res_tol
      real(WP) :: logy_lim                            !< Upper limit on log(y)  (to prevent overflow in ceq_y)
      real(WP) :: srat_lim                            !< Lower limit on singular-value ratio (ceq_Sinv)
      real(WP) :: err_huge                            !< Upper limit on error (ceq_newt)
      real(WP) :: dec_min                             !< Minimum acceptable decrease in residual (ceq_newt)
      real(WP) :: eps_el                              !< Relative lower bound on element moles (ceq_perturb)
      real(WP) :: eps_sp                              !< Relative lower bound on species moles (ceq_perturb)
      real(WP) :: pert_tol                            !< Largest allowed normalized perturbation
      integer  :: pert_skip                           !< Set =1 to skip perturbing undetermined species

   contains
      procedure :: initialize
      procedure :: param_def=>chem_sys_param_def
      procedure :: param_set=>chem_sys_param_set
      procedure :: red_con
   end type chem_sys


   contains
      

      !> Chemical system initializer
      subroutine initialize(this,np,ns,ne,ncs,ng,P,Ein,CS,Bg,thermo_in,lu_op,diag)

         !  Specify a constrained equilibrium system.
      
         !   S.B. Pope 6/30/03

         ! Input:
         !   np    - number of phases
         !   ns    - number of species
         !   ne    - number of elements
         !   ncs   - number of constrained species
         !   ng    - number of general linear constraints
         !   Ein   - element matrix (ns x ne).  A molecule of species k contains
         !           Ein(k,j) atoms of element j.
         !   CS    - array containing indexes of constrained species (ncs x 1)
         !   Bg    - general linear constraint matrix (ns x ng)
         !   thermo_in  - thermodynamic data for species (ns x 2*ncof+1)
         !              - see below for details
         !   P     - Phase summation matrix     
         !   lu_op - logical unit for output 
         !   diag  - level of diagnostic output (0=none,1=severe errors,...5=full)
         !   (The values of lu_op and diag are stored in chem_sys,and used in calls to ceq_state.
         !    These values can be changed by,for example,call chem_sys_param_set(lu_op=3,diag=4).)
      
         ! Unconstrained equilibrium:
         !   For unconstrained equilibrium calculations (i.e.,in which the only 
         !   constraint is on the elements),set ncs=0,ng=0,CS=[] and Bg=[].
      
         ! Details of thermo_in:
         !   As in Chemkin,the thermodynamic properties of each species are given
         !   in terms of 7 non-dimensional coefficients,a(1:7).  Different values 
         !   of a(1:7) are given for two temperature ranges.  T* [K] denotes the
         !   upper limit of the lower temperature range,and the lower limit of the 
         !   upper temperature range.  The contents of the array thermo_in are: 
         !   thermo_in(k,1)   =T* for species k
         !   thermo_in(k,2:8) =a(1:7) for species k in the lower temperature range
         !   thermo_in(k,9:15)=a(1:7) for species k in the upper temperature range
         !   The relevant thermodynamic variables are:
         !   T   - temperature [K]
         !   H   - molar specific enthalpy [J/kmol]
         !   S   - molar specific entropy [J/(kmol K)]
         !   R   - universal gas constant [J/(kmol K)]
         !   Then for each species,H and S are given by:
         !   H/(RT)=a(6)/T+sum_{n=1}^{n=5} a(n) T^{n-1}/n
         !   S/R=a(1)log(T)+a(7)+sum_{n=2}^{n=5} a(n) T^{n-1}/(n-1)
      
         !  Notes:
         !    The values of lu_op and diag are stored in chem_sys,and used in subsequent calls to ceq_state.
         !    These values can be changed by,for example,call chem_sys_param_set(lu_op=3,diag=4).

         use messager,  only: die
         use mathtools, only: reorder_rows
         implicit none
         class(chem_sys), intent(inout) :: this
         integer,  intent(in) :: np,ns,ne,ncs,ng,CS(ncs),lu_op,diag
         real(WP), intent(in) :: Ein(ns,ne),Bg(ns,ng),thermo_in(ns,2*ncof+1),P(ns,Lphase:Gphase)
         real(WP) :: BR(ns,ne+ncs+ng),A(ne+ncs+ng,ne+ncs+ng)
         integer  :: nc,j,k
         real(WP) :: cpu_redcon0,cpu_redcon
      
         ! Define parameters
         call this%param_def()

         ! Check input
         if(np.ge.1) then
            this%np=np
         else
            if(diag.ge.1) then
               write(lu_op,*)'chem_sys initialize: bad input,np= ',np
               call die('np needs to be greater than zero')
            end if
         endif

         if(ns.ge.1) then
            this%ns=ns
         else
            if(diag.ge.1) then
               write(lu_op,*)'chem_sys initialize: bad input,ns= ',ns
               call die('ns needs to be greater than zero')
            end if
         endif

         if(ne.ge.1) then
            this%ne=ne
         else
            if(diag.ge.1) then
               write(lu_op,*)'chem_sys initialize: bad input,ne= ',ne
               call die('ne needs to be greater than zero')
            end if
         endif
      
         if(ncs.ge.0) then
            this%ncs=ncs
            if(ncs.ge.1) then
               allocate(this%CS(ncs))
               this%CS=CS(1:ncs)
            endif
         else
            if(diag.ge.1) then
               write(lu_op,*)'chem_sys initialize: bad input,ncs= ',ncs
               call die('ncs cannot be negative')
            end if
         endif
      
         if(ng.ge.0) then
            this%ng=ng
            if(ng.ge.1) then
               allocate(this%Bg(ns,ng))
               this%Bg=Bg(1:ns,1:ng)
            endif
         else
            if(diag.ge.1) then
               write(lu_op,*)'chem_sys initialize: bad input,ng= ',ng
               call die('ng cannot be negative') 
            end if
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
         call this%red_con(BR,A,diag,lu_op)
         
         allocate(this%BR(this%nsu,this%nrc))
         allocate(this%A(this%nb,nc))
      
         this%BR=BR(1:this%nsu,1:this%nrc)
         this%A =A(1:this%nb,1:nc)
      
         ! Re-order thermo
         call reorder_rows(thermo_in,this%sp_order,this%thermo)
      
         ! Re-order P
         call reorder_rows(P,this%sp_order,this%P)
      
         this%lu  =lu_op
         this%diag=diag
         this%initialized=.true.  ! indicate successful initialization

      end subroutine initialize


      !> Reset values of parameters in chem_sys: apart from chem_sys,all arguments are optional.
      subroutine chem_sys_param_set(this,diag,lu,T_low,T_high,frac_zm,T_tol,res_tol,ires_fac,logy_lim,srat_lim,err_huge,dec_min,eps_el,eps_sp,pert_tol,pert_skip)
         ! If this%diag=5,the values of the parameters are output.
         implicit none
         class(chem_sys), intent(inout) :: this
         integer,  intent(in), optional :: diag,lu,pert_skip
         real(WP), intent(in), optional :: T_low,T_high,frac_zm,T_tol,res_tol,ires_fac,logy_lim,srat_lim,err_huge,dec_min,eps_el,eps_sp,pert_tol
         integer :: luu
         if(present(diag))     this%diag    =diag
         if(present(lu))       this%lu      =lu
         if(present(T_low))    this%T_low   =T_low
         if(present(T_high))   this%T_high  =T_high
         if(present(frac_zm))  this%frac_zm =frac_zm
         if(present(T_tol))    this%T_tol   =T_tol
         if(present(res_tol))  this%res_tol =res_tol
         if(present(ires_fac)) this%ires_fac=ires_fac
         if(present(logy_lim)) this%logy_lim=logy_lim
         if(present(srat_lim)) this%srat_lim=srat_lim
         if(present(err_huge)) this%err_huge=err_huge
         if(present(dec_min))  this%dec_min =dec_min
         if(present(eps_el))   this%eps_el  =eps_el
         if(present(eps_sp))   this%eps_sp  =eps_sp
         if(present(pert_tol)) this%pert_tol=pert_tol
         if(present(pert_skip)) this%pert_skip=pert_skip
         if(this%diag.lt.5) return
         luu=this%lu
         write(luu,*)' '
         write(luu,*)' Parameters output from chem_sys_param_set '
         write(luu,*)' '
         write(luu,'(a,i4)') 'diag     =',this%diag
         write(luu,'(a,i4)') 'lu       =',this%lu  
         write(luu,'(a,i4)') 'pert_skip=',this%pert_skip 
         write(luu,'(a,1p,e13.4)') 'T_low   =',this%T_low
         write(luu,'(a,1p,e13.4)') 'T_high  =',this%T_high
         write(luu,'(a,1p,e13.4)') 'frac_zm =',this%frac_zm
         write(luu,'(a,1p,e13.4)') 'T_tol   =',this%T_tol
         write(luu,'(a,1p,e13.4)') 'res_tol =',this%res_tol
         write(luu,'(a,1p,e13.4)') 'ires_fac=',this%ires_fac
         write(luu,'(a,1p,e13.4)') 'logy_lim=',this%logy_lim
         write(luu,'(a,1p,e13.4)') 'srat_lim=',this%srat_lim
         write(luu,'(a,1p,e13.4)') 'err_huge=',this%err_huge
         write(luu,'(a,1p,e13.4)') 'dec_min =',this%dec_min
         write(luu,'(a,1p,e13.4)') 'eps_el  =',this%eps_el
         write(luu,'(a,1p,e13.4)') 'eps_sp  =',this%eps_sp
         write(luu,'(a,1p,e13.4)') 'pert_tol=',this%pert_tol
      end subroutine chem_sys_param_set


      !> Reset parameters in chem_sys to their default settings.
      subroutine chem_sys_param_def(this)
         implicit none
         class(chem_sys), intent(inout) :: this
         this%initialized=.false. !=.true. when chem_sys has been initialized
         this%diag    =1          ! >0 for diagnostics
         this%lu      =0          ! logical unit number for diagnostics
         this%T_low   =250.0_WP   ! lowest allowed temperature
         this%T_high  =5000.0_WP  ! highest allowed temperature
         this%frac_zm =1e-1       ! fraction of zm used in initial guess
         this%T_tol   =1e-6       ! convergence tolerance on log(T)
         this%res_tol =1e-9       ! convergence tolerance for residual (ceq_fixed_T)
         this%ires_fac=1e2        ! factor by which the irreducible residual can exceed res_tol
         this%logy_lim=120.0_WP   ! upper limit on log(y)  (to prevent overflow in ceq_y)
         this%srat_lim=1e-9       ! lower limit on singular-value ratio (ceq_Sinv)
         this%err_huge=1e6        ! upper limit on error (ceq_newt)
         this%dec_min =0.5_WP     ! minimum acceptable decrease in residual (ceq_newt)
         this%eps_el  =1e-9       ! relative lower bound on element moles (ceq_perturb)
         this%eps_sp  =1e-9       ! relative lower bound on species moles (ceq_perturb)
         this%pert_tol=1e-4       ! largest allowed normalized perturbation
         this%pert_skip= 0        ! set =1 to skip perturbing undetermined species
      end subroutine chem_sys_param_def


      !> Set the reduced constraint matrix BR,and determine the ordering of the elements and species
      subroutine red_con(this,BR,A,ifop,lud)
         use messager,  only: die
         use mathtools, only: ind_col,reorder_rows
         implicit none
         class(chem_sys), intent(inout) :: this
         integer,  intent(in)  :: ifop,lud
         real(WP), intent(out) :: BR(this%ns,this%nc),A(this%nc,this%nc)
         integer  :: k,info,lwork,rank_def,sp_det(this%ns),kk,ndebg,el_det(this%ne),ngi,is,js,ie,je,j,indcol(this%ns+this%ne+this%ng)
         real(WP) :: Sc(this%ns,this%ncs),S(this%ns+this%nc),U(this%ns,this%ns),VT(this%nc,this%nc),work(4*(this%ns+this%nc)*(this%ns+this%nc))
         real(WP) :: DEBg(this%ns,this%ns+this%ne+this%ng),EU(this%ns,this%ne),BS(this%ns,this%nc),PU(this%ns,this%nc),BTTPU(this%nc,this%nc),BB(this%ns,this%nc)
         real(WP), parameter :: thresh=1e-10 ! threshold for singular values and vectors

         ! print out dimensions and input
         if(ifop.ge.5) then
            write(lud,*)' '
            write(lud,*)'Input to red_con'
            write(lud,*)'ne=',this%ne
            write(lud,*)'ns=',this%ns
            write(lud,*)'ncs= ',this%ncs
            write(lud,*)'ng=  ',this%ng
            write(lud,*)'nc=',this%nc
            write(lud,*)'Constrained species'
            write(lud,'((20i4))')this%CS
         endif
         
         ! Check input
         call check_input()

         if(.true.) then !  new test for constrained species being distinct SBP 2/6/09 

            do j=1,this%ncs-1
               do k=j+1,this%ncs
                  if(this%CS(j).eq.this%CS(k)) then
                     write(0,*)'red_con: CS not distinct'
                     call die('red_con: CS not distinct')
                  endif
               end do
            end do

            else  ! original test which inappropriately assumes ordeering

               if(this%ncs.gt.1) then
                  if(minval(this%CS(2:this%ncs)-this%CS(1:this%ncs-1)).lt.1) then
                     write(0,*)'red_con: CS not distinct'
                     call die('red_con: CS not distinct')
                  endif
               endif

         endif  !  end of distinct test

         this%B=0.0_WP     ! set basic constraint matrix
         ! set constrained species matrix
         if(this%ncs.gt.0) then
            Sc=0.0_WP
            do k=1,this%ncs 
               Sc(this%CS(k),k)=1.0_WP
            end do 
            this%B(:,1:this%ncs)=Sc              ! first ncs columns - constrained species
         endif
         
         this%B(1:this%ns,this%ncs+1:this%ncs+this%ne)=this%Ein(1:this%ns,1:this%ne)      ! next ne columns - element vectors

         if(this%ng.gt.0) &
         this%B(1:this%ns,this%ne+this%ncs+1:this%nc)=this%Bg(1:this%ns,1:this%ng)        ! last ng columns - general constraints

         BB(1:this%ns,1:this%nc)=this%B(1:this%ns,1:this%nc)
         
         lwork=size(work)
         call dgesvd('A','A',this%ns,this%nc,BB(1:this%ns,1:this%nc),this%ns,S(1:this%ns+this%nc),U(1:this%ns,1:this%ns),this%ns,&
                     VT(1:this%nc,1:this%nc),this%nc,work(1:lwork),lwork,info)
         
         if(info.ne.0) then
            if(ifop.ge.1) then
               write(lud,*)'red_con: svd failed,info= ',info
               call die('red_con: svd failed')
            endif
         end if

         rank_def=0     !  determine rank deficiency of B
         do j=1,this%nc
            if(S(j)<S(1)*thresh) then
               rank_def=this%nc-j+1
               exit
            endif
         end do

         if(rank_def.gt.0.and.ifop.ge.3) then
            write(lud,*)' '
            write(lud,*)'red_con: basic constraint matrix is singular'
            write(lud,*)'rank deficiency= ',rank_def
         endif

         this%nb=this%nc-rank_def

         ! Sdentify all determined species
         ! Species k is determined if U(k,nb+1:ns)=0
         sp_det=0
         this%sp_order=0
         kk=0
         do k=1,this%ns
            if(norm2(U(k,this%nb+1:this%ns)).lt.thresh) then
               sp_det(k)=1	! species k is determined
               kk=kk+1
               this%sp_order(kk)=k	! determined species are first in ordering
            endif
         end do
         this%nsd=sum(sp_det)
         this%nsu=this%ns-this%nsd

         do k=1,this%ns
            if(sp_det(k).eq.0) then
               kk=kk+1
               this%sp_order(kk)=k	! undetermined species are last in ordering
            endif
         end do

         if(ifop.ge.3) then
            write(lud,*)' '
            write(lud,'(a,i5)')'nsd=',this%nsd
            write(lud,'(a,i5)')'nsu=',this%nsu
            if(ifop.ge.4) then
               write(lud,*)'sp_order='
               write(lud,'((20i4))')this%sp_order
            endif
         endif
         
         !  form DEBg=[D E Bg]
         ndebg=this%nsd+this%ne+this%ng
         DEBg=0.0_WP
         DEBg(1:this%ns,this%nsd+1:ndebg)=this%B(1:this%ns,this%ncs+1:this%nc)  ! set [E Bg]
         do kk=1,this%nsd
            k=this%sp_order(kk)
            DEBg(k,kk)=1.0_WP   ! set [D]
            DEBg(k,this%nsd+1:ndebg)=0.0_WP  ! set rows of [E Bg] to zero for determined species
         end do

         ! determine independent columns of DEBg
         call ind_col(this%ns,ndebg,this%nsd,DEBg(1:this%ns,1:ndebg),thresh,indcol(1:ndebg),info)

         if(info.lt.0) then
            if(ifop.ge.1) then
               write(lud,*) 'red_con: ind_col failed; info=',info
               call die('red_con: ind_col failed')
            end if
         endif

         !  identify all determined elements
         el_det=0
         this%el_order=0
         this%neu=0
         this%ned=0

         do k=1,this%ne
            if(indcol(this%nsd+k).eq.0) then	!	element k is determined
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
            if(el_det(k).eq.0) then
               kk=kk+1
               this%el_order(kk+this%ned)=k		! undetermined elements are last in ordering
               EU(:,kk)=this%Ein(:,k)
            endif
         end do

         if(ifop.ge.3) then
            write(lud,*)' '
            write(lud,*)'ned=',this%ned
            write(lud,*)'neu=',this%neu
            if(ifop.ge.4) then
               write(lud,*)'el_order='
               write(lud,'((20i4))')this%el_order
            endif
         endif

         !  identify any linearly dependent columns of Bg
         if(this%ng.gt.0) then
            ngi=sum(indcol(this%nsd+this%ne+1:ndebg))
         else
            ngi=0
         endif

         !  assemble BS (BR prior to species re-ordering)
         this%nrc=this%neu+ngi		! number of reduced constraints
         BS=0.0_WP
         kk=0
         do k=1,this%ne
            if(el_det(k).eq.0) then
               kk=kk+1
               BS(:,kk)=this%Ein(:,k)
            endif
         end do

         do k=1,this%ng
            if(indcol(this%nsd+this%ne+k).eq.1) then
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
         if(this%nsd.gt.0) BTTPU(1:this%nsd,1:this%nb)=PU(1:this%nsd,1:this%nb)
         if(this%nsu.gt.0) BTTPU(this%nsd+1:this%nb,1:this%nb)=matmul(transpose(BR(1:this%nsu,1:this%nrc)),PU(this%nsd+1:this%ns,1:this%nb))
         do j=1,this%nb
            BTTPU(1:this%nb,j)=BTTPU(1:this%nb,j)/S(j)
         end do
         A=0.0_WP
         A(1:this%nb,1:this%nc)=matmul(BTTPU(1:this%nb,1:this%nc),VT)

         ! Set to zero near-zero components of A
         do j=1,this%nb
            do k=1,this%nc
               if(abs(A(j,k))<thresh) A(j,k)=0.0_WP
            end do
         end do

         if(ifop<3) return

         write(lud,*)'number of elements,             ne= ' ,this%ne
         write(lud,*)'number of species ,             ns= ' ,this%ns
         write(lud,*)'number of constrained species ,ncs= ' ,this%ncs
         write(lud,*)'number of general constraints ,ng=  ' ,this%ng
         write(lud,*)'number of determined elements ,ned= ' ,this%ned
         write(lud,*)'number of determined species , nsd= ' ,this%nsd
         write(lud,*)'  '
         write(lud,*)'number of undetermined species , nsu= ' ,this%nsu
         write(lud,*)'number of undetermined elements ,neu= ' ,this%neu
         write(lud,*)'number of indep. gen. constr.,   ngi =' ,ngi
         write(lud,*)'number of reduced constraints,   nrc= ' ,this%nrc
         write(lud,*)' '

      contains
         subroutine check_input()
            integer :: myerr(5),i
            myerr=0
            if(this%ns.lt.1) myerr(1)=1
            if(this%ne.lt.1) myerr(2)=1
            if(this%ncs.lt.0.or.this%ncs.gt.this%ns) myerr(3)=1
            if(this%ng.lt.0) myerr(4)=1
            if(this%nc.ne. this%ne+this%ncs+this%ng) myerr(5)=1
            if(maxval(myerr).eq.0) return
            if(ifop.ge.1) then
               write(lud,*)'red_con: error in input'
               do i=1,5
                  if(myerr(i).gt.0) write(lud,'(a,i3)') '   error in argument ',i
                  call die('')
               end do
            endif
         end subroutine check_input
      end subroutine red_con


end module chem_sys_class