!> Example to read in YAML thermodynamic file
module simulation
   use precision,        only: WP
   use string,           only: str_short,str_medium
   use monitor_class,    only: monitor
   use YAMLRead,         only: YAMLElement
   use chem_sys_class,   only: chem_sys,Lphase,Gphase
   use chem_state_class, only: chem_state
   implicit none
   private
   
   !> The array of the species. Eeach stored as a YAMLElement object
   type(YAMLElement), dimension(:), allocatable :: species

   !> Number of species, elements, constraints, and phases
   integer :: ns,ne,nsd,nsu,ned,neu,nc,nrc,np=2,ncs

   !> Species and elements names
   character(len=str_medium), dimension(:), allocatable :: sp_names
   character(len=str_short),  dimension(:), allocatable :: e_names

   !> Elements matrix
   real(WP), dimension(:,:), allocatable :: elem_mat

   !> Phase summation matrix
   real(WP), dimension(:,:), allocatable :: phse_mat

   !> Thermodynamic quantities
   real(WP) :: T
   real(WP) :: PoPref

   !> Chemical system and state
   type(chem_sys)   :: sys
   type(chem_state) :: state

   !> Newton solver
   integer  :: iter_max
   real(WP) :: tol
   integer  :: iter
   real(WP) :: Rnorm,dxnorm

   !> Private work arrays
   real(WP), dimension(:),   allocatable :: N,Nbar,Nu,Nd,Neq ! Mole numbers
   real(WP), dimension(:),   allocatable :: R,Rd             ! Residuals
   real(WP), dimension(:),   allocatable :: x,x0,dx          ! Chemical state unknowns
   real(WP), dimension(:),   allocatable :: gort             ! Normalized Gibbs free energy (Molar G over RT)
   real(WP), dimension(:,:), allocatable :: BtildeT,PtildeT  ! Coefficient matrices

   !> Simulation monitor file
   type(monitor) :: nfile

   !> Simulation subroutines
   public :: simulation_init,simulation_run,simulation_final
   

contains


   !> Get the square root of moles
   function get_y(xin,g) result(y)
      real(WP), dimension(nrc+np), intent(in) :: xin
      real(WP), dimension(nsu),    intent(in) :: g
      real(WP), dimension(nsu) :: y
      y=exp(0.5_WP*(-g+matmul(sys%BR,xin(1:nrc))+matmul(sys%P(nsd+1:ns,:),xin(nrc+1:nrc+np))))
   end function get_y


   !> Get the residual vector
   function get_res(y) result(res)
      real(WP), dimension(nsu), intent(in) :: y
      real(WP), dimension(nrc+np) :: res
      res(1:nrc)=matmul(BtildeT,y)-state%cr
      res(nrc+1:nrc+np)=matmul(PtildeT,y)-Nbar+Rd
   end function get_res


   !> Initialization of problem solver
   subroutine simulation_init
      use param,     only: param_read,param_getsize
      use messager,  only: die
      use mathtools, only: lss
      implicit none
      real(WP), allocatable :: nasa_coef(:,:)
      character(len=str_medium), dimension(:), allocatable :: const_sp
      integer,  dimension(:), allocatable :: CS
      real(WP), dimension(:), allocatable :: N_init

      ! Parse the mechanism file
      parse_mech: block
         use YAMLRead, only: YAMLHandler,YAMLSequence,YAMLMap,yaml_open_file,yaml_start_from_sequence,yaml_close_file
         character(len=str_medium) :: mch_file
         character(len=str_short), dimension(:), allocatable :: sp_names_copy,const_sp_copy
         real(WP), dimension(:), allocatable :: N_init_copy
         type(YAMLHandler)  :: domain
         type(YAMLSequence) :: sp_list,phases,elements
         type(YAMLElement)  :: sp,gas
         type(YAMLMap)      :: thermo,comp
         integer :: isc,nn,i,j,k,e,code
         character(len=:), allocatable :: name_arr(:)
         character(len=:), allocatable :: name
         real(WP), allocatable :: T_range(:)
         real(WP), dimension(:,:), allocatable :: a
         logical :: new_elem
         ! Get the target species from input
         ns=param_getsize('Species')
         nn=param_getsize('Initial moles')
         ncs=param_getsize('Constrained species')
         if (ns.ne.nn) call die('Unequal number of species and moles in the input file.')
         allocate(sp_names(1:ns))
         allocate(sp_names_copy(1:ns))
         allocate(N_init(1:ns))
         allocate(N_init_copy(1:ns))
         allocate(const_sp(1:ncs))
         allocate(const_sp_copy(1:ncs))
         allocate(CS(ncs));    CS=[1,4]
         call param_read('Species',sp_names)
         call param_read('Initial moles',N_init)
         call param_read('Constrained species',const_sp)
         sp_names_copy=sp_names
         N_init_copy=N_init
         const_sp_copy=const_sp
         ! Read the mechanism file path
         call param_read('Mechanism file',mch_file)
         ! Open the mechanism
         domain=yaml_open_file(trim(mch_file))
         ! Get the list of all species
         sp_list=yaml_start_from_sequence(domain,'species')
         ! Extract the target species from the mechanism
         allocate(species(1:ns))
         nn=0
         k=0
         do isc=0,sp_list%size-1 ! Index in YAMLSequence starts from 0
            sp=sp_list%element(isc)
            name_arr=sp%value_str('name',code)
            name=''
            do i=1,size(name_arr)
               name=trim(name//name_arr(i))
            end do
            do i=1,ns
               if (sp_names_copy(i).eq.name) then
                  nn=nn+1
                  species(nn)=sp
                  sp_names(nn)=name
                  N_init(nn)=N_init_copy(i)
                  do j=1,ncs
                     if (const_sp_copy(j).eq.name) then
                        k=k+1
                        const_sp(k)=name
                        CS(k)=nn
                     end if
                  end do
               end if
            end do
            call sp%destroy()
         end do
         if(nn.ne.ns) call die('Some species are missing in the mechanism file.')
         ! Get the elements that exist in the target species
         phases=yaml_start_from_sequence(domain,'phases')
         gas=phases%element(0)
         elements=gas%value_sequence('elements',code) ! For some reason I couldn't directly get the element names using yaml-fortran
         allocate(e_names(elements%size)); e_names=''
         ne=0
         do isc=1,ns
            sp=species(isc)
            comp=sp%value_map('composition')
            do e=1,size(comp%labels)
               name=''
               do i=1,size(comp%labels(e)%str)
                  name=name//trim(comp%labels(e)%str(i))
               end do
               new_elem=.true.
               do i=1,ne
                  if (trim(e_names(i)).eq.trim(name)) then
                     new_elem=.false.
                     exit
                  end if
               end do
               if (new_elem) then
                  ne=ne+1
                  e_names(ne)=trim(name)
               end if
            end do
         end do
         e_names=e_names(1:ne)
         ! Form the element matrix
         allocate(elem_mat(ns,ne)); elem_mat=0.0_WP
         do isc=1,ns
            sp=species(isc)
            comp=sp%value_map('composition')
            do i=1,size(comp%labels)
               name=''
               do nn=1,size(comp%labels(i)%str)
                  name=name//trim(comp%labels(i)%str(nn))
               end do
               do e=1,ne
                  if (trim(e_names(e)).eq.trim(name)) elem_mat(isc,e)=real(comp%value_int(name,code),WP)
               end do
            end do
         end do
         ! Read the NASA-7 polynomials
         allocate(nasa_coef(1:ns,15)); nasa_coef=0.0_WP
         allocate(a(1:2,1:7)); a=0.0_WP
         do isc=1,ns
            sp=species(isc)
            thermo=sp%value_map('thermo')
            T_range=thermo%value_double_1d('temperature-ranges',code)
            select case (size(T_range))
            case (3)
               a=thermo%value_double_2d('data',code)
            case (2)
               a(1,:)=thermo%value_double_1d('data',code)
            case default
               call die('Invalid temperature range')
            end select
            nasa_coef(isc,1   )=T_range(2)
            nasa_coef(isc,2:8 )=a(1,:)
            nasa_coef(isc,9:15)=a(2,:)
         end do
         ! Form the phase summation matrix
         allocate(phse_mat(ns,Lphase:Gphase)); phse_mat(:,Lphase)=0.0_WP; phse_mat(:,Gphase)=1.0_WP
         do isc=1,ns
            if (len_trim(sp_names(isc)).ge.3) then
               if (sp_names(isc)(len_trim(sp_names(isc))-2:len_trim(sp_names(isc))).eq.'(L)') then
                  phse_mat(isc,Lphase)=1.0_WP
                  phse_mat(isc,Gphase)=0.0_WP
               end if
            end if
         end do
         ! Close the mechanism file and clean up
         call yaml_close_file(domain)
         call sp_list%destroy()
         call sp%destroy()
         call comp%destroy()
         call thermo%destroy()
         deallocate(sp_names_copy,N_init_copy,const_sp_copy)
      end block parse_mech

      ! Initialize the chemical system
      ceq_init: block
         integer :: ng=1
         integer :: lu,iostat,iret,info
         real(WP), dimension(:,:), allocatable :: Bg
         real(WP), dimension(:),   allocatable :: stats
         real(WP) :: HoR,T_eq
         integer :: isc
         ! Allocate arrays
         allocate(Bg(ns,ng));  Bg=0.0_WP
         allocate(N(ns))
         allocate(stats(20))
         ! Create the general constraints
         do isc=1,ns
            if (sp_names(isc).eq.'H2O')    Bg(isc,1)=1.0_WP
            if (sp_names(isc).eq.'H2O(L)') Bg(isc,1)=1.0_WP
         end do
         ! Read inputs
         call param_read('Temperature',T)
         call param_read('Pressure over reference pressure',PoPref)
         ! Print the initial conditions
         print*,'Initial moles:'
         do isc=1,ns
            print*,trim(sp_names(isc)),': ',N_init(isc)
         end do
         ! Open CEQ file
         open(unit=lu,file='ceq_out',status='replace',action='write',iostat=iostat)
         ! Inizialize the chemical system
         call sys%initialize(np=np,ns=ns,ne=ne,ncs=ncs,ng=ng,P=phse_mat,Ein=elem_mat,CS=CS,Bg=Bg,thermo_in=nasa_coef,lu_op=lu,diag=5)
         if (iret.lt.0) call die('chem_sys initialization failed.')
         nsd=sys%nsd
         nsu=sys%nsu
         ned=sys%ned
         neu=sys%neu
         nc =sys%nc
         nrc=sys%nrc
         ! Initialize the chemical state
         call state%initialize(sys=sys,N=N_init,p_Pa=101325.0_WP,T=T,N_eq=N,T_eq=T_eq,HoR_eq=HoR,stats=stats,info=info)
         ! Close CEQ file
         close(lu)
         ! Error check
         if (info.lt.0) call die('chem_state failed.')
         ! CEQ initialization of moles
         print*,'CEQ initialization of moles:'
         do isc=1,ns
            print*,trim(sp_names(isc)),': ',N(isc)
         end do
         print*,'Equilibrium temperature = ',T_eq, '(k)'
         ! Deallocate arrays
         deallocate(Bg,stats)
      end block ceq_init

      ! Initialize the chemical state solution vector
      x_init: block
         integer :: info
         real(WP), dimension(:), allocatable :: rhs,lam
         ! Allocate arrays
         allocate(x   (nrc+np)); x=0.0_WP
         allocate(x0  (nrc+np)); x0=0.0_WP
         allocate(dx  (nrc+np)); dx=0.0_WP
         allocate(Nbar(np));     Nbar=0.0_WP
         allocate(Nu  (nsu));    Nu=0.0_WP
         allocate(Nd  (nsd));    Nd=0.0_WP
         allocate(Neq (ns));     Neq=0.0_WP
         allocate(R   (nrc+np)); R=0.0_WP
         allocate(Rd  (np));     Rd=0.0_WP
         allocate(gort(nsu));    gort=0.0_WP
         allocate(rhs (nrc))
         allocate(lam (nrc))
         ! Get the species and phase moles
         Nd=state%Nd
         Nu=state%Nu
         Nbar=matmul(transpose(sys%P),[Nd,Nu])
         ! Calculate the contribution of Nd in the residual
         Rd=matmul(transpose(sys%P(1:nsd,:)),Nd)
         ! Calculate the Lagrange multipliers
         call state%get_gort(nsu,T,PoPref,sys%thermo(nsd+1:ns,:),sys%P(nsd+1:ns,Gphase),gort)
         rhs=log(Nu)-matmul(sys%P(nsd+1:ns,:),log(Nbar))+gort
         call lss(nsu,nrc,sys%BR,rhs,lam,info)
         if (info.ne.0) call die('Least squares for lambda initialization failed.')
         ! Set the initial solution vector
         x(1:nrc)=lam
         x(nrc+1:nrc+np)=log(Nbar)
         ! Deallocate arrays
         deallocate(rhs,lam)
      end block x_init

      ! Read in newton iterations inputs
      call param_read('Tolerance',tol)
      call param_read('Max number of iterations',iter_max)

      ! Create a monitor file
      create_monitor: block
         nfile=monitor(.true.,'newton')
         call nfile%add_column(iter,'Iter number')
         call nfile%add_column(Rnorm,'R norm')
         call nfile%add_column(dxnorm,'dx norm')
         ! call nfile%write()
      end block create_monitor

      ! Deallocate arrays
      deallocate(nasa_coef,const_sp,CS,N_init)

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use messager, only: die
      implicit none
      integer :: i,j,iJ,jJ,isc,info
      real(WP), dimension(:,:), allocatable :: Jac,BTB,PTP,BTP
      real(WP), dimension(:,:), allocatable :: Btilde,Ptilde
      real(WP), dimension(:),   allocatable :: y
      real(WP), dimension(:),   allocatable :: S,work
      real(WP) :: rcond
      integer  :: rank,lwork

      ! Allocate arrays
      allocate(Jac(nrc+np,nrc+np)); Jac=0.0_WP
      allocate(BTB(nrc,nrc));       BTB=0.0_WP
      allocate(BTP(nrc,np));        BTP=0.0_WP
      allocate(PTP(np,np));         PTP=0.0_WP
      allocate(Btilde(nsu,nrc));    Btilde=0.0_WP
      allocate(Ptilde(nsu,np));     Ptilde=0.0_WP
      allocate(BtildeT(nrc,nsu));   BtildeT=0.0_WP
      allocate(PtildeT(np,nsu));    PtildeT=0.0_WP
      allocate(y(nsu));             y=0.0_WP
      allocate(S(nrc+np))
      lwork=10*(nrc+np)
      allocate(work(lwork))
      rcond=-1.0_WP

      ! Newton-Raphson
      iter=0
      y=get_y(x,gort)
      Rnorm=10*tol
      do while(Rnorm.ge.tol)
         ! Increment iteration number
         iter=iter+1
         if (iter.gt.iter_max) then
            iter=iter-1
            exit
         end if
         ! Remember the old solution
         x0=x
         ! Build the Jacobian matrix
         do j=1,nrc
            Btilde(:,j)=y*sys%BR(:,j)
         end do
         do j=1,np
            Ptilde(:,j)=y*sys%P(nsd+1:ns,j)
         end do
         BtildeT=transpose(Btilde)
         PtildeT=transpose(Ptilde)
         BTB=matmul(BtildeT,Btilde)
         PTP=matmul(PtildeT,Ptilde)
         BTP=matmul(BtildeT,Ptilde)
         do j=1,nrc
            jJ=j
            do i=1,j
               iJ=i
               Jac(iJ,jJ)=BTB(i,j)
            end do
         end do
         do j=1,np
            jJ=nrc+j
            do i=1,nrc
               iJ=i
               Jac(iJ,jJ)=BTP(i,j)
            end do
         end do
         do j=1,np
            jJ=nrc+j
            do i=1,j
               iJ=nrc+i
               Jac(iJ,jJ)=PTP(i,j)
            end do
         end do
         do j=1,nrc+np-1
            do i=j+1,nrc+np
               Jac(i,j)=Jac(j,i)
            end do
         end do
         do i=1,np
            iJ=nrc+i
            jJ=nrc+i
            Jac(iJ,jJ)=Jac(iJ,jJ)-Nbar(i)
         end do
         ! print*,'J = '
         ! do i=1,np+nrc
         !    print*,Jac(i,:)
         ! end do
         ! call die('')
         ! Get the residual error
         R=get_res(y)
         Rnorm=norm2(R)
         ! Solve for dx
         dx=-R
         call dgelss(nrc+np,nrc+np,1,Jac,nrc+np,dx,nrc+np,S,rcond,rank,work,lwork,info)
         if (rank.ne.nrc+np) call die('Jacobian is not full rank')
         if (info.ne.0) call die('Least-squares solver failed')
         dxnorm=norm2(dx)
         ! Update the solution
         x=x0+dx
         ! Get the species and phase moles
         y=get_y(x,gort)
         Nu=y*y
         Nbar=exp(x(nrc+1:nrc+np))
         ! Perform and output monitoring
         call nfile%write()
      end do

      ! Output
      Neq=[Nd,Nu]
      do isc=1,ns
         N(sys%sp_order(isc))=Neq(isc)
      end do
      print*,'Number of iterations = ', iter
      print*,'Residal error = ', R
      print*,'Equilibrium moles:'
      do isc=1,ns
         print*,trim(sp_names(isc)),': ',N(isc)
      end do

      ! Deallocate arrays
      deallocate(Jac,BTB,BTP,PTP,Btilde,Ptilde,y,S,work)

   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects-need destructors
      deallocate(species,sp_names,e_names,elem_mat,phse_mat)
      deallocate(N,Nbar,Nu,Nd,Neq,R,Rd,x,x0,dx,gort,BtildeT,PtildeT)

   end subroutine simulation_final

   
end module simulation
