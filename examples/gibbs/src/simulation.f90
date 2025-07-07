!> Example to read in YAML thermodynamic file
module simulation
   use precision, only: WP
   use YAMLRead,  only: YAMLElement
   use string,    only: str_short,str_medium
   use ceq_types, only: sys_type
   implicit none
   private
   
   ! The array of the species. Eeach stored as a YAMLElement object
   type(YAMLElement), dimension(:), allocatable :: species

   ! Number of species, elements, constraints, and phases
   integer :: ns,ne,nc,np=2

   ! Species and elements names
   character(len=str_medium), dimension(:), allocatable :: sp_names
   character(len=str_short),  dimension(:), allocatable :: e_names

   ! Elements matrix
   real(WP), dimension(:,:), allocatable :: elem_mat

   ! Phase summation matrix
   real(WP), dimension(:,:), allocatable :: phse_mat

   ! Thermodynamic quantities
   real(WP) :: T
   real(WP) :: PoPref=1.0_WP

   ! Phase indices (need to get these from two_phase classes)
   integer :: Lphase=0,Gphase=1

   ! Chemical system
   type(sys_type), pointer :: sys

   ! Working arrays
   real(WP), dimension(:),   allocatable :: N,Nbar          ! Mole numbers
   real(WP), dimension(:),   allocatable :: R               ! Residuals
   real(WP), dimension(:),   allocatable :: constraints     ! Residuals
   real(WP), dimension(:),   allocatable :: x,x0,dx         ! Chemical state unknowns
   real(WP), dimension(:),   allocatable :: gort            ! Normalized Gibbs free energy (Molar G over RT)
   real(WP), dimension(:,:), allocatable :: BtildeT,PtildeT

   ! Simulation subroutines
   public :: simulation_init,simulation_run,simulation_final
   

contains


   !> Get the residual vector
   function get_y(xin,g) result(y)
      real(WP), dimension(nc+np), intent(in) :: xin
      real(WP), dimension(ns),    intent(in) :: g
      real(WP), dimension(ns) :: y
      y=exp(0.5_WP*(-g+matmul(sys%B,xin(1:nc))+matmul(sys%P,xin(nc+1:nc+np))))
   end function get_y


   !> Get the residual vector
   function get_res(y) result(res)
      real(WP), dimension(ns), intent(in) :: y
      real(WP), dimension(nc+np) :: res
      res(1:nc)=matmul(BtildeT,y)-constraints
      res(1:np)=matmul(PtildeT,y)-Nbar
   end function get_res


   !> Initialization of problem solver
   subroutine simulation_init
      use param,    only: param_read,param_getsize
      use messager, only: die
      implicit none
      real(WP), allocatable :: nasa_coef(:,:)

      ! Parse the mechanism file
      parse_mech: block
         use YAMLRead, only: YAMLHandler,YAMLSequence,YAMLMap,yaml_open_file,yaml_start_from_sequence,yaml_close_file
         character(len=str_medium) :: yaml_file
         character(len=str_short), dimension(:), allocatable :: sp_names_copy
         type(YAMLHandler)  :: domain
         type(YAMLSequence) :: sp_list,phases,elements
         type(YAMLElement)  :: sp,gas
         type(YAMLMap)      :: thermo,comp
         integer :: isc,nn,i,e,code
         character(len=:), allocatable :: name_arr(:)
         character(len=:), allocatable :: name
         real(WP), allocatable :: T_range(:)
         real(WP), dimension(:,:), allocatable :: a
         logical :: new_elem
         ! Get the target species from input
         ns=param_getsize('Species')
         allocate(sp_names(1:ns))
         allocate(sp_names_copy(1:ns))
         call param_read('Species',sp_names)
         sp_names_copy=sp_names
         ! Read the mechanism file path
         call param_read('YAML file',yaml_file)
         ! Open the mechanism
         domain=yaml_open_file(trim(yaml_file))
         ! Get the list of all species
         sp_list=yaml_start_from_sequence(domain,'species')
         ! Extract the target species from the mechanism
         allocate(species(1:ns))
         nn=0
         do isc=0,sp_list%size-1 ! Index in YAMLSequence starts from 0
            sp=sp_list%element(isc)
            name_arr=sp%value_str('name',code)
            name=''
            do i=1,size(name_arr)
               name=trim(name//name_arr(i))
            end do
            if (any(sp_names_copy.eq.name)) then
               nn=nn+1
               species(nn)=sp
               sp_names(nn)=name
            end if
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
      end block parse_mech

      ! Initialize the chemical system
      sys_init: block
         use ceq_system,  only: ceq_sys_init
         use ceq_state_m, only: ceq_state
         integer :: ncs=0,ng=0
         integer :: lu,iostat,iret,info
         real(WP), dimension(:,:), allocatable :: Bg
         integer,  dimension(:),   allocatable :: CS
         real(WP), dimension(:),   allocatable :: N_init,c,stats
         real(WP) :: HoR
         integer :: isc
         ! Allocate arrays
         allocate(CS(ncs))
         allocate(Bg(ns,ng))
         allocate(N(ns))
         allocate(N_init(ns)); N_init=[0.0_WP,1.0_WP]
         allocate(c(ng))
         allocate(stats(20))
         ! Print the initial conditions
         print*,'Initial moles:'
         do isc=1,ns
            print*,trim(sp_names(isc)),': ',N_init(isc)
         end do
         ! Open CEQ file
         open(unit=lu,file="ceq_out",status="replace",action="write",iostat=iostat)
         ! Inizialize the system
         call ceq_sys_init(ns=ns,ne=ne,ncs=ncs,ng=ng,Ein=elem_mat,CS=CS,Bg=Bg,thermo_in=nasa_coef,P=phse_mat,lu_op=lu,diag=5,sys=sys,iret=iret)
         if (iret.lt.0) call die('System initialization failed.')
         ! Get the equilibrium state (I commented out the equilibrium calculcations in ceq_state and made it output the initial mole numberes found by maxmin and ming. See parts with comment "DEBUG")
         call ceq_state(sys=sys,N=N_init,p_Pa=101325.0_WP,T=373.15_WP,N_eq=N,T_eq=T,HoR_eq=HoR,stats=stats,info=info)
         ! call ceq_state(sys=sys,N=N_init,p_Pa=101325.0_WP,N_h=N_init,T_h=380.0_WP,N_eq=N,T_eq=T,HoR_eq=HoR,stats=stats,info=info)
         ! Close CEQ file
         close(lu)
         ! Error check
         if (info.lt.0) call die('ceq_state failed.')
         ! CEQ initialization of moles
         print*,'CEQ initialization of moles:'
         do isc=1,ns
            print*,trim(sp_names(isc)),': ',N(isc)
         end do
         print*,'Equilibrium temperature = ',T, '(k)'
         nc=sys%nc
         allocate(constraints(nc))
         constraints=matmul(N_init,sys%B)
         N=[0.5_WP,0.5_WP]
      end block sys_init

      ! Initialize the chemical state solution vector
      x_init: block
         integer :: info
         real(WP), dimension(:),   allocatable :: rhs,lam
         ! Allocate arrays
         allocate(x (nc+np)); x=0.0_WP
         allocate(x0(nc+np)); x0=0.0_WP
         allocate(dx(nc+np)); dx=0.0_WP
         allocate(Nbar(np));  Nbar=0.0_WP
         allocate(R (nc+np)); R=0.0_WP
         allocate(gort(ns));  gort=0.0_WP
         allocate(rhs(nc))
         allocate(lam(nc))
         ! Check phase matrix
         if (maxval(abs(sys%P-phse_mat)).gt.0.0_WP) call die('P has changed after reordering')
         ! Get the phase moles
         Nbar=matmul(transpose(sys%P),N)
         ! Calculate the Lagrange multipliers
         call ceq_g(ns,T,PoPref,sys%thermo,sys%P(:,Gphase),gort)
         rhs=log(N)-matmul(sys%P,log(Nbar))+gort
         call ceq_lss(ns,nc,sys%B,rhs,lam,info)
         if (info.ne.0) call die('Least squares for lambda initialization failed.')
         ! Set the initial solution vector
         x(1:nc)=lam
         x(nc+1:nc+np)=log(Nbar)
      end block x_init

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use messager, only: die
      implicit none
      integer  :: i,j,iJ,jJ,isc,info
      integer  :: iter,iter_max
      real(WP) :: err,tol
      real(WP), dimension(:,:), allocatable :: Jac,BTB,PTP,BTP
      real(WP), dimension(:,:), allocatable :: Btilde,Ptilde
      real(WP), dimension(:),   allocatable :: y

      ! Allocate arrays
      allocate(Jac(nc+np,nc+np)); Jac=0.0_WP
      allocate(BTB(nc,nc));       BTB=0.0_WP
      allocate(PTP(np,np));       PTP=0.0_WP
      allocate(Btilde(ns,nc));    Btilde=0.0_WP
      allocate(Ptilde(ns,np));    Ptilde=0.0_WP
      allocate(BtildeT(nc,ns));   BtildeT=0.0_WP
      allocate(PtildeT(np,ns));   PtildeT=0.0_WP
      allocate(y(ns));            y=0.0_WP
      
      ! Newton-Raphson
      iter_max=100
      tol=1e-6
      y=get_y(x,gort)
      do iter=1,iter_max
         ! Remember the old solution
         x0=x
         ! Build the Jacobian matrix
         do j=1,nc
            Btilde(:,j)=y*sys%B(:,j)
         end do
         do j=1,np
            Ptilde(:,j)=y*sys%P(:,j)
         end do
         BtildeT=transpose(Btilde)
         PtildeT=transpose(Ptilde)
         do j=1,nc
            Btilde(:,j)=y*sys%B(:,j)
         end do
         do j=1,np
            Ptilde(:,j)=y*sys%P(:,j)
         end do
         BtildeT=transpose(Btilde)
         PtildeT=transpose(Ptilde)
         BTB=matmul(BtildeT,Btilde)
         PTP=matmul(PtildeT,Ptilde)
         BTP=matmul(BtildeT,Ptilde)
         do j=1,nc
            jJ=j
            do i=1,j
               iJ=i
               Jac(iJ,jJ)=BTB(iJ,jJ)
            end do
         end do
         do j=1,np
            jJ=nc+j
            do i=1,nc
               iJ=i
               Jac(iJ,jJ)=BTP(iJ,jJ)
            end do
         end do
         do j=1,np
            jJ=nc+j
            do i=1,j
               iJ=nc+i
               Jac(iJ,jJ)=PTP(iJ,jJ)
            end do
         end do
         do j=1,nc+np-1
            do i=j+1,nc+np
               Jac(i,j)=Jac(j,i)
            end do
         end do
         do i=1,np
            iJ=nc+i
            jJ=nc+i
            Jac(iJ,jJ)=Jac(iJ,jJ)-Nbar(i)
         end do
         ! Evaluate the residual error
         R=get_res(y)
         err=norm2(R)
         if (err.lt.tol) exit
         ! Solve for dx
         dx=-R
         call dpotrf('L',nc+np,Jac,nc+np,info)
         if (info.ne.0) call die('Cholesky factorization failed')
         call dpotrs('L',nc+np,1,Jac,nc+np,dx,nc+np,info)
         if (info.ne.0) call die('Linear solver failed')
         ! Update the solution
         x=x0+dx
         ! Get the species and phase moles
         y=get_y(x,gort)
         N=y*y
         Nbar=exp(x(nc+1:nc+np))
      end do

      ! Output
      print*,'Number of iterations = ', iter
      print*,'Residal error = ', err
      print*,'Equilibrium moles:'
      do isc=1,ns
         print*,trim(sp_names(isc)),': ',N(isc)
      end do

      ! Deallocate arrays
      deallocate(Jac,BTB,PTP,Btilde,Ptilde,y)

   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects-need destructors
      deallocate(species,sp_names,e_names,elem_mat,phse_mat)
      deallocate(N,Nbar,R,constraints,x,x0,dx,gort,BtildeT,PtildeT)

   end subroutine simulation_final


   
   
end module simulation
