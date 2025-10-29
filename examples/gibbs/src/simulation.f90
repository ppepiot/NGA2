!> Example to read in YAML thermodynamic file
module simulation
   use precision,        only: WP
   use string,           only: str_short,str_medium
   use YAMLRead,         only: YAMLElement
   use chem_sys_class,   only: chem_sys,Lphase,Gphase
   use chem_state_class, only: chem_state,fixed_PT,fixed_PH
   implicit none
   private
   
   !> The array of the species. Eeach stored as a YAMLElement object
   type(YAMLElement), dimension(:), allocatable :: species

   !> Species names
   character(len=str_medium), dimension(:), allocatable :: sp_names
   integer, dimension(:), allocatable :: inpt2mch_sp_order

   !> Thermodynamic quantities
   real(WP) :: T,p
   real(WP), dimension(:), allocatable :: N_init

   !> Chemical system and state
   type(chem_sys)   :: sys
   type(chem_state) :: state

   !> Problem definition
   logical  :: scale
   real(WP) :: Nsum

   !> Simulation subroutines
   public :: simulation_init,simulation_run,simulation_final
   

contains


   !> Initialization of problem solver
   subroutine simulation_init
      use param,     only: param_exists,param_read,param_getsize
      use messager,  only: die
      implicit none
      integer :: np=2,ns,ne,ncs
      character(len=str_short), dimension(:), allocatable :: e_names
      real(WP), dimension(:,:), allocatable :: elem_mat
      real(WP), dimension(:,:), allocatable :: phse_mat
      real(WP), allocatable :: nasa_coef(:,:)
      character(len=str_medium), dimension(:), allocatable :: const_sp
      integer,  dimension(:), allocatable :: CS

      ! Parse the mechanism file
      parse_mech: block
         use mathtools,      only: reorder_rows
         use chem_sys_class, only: ncof
         use YAMLRead,       only: YAMLHandler,YAMLSequence,YAMLMap,yaml_open_file,yaml_start_from_sequence,yaml_close_file
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
         if (param_exists('Constrained species')) then
            ncs=param_getsize('Constrained species')
            allocate(const_sp(1:ncs))
            allocate(const_sp_copy(1:ncs))
            call param_read('Constrained species',const_sp)
            const_sp_copy=const_sp
         else
            ncs=0
         end if
         allocate(sp_names(1:ns))
         allocate(sp_names_copy(1:ns))
         allocate(CS(ncs))
         call param_read('Species',sp_names)
         sp_names_copy=sp_names
         ! Get the target species from input
         nn=param_getsize('Initial moles')
         if (ns.ne.nn) call die('Unequal number of species and moles in the input file.')
         allocate(inpt2mch_sp_order(1:ns))
         allocate(N_init(1:ns))
         allocate(N_init_copy(1:ns))
         call param_read('Initial moles',N_init)
         N_init_copy=N_init
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
                  inpt2mch_sp_order(nn)=i
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
         call reorder_rows(N_init_copy,inpt2mch_sp_order,N_init)
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
         allocate(nasa_coef(1:ns,2*ncof+1)); nasa_coef=0.0_WP
         allocate(a(1:2,1:ncof)); a=0.0_WP
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
            nasa_coef(isc,1)=T_range(2)
            nasa_coef(isc,2:  ncof+1)=a(1,:)
            nasa_coef(isc,9:2*ncof+1)=a(2,:)
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
         deallocate(sp_names_copy,N_init_copy)
         if (allocated(const_sp_copy)) deallocate(const_sp_copy)
      end block parse_mech

      ! Initialize the chemical equilibrium framework
      ceq_init: block
         use param,            only: param_exists
         use mathtools,        only: reorder_rows
         use messager,         only: die
         use chem_state_class, only: BS,NR,FD,LS
         integer :: ng=1,PH_method,dNdT_method
         real(WP), dimension(:,:), allocatable :: Bg
         real(WP), dimension(:),   allocatable :: N_h,N_h_c
         real(WP) :: T_h,T_g
         character(len=2) :: eq_cond,PH_alg,dNdT_alg
         integer :: isc
         ! Read inputs
         call param_read('Temperature',T)
         call param_read('Pressure',p)
         call param_read('Equilibrium condition',eq_cond)
         call param_read('Scale mole numbers',scale)
         ! Allocate arrays
         allocate(Bg(ns,ng));  Bg=0.0_WP
         allocate(N_h(ns))
         allocate(N_h_c(ns))
         ! Create the general constraints
         do isc=1,ns
            if (sp_names(isc).eq.'H2O')    Bg(isc,1)=1.0_WP
            if (sp_names(isc).eq.'H2O(L)') Bg(isc,1)=1.0_WP
         end do
         ! Print the initial conditions
         print*,'Initial moles:'
         do isc=1,ns
            print*,trim(sp_names(isc)),': ',N_init(isc)
         end do
         ! Inizialize the chemical system
         call sys%initialize(np=np,ns=ns,ne=ne,ncs=ncs,ng=ng,P=phse_mat,Ein=elem_mat,CS=CS,Bg=Bg,thermo_in=nasa_coef,diag=5)
         ! Initialize the chemical state
         Nsum=1.0_WP
         if (scale) then
            Nsum=sum(N_init)
            N_init=N_init/Nsum
         end if
         select case (eq_cond)
            case ('PT')
               call state%initialize(sys=sys,cond=fixed_PT,p=p)
               call state%N_init(T=T,N=N_init)
            case ('PH')
               call param_read('PH algorithm',PH_alg)
               if (PH_alg.eq.'BS') then
                  PH_method=BS
               else if (PH_alg.eq.'NR') then
                  PH_method=NR
                  call param_read('dNdT algorithm',dNdT_alg)
                  if (dNdT_alg.eq.'FD') then
                     dNdT_method=FD
                  else if (dNdT_alg.eq.'LS') then
                     dNdT_method=LS
                  else
                     call die('Wrong dNdT method')
                  end if
               else
                  call die('Wrong PH method')
               end if
               call state%initialize(sys=sys,cond=fixed_PH,PH_method=PH_method,dNdT_method=dNdT_method,p=p)
               call param_read('Temperature for enthalpy calculation',T_h)
               call param_read('Composition for enthalpy calculation',N_h)
               N_h_c=N_h
               call reorder_rows(N_h_c,inpt2mch_sp_order,N_h)
               if (scale) N_h=N_h/Nsum
               if (param_exists('Temperature initial guess')) then
                  call param_read('Temperature initial guess',T_g)
                  call state%N_init(N=N_init,N_h=N_h,T_h=T_h,T_g=T_g)
                  ! call state%N_init(N=N_init,HoR=-33690.099223276899_WP,T_g=T_g)
                  print*,'state%HoR = ',state%HoR
               else
                  call state%N_init(N=N_init,N_h=N_h,T_h=T_h)
               end if
            case default
               call die('Equilibrium condition must be either PT or PH')
         end select
         if (.not.state%success) call die('chem state N_init failed')
         print*,'Equilibrium condition: Constant ',eq_cond
         ! Read in numerical inputs
         call param_read('Newton tolerance',state%tol_N)
         call param_read('Newton max iterations',state%iter_N_max)
         if (state%cond.eq.fixed_PH) then
            if (state%PH_method.eq.BS) then
               call param_read('H tolerance',state%tol_H)
            else if (state%PH_method.eq.NR) then
               call param_read('T tolerance',state%tol_T)
            end if
            call param_read('T max iterations',state%iter_T_max)
         end if
         ! Re-initialization of moles
         if (scale) then 
            print*,'Re-initialization of moles (Scaled):'
         else
            print*,'Re-initialization of moles:'
         end if
         do isc=1,sys%ns
            print*,trim(sp_names(isc)),': ',state%N(isc)
         end do
         ! Deallocate arrays
         deallocate(Bg,N_h,N_h_c)
      end block ceq_init

      ! Deallocate arrays
      deallocate(nasa_coef,CS,e_names,elem_mat,phse_mat)
      if (allocated(const_sp)) deallocate(const_sp)
      
   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use messager, only: die
      use chem_state_class, only: BS,NR
      implicit none
      integer :: isc

      ! Obtain the chemical equilibrium state
      if (state%PH_method.eq.BS) then
         state%Tlo=350.0_WP
         state%Thi=450.0_WP
      end if
      call state%equilibrate()
      if (.not.state%success) call die('chem state equilibrate failed')

      ! Output
      if (state%cond.eq.fixed_PH) then
         print*,'Number of temperature iterations = ',state%iter_T
         if (state%PH_method.eq.NR) then
            print*,'Relative residual error of T = ',state%dT/state%T
         else if (state%PH_method.eq.BS) then
            print*,'Relative residual error of H = ',state%RH/state%HoR
         end if
      else
         print*,'Number of Newton iterations = ',state%iter_N
         print*,'Residal error = ', norm2(state%R)
      end if
      print*,'Equilibrium temperature = ',state%T,' (K)'
      print*,'Equilibrium moles:'
      if (scale) state%N=Nsum*state%N
      do isc=1,sys%ns
         print*,trim(sp_names(isc)),': ',state%N(isc)
      end do

   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects-need destructors
      deallocate(species,sp_names,inpt2mch_sp_order)

   end subroutine simulation_final

   
end module simulation
