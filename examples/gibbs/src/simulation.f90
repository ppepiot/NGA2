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

   !> Species and elements names
   character(len=str_medium), dimension(:), allocatable :: sp_names

   !> Thermodynamic quantities
   real(WP) :: T,p

   !> Chemical system and state
   type(chem_sys)   :: sys
   type(chem_state) :: state

   !> Simulation subroutines
   public :: simulation_init,simulation_run,simulation_final
   

contains


   !> Initialization of problem solver
   subroutine simulation_init
      use param,     only: param_read,param_getsize
      use messager,  only: die
      use mathtools, only: lss
      implicit none
      integer :: np=2,ns,ne,ncs
      character(len=str_short), dimension(:), allocatable :: e_names
      real(WP), dimension(:,:), allocatable :: elem_mat
      real(WP), dimension(:,:), allocatable :: phse_mat
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

      ! Initialize the chemical equilibrium framework
      ceq_init: block
         use messager, only: die
         integer :: ng=1
         real(WP), dimension(:,:), allocatable :: Bg
         character(len=2) :: eq_cond
         integer :: isc
         ! Read inputs
         call param_read('Temperature',T)
         call param_read('Pressure',p)
         call param_read('Equilibrium condition',eq_cond)
         ! Allocate arrays
         allocate(Bg(ns,ng));  Bg=0.0_WP
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
         select case (eq_cond)
            case ('PT')
               call state%initialize(sys=sys,cond=fixed_PT,p=p,T=T,N=N_init)
            case ('PH')
               call state%initialize(sys=sys,cond=fixed_PH,p=p,T=T,N=N_init,N_h=[1.0_WP,0.17260218510310191_WP,0.82739781489689534_WP,3.71_WP],T_h=300.0_WP)
            case default
               call die('Equilibrium condition must be either PT or PH')
         end select
         print*,'Equilibrium condition: Constant ',eq_cond
         ! Re-initialization of moles
         print*,'Re-initialization of moles:'
         do isc=1,sys%ns
            print*,trim(sp_names(isc)),': ',state%N(isc)
         end do
         ! Initialize the chemical state solution vector
         call state%x_init()
         ! Deallocate arrays
         deallocate(Bg)
      end block ceq_init

      ! Read in Newton and temperature iterations inputs
      call param_read('Newton tolerance',state%tol_N)
      call param_read('Newton max iterations',state%iter_N_max)
      if (state%cond.eq.fixed_PH) then
         call param_read('T tolerance',state%tol_T)
         call param_read('T max iterations',state%iter_T_max)
      end if

      ! Deallocate arrays
      deallocate(nasa_coef,const_sp,CS,N_init,e_names,elem_mat,phse_mat)

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use messager, only: die
      implicit none
      integer :: isc

      ! Obtain the chemical equilibrium state
      call state%equilibrate()

      ! Output
      if (state%cond.eq.fixed_PH) then
         print*,'Number of temperature iterations = ',state%iter_T
         print*,'Relative residual error of T = ',state%dT/state%T
      else
         print*,'Number of Newton iterations = ',state%iter_N
         print*,'Residal error = ', norm2(state%R)
      end if
      print*,'Equilibrium temperature = ',state%T,' (K)'
      print*,'Equilibrium moles:'
      do isc=1,sys%ns
         print*,trim(sp_names(isc)),': ',state%N(isc)
      end do

   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects-need destructors
      deallocate(species,sp_names)

   end subroutine simulation_final

   
end module simulation
