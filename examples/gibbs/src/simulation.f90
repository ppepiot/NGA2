!> Example to read in YAML thermodynamic file
module simulation
   use precision, only: WP
   use YAMLRead,  only: YAMLElement
   use string,    only: str_short,str_medium
   implicit none
   private
   
   ! The array of the species. Eeach stored as a YAMLElement object
   type(YAMLElement), dimension(:), allocatable :: species

   ! Number of species and elements
   integer :: ns,ne

   ! Species and elements names
   character(len=str_medium), dimension(:), allocatable :: sp_names
   character(len=str_short),  dimension(:), allocatable :: e_names

   ! Elements matrix
   real(WP), dimension(:,:), allocatable :: elem_mat

   ! Simulation subroutines
   public :: simulation_init,simulation_run,simulation_final
   
contains


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_getsize
      implicit none
      real(WP), allocatable :: nasa_coef(:,:)
      ! Parse the mechanism file
      parse_mech: block
         use messager, only: die
         use YAMLRead, only: YAMLHandler,YAMLSequence,YAMLMap,yaml_open_file,yaml_start_from_sequence,yaml_close_file
         character(len=str_medium) :: yaml_file
         character(len=str_short), dimension(:), allocatable :: sp_names_copy
         type(YAMLHandler)  :: domain
         type(YAMLSequence) :: sp_list,phases,elements
         type(YAMLElement)  :: sp,gas
         type(YAMLMap)      :: thermo,comp
         integer :: isc,n,i,e,code
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
         n=0
         do isc=0,sp_list%size-1 ! Index in YAMLSequence starts from 0
            sp=sp_list%element(isc)
            name_arr=sp%value_str('name',code)
            name=''
            do i=1,size(name_arr)
               name=trim(name//name_arr(i))
            end do
            if (any(sp_names_copy.eq.name)) then
               n=n+1
               species(n)=sp
               sp_names(n)=name
            end if
            call sp%destroy()
         end do
         if(n.ne.ns) call die('Some species are missing in the mechanism file.')
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
               do n=1,size(comp%labels(i)%str)
                  name=name//trim(comp%labels(i)%str(n))
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
         ! Close the mechanism file and clean up
         call yaml_close_file(domain)
         call sp_list%destroy()
         call sp%destroy()
         call comp%destroy()
         call thermo%destroy()
      end block parse_mech

      ! Perform chemical equilibrium analysis
      ceq: block
         use ceq_types,   only: sys_type
         use ceq_system,  only: ceq_sys_init
         use ceq_state_m, only: ceq_state
         use messager,    only: die
         type(sys_type), pointer :: sys
         integer :: ncs=0,ng=0
         integer :: lu,iret,info
         real(WP), dimension(:,:), allocatable :: Bg
         integer,  dimension(:),   allocatable :: CS
         real(WP), dimension(:),   allocatable :: N_init,N,c,stats
         real(WP) :: T,HoR
         integer :: isc
         ! Allocate memory
         allocate(CS(ncs))
         allocate(Bg(ns,ng))
         allocate(N(ns))
         allocate(N_init(ns)); N_init=[2.0_WP,0.0_WP,1.0_WP,0.0_WP]
         allocate(c(ng))
         allocate(stats(20))
         ! Print the initial conditions
         print*,'Initial number of moles:'
         do isc=1,ns
            print*,trim(sp_names(isc)),': ',N_init(isc)
         end do
         ! Inizialize the system
         call ceq_sys_init(ns=ns,ne=ne,ncs=ncs,ng=ng,Ein=elem_mat,CS=CS,Bg=Bg,thermo_in=nasa_coef,lu_op=lu,diag=5,sys=sys,iret=iret)
         if (iret.lt.0) call die('System initialization failed.')
         ! Get the equilibrium state
         call ceq_state(sys=sys,N=N_init,p_Pa=101325.0_WP,T=300.0_WP,N_eq=N,T_eq=T,HoR_eq=HoR,stats=stats,info=info)
         if (info.lt.0) call die('ceq_state failed.')
         ! Print the solution
         print*,'Final number of moles:'
         do isc=1,ns
            print*,trim(sp_names(isc)),': ',N(isc)
         end do
      end block ceq

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none

      ! Deallocate
      ! deallocate(species,sp_names)
      
   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects-need destructors

   end subroutine simulation_final


   
   
end module simulation
