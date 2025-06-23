!> Example to read in YAML thermodynamic file
module simulation
   use precision, only: WP
   use YAMLRead,  only: YAMLElement
   use string,    only: str_short
   implicit none
   private
   
   ! The array of the species. Eeach stored as a YAMLElement object
   type(YAMLElement), dimension(:), allocatable :: species

   ! Simulation subroutines
   public :: simulation_init,simulation_run,simulation_final
   
contains


   !> Initialization of problem solver
   subroutine simulation_init
      use param,    only: param_read,param_getsize
      use string,   only: str_short,str_medium
      use messager, only: die
      use YAMLRead, only: YAMLHandler,YAMLSequence,YAMLMap,yaml_open_file,yaml_start_from_sequence,yaml_close_file
      implicit none
      character(len=str_medium) :: yaml_file
      character(len=str_short), dimension(:), allocatable :: sp_names
      type(YAMLHandler)  :: domain
      type(YAMLSequence) :: sp_list
      type(YAMLElement)  :: sp
      type(YAMLMap)      :: thermo
      real(WP), allocatable :: thermo_in(:,:)
      integer :: n_species,nsc,n,i,code
      character(len=:), allocatable :: name_arr(:)
      character(len=:), allocatable :: name
      real(WP), allocatable :: T_range(:)
      real(WP), dimension(:,:), allocatable :: a

      ! Get the target species from input
      n_species=param_getsize('Species')
      allocate(sp_names(1:n_species))
      call param_read('Species',sp_names)
      
      ! Read the mechanism file
      call param_read('YAML file',yaml_file)

      ! Get all the species from the mechanism
      domain=yaml_open_file(trim(yaml_file))
      sp_list=yaml_start_from_sequence(domain,'species')

      ! Extract the target species from the mechanism
      allocate(species(1:n_species))
      n=0
      do nsc=0,sp_list%size-1 ! Index in YAMLSequence starts from 0
         sp=sp_list%element(nsc)
         name_arr=sp%value_str('name',code)
         name=''
         do i=1,size(name_arr)
            name=trim(name//name_arr(i))
         end do
         if (any(sp_names.eq.name)) then
            n=n+1
            species(n)=sp
         end if
         call sp%destroy()
      end do
      if(n.ne.n_species) call die('Some species are missing in the mechanism file.')

      ! Read the NASA-7 polynomials
      allocate(thermo_in(n_species,15))
      allocate(a(2,7)); a=0.0_WP
      thermo_in=0.0_WP
      do nsc=1,n_species
         sp=species(nsc)
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
         thermo_in(nsc,1   )=T_range(2)
         thermo_in(nsc,2:8 )=a(1,:)
         thermo_in(nsc,9:15)=a(2,:)
      end do

      ! Close the mechanism file and clean up
      call yaml_close_file(domain)
      call sp_list%destroy()
      call sp%destroy()
      call thermo%destroy()

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none

      ! Deallocate
      deallocate(species)
      
   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects-need destructors
      
   end subroutine simulation_final


   
   
end module simulation
