!> Example to read in YAML thermodynamic file
module simulation
   use precision, only: WP
   use YAMLRead,  only: YAMLElement
   implicit none
   private
   
   ! The array of the species. Each stored as a YAMLElement object
   type(YAMLElement), dimension(:), allocatable :: species

   ! Simulation subroutines
   public :: simulation_init,simulation_run,simulation_final
   
contains


   !> Initialization of problem solver
   subroutine simulation_init
      use param,    only: param_read,param_getsize
      use string,   only: str_short,str_medium
      use messager, only: die
      use YAMLRead, only: YAMLHandler,YAMLSequence,yaml_open_file,yaml_start_from_sequence,yaml_close_file
      implicit none
      character(len=str_medium) :: yaml_file
      character(len=str_short), dimension(:), allocatable :: sp_names
      type(YAMLHandler)  :: domain
      type(YAMLSequence) :: sp_list
      type(YAMLElement)  :: sp
      real(WP), allocatable :: thermo_in(:,:)
      integer :: n_species,nsc,n,i,code
      character(len=:), allocatable :: name_arr(:)
      character(len=:), allocatable :: name

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

      ! Close the mechanism file
      call yaml_close_file(domain)

      ! Read the NASA-7 polynomials
      allocate(thermo_in(n_species,15))
      thermo_in=0.0_WP

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none

      ! Nothing to do
      
   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects-need destructors
      
   end subroutine simulation_final


   
   
end module simulation
