!> Example to read in YAML thermodynamic file
module simulation
   use precision, only: WP
   !use ceq_types, only: sys_type
   use YAMLRead
   implicit none
   private
   
   ! YAML-related variables   
   character(len=:), allocatable :: fpath
   type(YAMLHandler) :: domain  ! be sure to close
   type(YAMLMap) :: parent
   type(YAMLMap) :: mapper_1
   type(YAMLMap) :: mapper_2
   integer :: n = -9999
   integer :: i = -9999
   integer :: code

   public :: simulation_init,simulation_run,simulation_final
   
contains

   !> Nothing to initialize
   subroutine simulation_init
      use param, only: param_read
      implicit none  

      call param_read('YAML file',fpath)
      domain = yaml_open_file(fpath)
      parent = yaml_start_from_map(domain, "parent")

      mapper_1 = parent%value_map("child1")
      mapper_2 = parent%value_map("child_22222")

      write(*,*) "mapper 1 = ", mapper_1%value_str("boolean_flag", code)
      write(*,*) "mapper 1 = ", mapper_1%value_double_1d("an_aray", code)
      write(*,*) "mapper 1 = ", mapper_1%value_int_2d("array_int_2d", code)
      write(*,*) ""
      write(*,*) "mapper 2 = ", mapper_2%value_str("boolean_flag", code)
      write(*,*) "mapper 2 = ", mapper_2%value_double_1d("an_aray", code)
      write(*,*) "mapper 2 = ", mapper_2%value_double_2d("array_float_2d", code)

      n = size(mapper_1%labels)
      write(*,*) ""//achar(13)//achar(10)//"Mapper 1 keys:"
      do i = 1,n
         write(*,*) "  ", mapper_1%labels(i)%str
      end do

      call mapper_1%destroy()

      n = size(mapper_2%labels)
      write(*,*) ""//achar(13)//achar(10)//"Mapper 2 keys:"
      do i = 1,n
         write(*,*) "  ", mapper_2%labels(i)%str
      end do

      call parent%destroy()
      call mapper_2%destroy()
      call yaml_close_file(domain)
      deallocate(fpath)

   end subroutine simulation_init


   !> NOthing here yet
   subroutine simulation_run
      implicit none

      ! Nothing to do
      
   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      
   end subroutine simulation_final


   
   
end module simulation
