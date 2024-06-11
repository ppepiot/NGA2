!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Simulation monitor file
   type(monitor) :: mfile
   
   !> Routines in simulation.f90
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays and variables
   ! Reactor characteristics
   real(WP) :: tend,vol,time,tau
   integer :: niter,nout

   ! Current number of molecules
   real(WP), dimension(:), allocatable :: N,C,RR,RP

   ! Reaction stoichiometry
   real(WP), dimension(:,:), allocatable :: stoic

   ! Chemistry model parameters
   integer :: nreac, nspec
   real(WP), dimension(:), allocatable :: kmacro, k
   integer, dimension(:), allocatable :: nor, delta

   ! Temporary arrays
   real(WP), dimension(:), allocatable :: tmpR
   integer, dimension(:), allocatable :: itmpR
  
   ! Avogadro constant
   real(WP), parameter :: NA = 6.02214076e23_WP 
      

contains   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Initial and current concentrations
      real(WP), dimension(:), allocatable :: Cinit
  
      ! Local variables
      integer :: i
      
      ! Read in input sizing parameters
      call param_read('KMC volume',vol)
      call param_read('Number of species', nspec)
      call param_read('Number of reactions', nreac)

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(kmacro(nreac))
         allocate(k(nreac))
         allocate(delta(nreac))
         allocate(nor(nreac))
         allocate(Cinit(nspec))
         allocate(N(nspec))
         allocate(C(nspec))
         allocate(RR(nreac))
         allocate(RP(nreac))
         allocate(tmpR(nreac),itmpR(nreac))
      end block allocate_work_arrays

      ! Read in rest of parameters
      call param_read('Reaction rates',kmacro)
      call param_read('Reaction orders',nor)
      call param_read('Reaction type',delta)
      call param_read('Initial concentrations',Cinit)
      call param_read('Max run time',tend)
      call param_read('Output iteration period',nout)

      ! Allocate reaction stoichiometry array
      allocate(stoic(nreac,nspec))
      stoic(1,1) = -1
      stoic(1,2) = 1
      stoic(1,3) = 0
      stoic(2,1) = 0
      stoic(2,2) = -1
      stoic(2,3) = 1

      ! Calculate number of molecules
      N = Cinit*vol*NA
      C = Cinit

      ! Convert macro reaction rates to micro ones
      do i=1,nreac
         k(i) = kmacro(i)*delta(i)/((vol*NA)**(nor(i)-1))
      end do

      ! Initialize time variables
      niter = 0
      time = 0.0_WP 
      tau = 0.0_WP
      
      ! Create a monitor file
      create_monitor: block
         use parallel, only : amRoot
         ! Create simulation monitor
         mfile=monitor(amRoot,'simulation')
         call mfile%add_column(niter,'Timestep number')
         call mfile%add_column(time,'Time')
         call mfile%add_column(tau,'Timestep size')
         call mfile%add_column(C(1),'N1')
         call mfile%add_column(C(2),'N2')
         call mfile%add_column(C(3),'N3')
         call mfile%write()
      end block create_monitor
      
   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
   use quicksort
      implicit none
         
      ! Local variables
      real(WP) :: RRsum,rand
      integer :: i,ir

      ! Loop in time
      niter = 0
      do while (time<tend .and. N(1)>0)
  
         ! Calculate current reaction rates
         RR(1) = k(1)*N(1)
         RR(2) = k(2)*N(2)
         RRsum = sum(RR)
         RP = RR/RRsum

         ! Stochastic reaction time
         call random_number(rand)
         tau = log(1.0_WP/rand)/RRsum

         ! Selection of reaction channel
         do i=1,nreac
            itmpR(i) = i
            tmpR(i) = RR(i)
         end do
         ! Sort reactions by increasing probability to happen
         call quick_sort(tmpR,itmpR) 
         ! Cumulative reaction rate array
         do i=2,nreac
            tmpR(i) = sum(tmpR(1:i))
         end do
         ! random number generation for channel selection
         call random_number(rand)
         ! Channel selection
         i = 1
         do while (tmpR(i).lt.rand*RRsum)
            i = i+1
         end do
         ! Save selected channel
         ir = itmpR(i)

         ! Advance reactor one step following selected reaction
         time = time + tau
         do i=1,nspec
            N(i) = N(i)+stoic(ir,i)
         end do
         C = N/(vol*NA)

         ! Increment number of iteration
         niter = niter+1

         ! Update monitor file every nout iterations
         if (modulo(niter,nout) == 0) call mfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      !deallocate()
      
   end subroutine simulation_final
   
end module simulation
