!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: geometry_init

contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read,param_exists
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,ir,il,k,nx,ny,nz,nu,nn
         real(WP) :: Lx,Lu,r,xL,dx,my_x
         real(WP), dimension(:), allocatable :: x,y,z,xu
         
         ! Read in grid definition
         call param_read('L',Lx)
         
         ! Create simple rectilinear grid
         if (param_exists('n')) then
            call param_read('n',nx)
            allocate(x(nx+1))
            do i=1,nx+1
               x(i)=real(i-1,WP)/real(nx,WP)*Lx
            end do
         end if
         ny=nx; nz=nx
         allocate(y(ny+1)); allocate(z(nz+1))
         y=x; z=x

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='bubble_vapor_3D_sym')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_cfg


      ! Create masks for this config
      create_walls: block
         use mathtools, only: twoPi
         integer :: i,j,k
         cfg%VF=1.0_WP
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry