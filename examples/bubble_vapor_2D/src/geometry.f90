!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: geometry_init,Lz

   real(WP) :: Lz

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
            Lz=Lx/nx
            allocate(x(nx+1))
            do i=1,nx+1
               x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
            end do
         else
            call param_read('Lu',Lu); call param_read('nu',nu); call param_read('r',r)
            allocate(xu(nu+1))
            do i=1,nu+1
               xu(i)=real(i-1,WP)/real(nu,WP)*Lu-0.5_WP*Lu
            end do
            xL=xu(nu+1)
            dx=xL-xu(nu)
            nn=0
            do while (xL.le.0.5_WP*Lx)
               nn=nn+1
               dx=r*dx
               xL=xL+dx
            end do
            Lx=2.0_WP*xL
            Lz=Lu/nu
            nx=nu+2*nn
            allocate(x(nx+1))
            x=0.0_WP
            x(nn+1:nn+nu+1)=xu
            dx=xu(nu+1)-xu(nu)
            do i=1,nn
               ir=nn+nu+1+i
               il=nn+1-i
               dx=r*dx
               x(ir)=x(ir-1)+dx
               x(il)=x(il+1)-dx
            end do
         end if
         ny=nx; nz=1
         allocate(y(ny+1)); allocate(z(nz+1))
         y=x
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.true.,name='bubble_vapor')
         
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