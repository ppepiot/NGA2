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
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid

      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         use param,       only: param_read
         integer :: i,j,k,nx,ny,nz,iter
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         real(WP), dimension(:), allocatable :: dx,dy,dz
         real(WP) :: Rmax
         real(WP), dimension(3) :: center,axis
         real(WP) :: xRefine, xpos
         real(WP) :: ypos, yrefine1, yrefine2, yrefine3, yrefine4
         real(WP) :: zpos, zrefine1, zrefine2, zrefine3, zrefine4

         real(WP) :: variance

         logical :: zper

         call param_read("variance",variance)

         ! Read in rotor disk properties
         call param_read('Rotor Disk max radius',Rmax)
         call param_read('Rotor Disk center',center)
         call param_read('Rotor Disk axis',axis)

         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1)); allocate(dx(nx))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1)); allocate(dy(ny))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1)); allocate(dz(nz))

         ! Define dx array, refine at rotor disk
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do iter=1, 3
            do i=1,nx
               xpos = (x(i)+x(i+1))/2
               dx(i) = (1 - 0.7_WP*exp(-variance*(xpos - center(1))**2))
            end do
            dx = dx/sum(dx)*Lx
            x(1) = 0.0_WP
            do i=2,nx+1
               x(i) = x(i-1) + dx(i-1)
            end do
         end do


         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do

         yrefine1 = center(2)-Rmax; yrefine2 = center(2)-0.5_WP*Rmax; yrefine3 = center(2)+0.5_WP*Rmax; yrefine4 = center(2) + Rmax
         zrefine1 = center(3)-Rmax; zrefine2 = center(3)-0.5_WP*Rmax; zrefine3 = center(3)+0.5_WP*Rmax; zrefine4 = center(3) + Rmax

         do iter=1, 3
            do j=1,ny
               ypos = (y(j)+y(j+1))/2
               dy(j) = 1 - 0.7_WP*exp(-variance*(ypos - yrefine1)**2) - 0.7_WP*exp(-variance*(ypos - yrefine4)**2)
            end do
            dy = dy/sum(dy)*Ly
            y(1) = -0.5_WP*Ly
            do j=2,ny+1
               y(j) = y(j-1) + dy(j-1)
            end do
         end do

         if (nz .gt. 1) then
            do iter=1, 3
               do k=1,nz
                  zpos = (z(k)+z(k+1))/2
                  dz(k) = 1 - 0.7_WP*exp(-variance*(zpos - zrefine1)**2) - 0.7_WP*exp(-variance*(zpos - zrefine4)**2)
               end do
               dz = dz/sum(dz)*Lz
               z(1) = -0.5_WP*Lz
               do k=2,nz+1
                  z(k) = z(k-1) + dz(k-1)
               end do
            end do
         end if

         zper = .false.
         if (nz .eq. 1) zper = .true.

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=zper,name='channel')


         deallocate(dx); deallocate(dy); deallocate(dz)

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
         ! Initialize to all wall
         cfg%VF=1.0_WP

         call cfg%sync(cfg%VF)


      end block create_walls

   end subroutine geometry_init


end module geometry
