!> Liquid gas phase change class:
module lgpc_class
   use precision,         only: WP
   use string,            only: str_medium
   use config_class,      only: config
   use timetracker_class, only: timetracker
   use vfs_class,         only: vfs,VFlo,VFhi
   use tpscalar_class,    only: Lphase,Gphase
   use iterator_class,    only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: lgpc,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: symmetry=1                              !< Symmetry condition

   ! Index shift
   integer, dimension(1:3,1:3) :: ind_shift=reshape([1,0,0,0,1,0,0,0,1], shape(ind_shift))

   ! Default parameters for the lgpc class
   integer, parameter :: ext_lvl =5                                      !< The extension level for the interface normal (needs to be smaller than the VOF solver nband for the shift_mdot3p to work properly)
   integer, parameter :: extp_stc=2                                      !< The extrapolation stencil

   type :: arr_ptr_4d
      real(WP), pointer :: arr(:,:,:,:)
   end type arr_ptr_4d

   !> Boundary conditions for the extended normal field
   type :: bcond
      type(bcond), pointer :: next                                       !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'                  !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                                    !< Bcond type
      character(len=1) :: face                                           !< Bcond face (x/y/z)
      integer :: dir                                                     !< Bcond direction (-1 or +1)
      type(iterator) :: itr                                              !< This is the iterator for the bcond
   end type bcond

   !> lgpc object definition
   type :: lgpc

      ! This is our config
      class(config), pointer :: cfg                                      !< This is the config the object is build for
      
      ! This is the name of the object
      character(len=str_medium) :: name='UNNAMED_LGPC'                   !< Object name (default=UNNAMED_LGPC)

      ! This is the corresponding VOF solver
      class(vfs), pointer :: vf

      ! Liquid and gas densities
      real(WP) :: rho_l,rho_g

      ! Boundary condition list
      integer :: nbc                                                     !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                                   !< List of bcond for our solver

      ! Phase change mass flux data
      real(WP), dimension(:,:,:),   allocatable :: mdot2p                !< Phase change mass flux (m dot double prime of dimension M/(T*L^2))
      real(WP), dimension(:,:,:),   allocatable :: mdot3p                !< Volumetric phase change mass flux (m dot triple prime of dimension M/(T*L^3))
      real(WP), dimension(:,:,:,:), allocatable :: mdot3pLG,mdot3pLG_old !< Volumetric phase change mass flux spread out into the liquid and gas
      real(WP), dimension(:,:,:),   allocatable :: div_vel,div_vel_old   !< Phase change induced velocity divergence
      real(WP), dimension(:,:,:),   pointer     :: Tl,Tg                 !< The liquid and gas one-sided temperatures
      real(WP), dimension(:,:,:),   allocatable :: Tl_grd,Tg_grd         !< The liquid and gas one-sided temperature gradients

      ! Pseudo time over which the mdot3p is being shifted
      type(timetracker), public :: pseudo_time

      ! Phase change and interfacial velocity
      real(WP), dimension(:,:,:,:), allocatable :: normal                !< Interface normal vector
      real(WP), dimension(:,:,:,:), allocatable :: pseudo_vel            !< Pseudo velocity for shifting mdot3p

      ! Metrics (point to flow solver metrics)
      type(arr_ptr_4d), dimension(:), allocatable :: itp                 !< Cell to face interpolation coefficients
      type(arr_ptr_4d), dimension(:), allocatable :: div                 !< Divergence operator (cell-centered)

      ! Number of cells in each direction
      integer, dimension(3) :: nCell
      
      ! Monitoring quantities
      real(WP) :: mdot3p_int,mdot3p_tol                                  !< Integral and tolerence of the scaled lgpc mass flux
      real(WP) :: mdot3pL_int,mdot3pL_err,mdot3pL_int_err                !< Liquid side scaled lgpc mass flux maximum, integral, and error
      real(WP) :: mdot3pG_int,mdot3pG_err,mdot3pG_int_err                !< Gas side scaled lgpc mass flux maximum, integral, and error
      
   contains

      procedure :: initialize                                            !< Class initializer
      procedure :: init_mdot3pLG                                         !< Initialize lgpc mass flux on the liquid and gas sides
      procedure :: get_normal                                            !< Get the interface normal vector
      procedure :: get_pseudo_vel                                        !< Get the face-centered normilized gradient of VOF
      procedure :: get_mdot3p                                            !< Get the volumetric lgpc mass fllux
      procedure :: get_dmdot3pdt                                         !< Get the time derivative of the lgpc mass fllux
      procedure :: shift_mdot3p                                          !< Shift the lgpc mass flux
      procedure :: get_div                                               !< Get the lgpc source term
      procedure :: get_cfl                                               !< Get the CFL
      procedure :: extend_normal                                         !< Extend the interface normal                                             
      procedure :: pure_interfacial_extp                                 !< Extrapolate from pure to interfacial cells
      procedure :: get_grad                                              !< Get the Gauss gradient of a scalar field
      procedure :: add_bcond                                             !< Add a boundary condition
      procedure :: get_bcond                                             !< Get a boundary condition
      procedure :: apply_bcond                                           !< Apply all boundary conditions
      procedure :: get_temperature_grad                                  !< Get the temperature gradient
      procedure :: get_one_sided_grad                                    !<
      procedure :: itp_ccplane                                           !<
      procedure :: filter                                                !< Conservatively volume filter a filed that is defined at the interface
      procedure :: filterG                                               !< Conservatively volume filter a filed that is defined in the gas
      procedure :: filterL                                               !< Conservatively volume filter a filed that is defined in the liquid

   end type lgpc

   
contains
   
   
   !> Object initializer
   subroutine initialize(this,cfg,vf,SC,iTl,iTg,itp_x,itp_y,itp_z,div_x,div_y,div_z,name)
      use messager, only: die
      implicit none
      class(lgpc), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      class(vfs), target, intent(in) :: vf
      real(WP), target, dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:,1:), intent(in) :: SC
      integer, intent(in) :: iTl,iTg
      real(WP), target, dimension(-1: 0,cfg%imino_+1:cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: itp_x
      real(WP), target, dimension(-1: 0,cfg%imino_  :cfg%imaxo_,cfg%jmino_+1:cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: itp_y
      real(WP), target, dimension(-1: 0,cfg%imino_  :cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_+1:cfg%kmaxo_), intent(in) :: itp_z
      real(WP), target, dimension( 0:+1,cfg%imino_  :cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: div_x
      real(WP), target, dimension( 0:+1,cfg%imino_  :cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: div_y
      real(WP), target, dimension( 0:+1,cfg%imino_  :cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: div_z
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Check the extension level
      if (ext_lvl.ge.vf%nband) call die('The normal extension level needs to be smaller than the VOF solver band')

      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to pgrid object
      this%cfg=>cfg

      ! Point to VOF solver
      this%vf=>vf

      ! Allocate variables
      allocate(this%mdot2p      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));               this%mdot2p      =0.0_WP
      allocate(this%mdot3p      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));               this%mdot3p      =0.0_WP
      allocate(this%mdot3pLG    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,Lphase:Gphase)); this%mdot3pLG    =0.0_WP
      allocate(this%mdot3pLG_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,Lphase:Gphase)); this%mdot3pLG_old=0.0_WP
      allocate(this%div_vel     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));               this%div_vel     =0.0_WP
      allocate(this%div_vel_old (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));               this%div_vel_old =0.0_WP
      allocate(this%normal      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3));           this%normal      =0.0_WP
      allocate(this%pseudo_vel  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3));           this%pseudo_vel  =0.0_WP
      allocate(this%Tl_grd      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));               this%Tl_grd      =0.0_WP
      allocate(this%Tg_grd      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));               this%Tg_grd      =0.0_WP
      allocate(this%itp(3))
      allocate(this%div(3))

      ! Point to the temperature fields
      this%Tl(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)=>SC(:,:,:,iTl)
      this%Tg(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)=>SC(:,:,:,iTg)

      ! Set metrics
      this%itp(1)%arr(-1:,this%cfg%imino_+1:,this%cfg%jmino_  :,this%cfg%kmino_  :)=>itp_x
      this%itp(2)%arr(-1:,this%cfg%imino_  :,this%cfg%jmino_+1:,this%cfg%kmino_  :)=>itp_y
      this%itp(3)%arr(-1:,this%cfg%imino_  :,this%cfg%jmino_  :,this%cfg%kmino_+1:)=>itp_z
      this%div(1)%arr( 0:,this%cfg%imino_  :,this%cfg%jmino_  :,this%cfg%kmino_  :)=>div_x
      this%div(2)%arr( 0:,this%cfg%imino_  :,this%cfg%jmino_  :,this%cfg%kmino_  :)=>div_y
      this%div(3)%arr( 0:,this%cfg%imino_  :,this%cfg%jmino_  :,this%cfg%kmino_  :)=>div_z

      ! Number of cells in each direction
      this%nCell(1)=this%cfg%nx
      this%nCell(2)=this%cfg%ny
      this%nCell(3)=this%cfg%nz

      ! Create a pseudo time
      this%pseudo_time=timetracker(amRoot=this%cfg%amRoot,name='Pseudo time',print_info=.false.)

   end subroutine initialize


   !> Initialize the liquid and gas side volumetric mass fluxes
   subroutine init_mdot3pLG(this)
      implicit none
      class(lgpc), intent(inout) :: this
      ! Initialize the liquid and gas side mdot3p
      this%mdot3pLG(:,:,:,Lphase)=-this%mdot3p
      this%mdot3pLG(:,:,:,Gphase)= this%mdot3p
      ! Initialize errors to zero
      this%mdot3pL_err    =0.0_WP
      this%mdot3pG_err    =0.0_WP
      this%mdot3pL_int_err=0.0_WP
      this%mdot3pG_int_err=0.0_WP
      ! Get the integrals
      call this%cfg%integrate(this%mdot3p,this%mdot3p_int)
      call this%cfg%integrate(this%mdot3pLG(:,:,:,Lphase),this%mdot3pL_int)
      call this%cfg%integrate(this%mdot3pLG(:,:,:,Gphase),this%mdot3pG_int)
   end subroutine init_mdot3pLG
   

   !> Get the interfacial normal vector
   subroutine get_normal(this)
      use irl_fortran_interface, only: calculateNormal,getNumberOfVertices
      implicit none
      class(lgpc), intent(inout) :: this
      real(WP), dimension(3) :: n1,n2
      integer :: i,j,k,dir
      ! Loop over cells
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               n1=calculateNormal(this%vf%interface_polygon(1,i,j,k))
               if (getNumberOfVertices(this%vf%interface_polygon(2,i,j,k)).gt.0) then
                  n2=calculateNormal(this%vf%interface_polygon(2,i,j,k))
                  n1=0.5_WP*(n1+n2)
               end if
               do dir=1,3
                  this%normal(i,j,k,dir)=n1(dir)
               end do
            end do
         end do
      end do
   end subroutine get_normal


   !> Extend the interface normal for a smoother transition to zero
   subroutine extend_normal(this)
      implicit none
      class(lgpc), intent(inout) :: this
      integer  :: n,dir,i,j,k,index,index_pure
      integer  :: stx,sty,stz
      real(WP) :: vol
      real(WP), dimension(:,:,:,:), allocatable :: normal_tmp
      ! Allocate memory for the temporary field
      allocate(normal_tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:3))
      ! Loop over the extension levels
      do n=1,ext_lvl
         ! Take a copy of the normal
         normal_tmp=this%normal
         ! Loop over directions
         do dir=1,3
            ! Check the dimensions
            if (this%nCell(dir).gt.1) then
               ! Loop over the pure cells within the band
               do index=1,sum(this%vf%band_count(1:ext_lvl))
                  ! Offset for the interfacial cells' index
                  index_pure=this%vf%band_count(0)+index
                  ! Get the cell indices
                  i=this%vf%band_map(1,index_pure)
                  j=this%vf%band_map(2,index_pure)
                  k=this%vf%band_map(3,index_pure)
                  ! Initialize with zero
                  vol=0.0_WP
                  this%normal(i,j,k,dir)=0.0_WP
                  ! Loop over the stencils
                  do stz=-1,1
                     if (k+stz.lt.this%cfg%kmin.or.k+stz.gt.this%cfg%kmax) cycle
                     do sty=-1,1
                        if (j+sty.lt.this%cfg%jmin.or.j+sty.gt.this%cfg%jmax) cycle
                        do stx=-1,1
                           if (i+stx.lt.this%cfg%imin.and.i+stx.gt.this%cfg%imax) cycle
                           vol=vol+this%cfg%vol(i+stx,j+sty,k+stz)
                           this%normal(i,j,k,dir)=this%normal(i,j,k,dir)+this%cfg%vol(i+stx,j+sty,k+stz)*normal_tmp(i+stx,j+sty,k+stz,dir)
                        end do
                     end do
                  end do
                  ! Scale by volume
                  if (vol.gt.0.0_WP) this%normal(i,j,k,dir)=this%normal(i,j,k,dir)/vol
               end do
               ! Sync
               call this%cfg%sync(this%normal(:,:,:,dir))
            end if
         end do
      end do
      ! Deallocate the normal copy
      deallocate(normal_tmp)
   end subroutine extend_normal


   !> Calculate the pseudo velocity used to shift the lgpc mass flux
   subroutine get_pseudo_vel(this)
      implicit none
      class(lgpc), intent(inout) :: this
      integer :: i,j,k,dir
      integer :: im,jm,km
      integer :: ip,jp,kp
      ! Initialize with zeros
      this%pseudo_vel=0.0_WP
      ! Loop over directions
      do dir=1,3
         ! Skip if needed
         if (this%nCell(dir).eq.1) cycle
         ! Loop over cell faces
         do k=this%cfg%kmin_,this%cfg%kmax_+ind_shift(3,dir)
            do j=this%cfg%jmin_,this%cfg%jmax_+ind_shift(2,dir)
               do i=this%cfg%imin_,this%cfg%imax_+ind_shift(1,dir)
                  ! Prepare indices for the adjacent cells
                  im=i-ind_shift(1,dir); ip=i
                  jm=j-ind_shift(2,dir); jp=j
                  km=k-ind_shift(3,dir); kp=k
                  this%pseudo_vel(i,j,k,dir)=-this%itp(dir)%arr(-1,i,j,k)*this%normal(im,jm,km,dir) &
                  &                          -this%itp(dir)%arr( 0,i,j,k)*this%normal(ip,jp,kp,dir)
               end do
            end do
         end do
      end do
   end subroutine get_pseudo_vel
   

   !> Calculate the volumetric mass flux
   subroutine get_mdot3p(this)
      implicit none
      class(lgpc), intent(inout) :: this
      ! call this%filter(F=this%mdot2p,lvl=2,stc=3)
      this%mdot3p=this%mdot2p*this%vf%SD
      ! call this%filter(F=this%mdot3p,lvl=2,stc=3)
   end subroutine get_mdot3p

   
   !> Calculate the explicit mdot3p time derivative
   subroutine get_dmdot3pdt(this,vel,mdot3p_old,dmdot3pdt)
      implicit none
      class(lgpc), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(in)  :: vel        !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:3)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:),    intent(in)  :: mdot3p_old  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:),    intent(out) :: dmdot3pdt   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,dir,index,n
      integer :: ic,jc,kc
      integer :: im,jm,km
      integer :: ip,jp,kp
      real(WP), dimension(:,:,:,:), allocatable :: F
      ! Zero out dmdot3p/dt array
      dmdot3pdt=0.0_WP
      ! Allocate flux arrays
      allocate(F(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:3)); F=0.0_WP
      ! Fluxes of mdot3p
      do index=1,sum(this%vf%band_count(0:ext_lvl+1))
         ! Get the cell indices
         ic=this%vf%band_map(1,index)
         jc=this%vf%band_map(2,index)
         kc=this%vf%band_map(3,index)
         ! Loop over directions
         do dir=1,3
            ! Loop over the plus and minus faces
            do n=0,1
               ! Get the face indices
               i=ic+ind_shift(1,dir)*n
               j=jc+ind_shift(2,dir)*n
               k=kc+ind_shift(3,dir)*n
               ! Prepare indices for the adjacent cells
               im=i-ind_shift(1,dir); ip=i
               jm=j-ind_shift(2,dir); jp=j
               km=k-ind_shift(3,dir); kp=k
               ! Compute the face flux
               if (F(i,j,k,dir).eq.0.0_WP) then
                  F(i,j,k,dir)=-0.5_WP*(vel(i,j,k,dir)+abs(vel(i,j,k,dir)))*mdot3p_old(im,jm,km) &
                  &            -0.5_WP*(vel(i,j,k,dir)-abs(vel(i,j,k,dir)))*mdot3p_old(ip,jp,kp)
               end if
            end do
         end do
      end do
      ! Time derivative of mdot3p
      do dir=1,3
         ! Loop over the cells
         do index=1,sum(this%vf%band_count(0:ext_lvl+1))
            ! Get the cell indices
            i=this%vf%band_map(1,index)
            j=this%vf%band_map(2,index)
            k=this%vf%band_map(3,index)
            ! Get the residual
            dmdot3pdt(i,j,k)=dmdot3pdt(i,j,k)                                                                          &
            &              +this%div(dir)%arr(0,i,j,k)*F(i,j,k,dir)                                                    &
            &              +this%div(dir)%arr(1,i,j,k)*F(i+ind_shift(1,dir),j+ind_shift(2,dir),k+ind_shift(3,dir),dir)
         end do
      end do
      ! Deallocate flux arrays
      deallocate(F)
      ! Sync residual
      call this%cfg%sync(dmdot3pdt)
   end subroutine get_dmdot3pdt


   !> Shift mdot3p away from the interface
   subroutine shift_mdot3p(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(lgpc), intent(inout) :: this
      real(WP), dimension(:,:,:), allocatable :: resmdot3pL,resmdot3pG
      integer  :: i,j,k,ierr
      real(WP) :: my_mdot3p_int,mdot3p_err

      ! Allocate memory for mdot3p residuals
      allocate(resmdot3pL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(resmdot3pG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Initialize the lgpc mass fluxes on the liquid and gas sides
      this%mdot3pLG(:,:,:,Lphase)=-this%mdot3p
      this%mdot3pLG(:,:,:,Gphase)= this%mdot3p

      ! Get the interface normal
      call this%get_normal()

      ! Extend the interface normal
      call this%extend_normal()

      ! Apply boundary conditions
      call this%apply_bcond()

      ! Get the normalized gradient of VOF
      call this%get_pseudo_vel()

      ! Get the CFL based on the gradient of the VOF
      call this%get_cfl(this%pseudo_time%dt,this%pseudo_vel(:,:,:,1),this%pseudo_vel(:,:,:,2),this%pseudo_vel(:,:,:,3),this%pseudo_time%cfl)
      
      ! Reset the pseudo time
      call this%pseudo_time%reset()
      
      ! Adjust the pseudo time step
      call this%pseudo_time%adjust_dt()

      ! Move the lgpc mass flux away from the interface
      do while (.not.this%pseudo_time%done())
         
         ! Remember old mdot3p
         this%mdot3pLG_old=this%mdot3pLG
         
         ! Increment pseudo time
         call this%pseudo_time%increment()
         
         ! Assemble explicit residual
         call this%get_dmdot3pdt(vel= this%pseudo_vel,mdot3p_old=this%mdot3pLG_old(:,:,:,Lphase),dmdot3pdt=resmdot3pL)
         call this%get_dmdot3pdt(vel=-this%pseudo_vel,mdot3p_old=this%mdot3pLG_old(:,:,:,Gphase),dmdot3pdt=resmdot3pG)
         
         ! Apply these residuals
         this%mdot3pLG(:,:,:,Lphase)=this%mdot3pLG_old(:,:,:,Lphase)+this%pseudo_time%dt*resmdot3pL
         this%mdot3pLG(:,:,:,Gphase)=this%mdot3pLG_old(:,:,:,Gphase)+this%pseudo_time%dt*resmdot3pG

         ! Re-apply boundary conditions
         call this%apply_bcond()

         ! Calculate the integral of the residual error of mdot3pL
         my_mdot3p_int=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%vf%VF(i,j,k).lt.VFhi) my_mdot3p_int=my_mdot3p_int+this%cfg%vol(i,j,k)*this%cfg%VF(i,j,k)*this%vf%VF(i,j,k)*abs(this%mdot3pLG(i,j,k,Lphase))
               end do
            end do
         end do
         call MPI_ALLREDUCE(my_mdot3p_int,this%mdot3pL_err,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         
         ! Calculate the integral of the residual error of mdot3pG
         my_mdot3p_int=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%vf%VF(i,j,k).gt.VFlo) my_mdot3p_int=my_mdot3p_int+this%cfg%vol(i,j,k)*this%cfg%VF(i,j,k)*(1.0_WP-this%vf%VF(i,j,k))*abs(this%mdot3pLG(i,j,k,Gphase))
               end do
            end do
         end do
         call MPI_ALLREDUCE(my_mdot3p_int,this%mdot3pG_err,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)

         ! Check convergence
         mdot3p_err=max(this%mdot3pL_err,this%mdot3pG_err)
         if (mdot3p_err.lt.this%mdot3p_tol) exit

      end do

      ! Take conservative spatial filter
      ! call this%filterL(F=this%mdot3pLG(:,:,:,Lphase),lvl=2,stc=1)
      ! call this%filterG(F=this%mdot3pLG(:,:,:,Gphase),lvl=2,stc=1)

      ! Integral of mdot3p
      call this%cfg%integrate(this%mdot3p,this%mdot3p_int)
      call this%cfg%integrate(this%mdot3pLG(:,:,:,Lphase),this%mdot3pL_int)
      call this%cfg%integrate(this%mdot3pLG(:,:,:,Gphase),this%mdot3pG_int)
      this%mdot3pL_int_err=abs(abs(this%mdot3pL_int)-this%mdot3p_int)
      this%mdot3pG_int_err=abs(abs(this%mdot3pG_int)-this%mdot3p_int)

      ! Deallocate mdot3p residuals
      deallocate(resmdot3pL,resmdot3pG)

   end subroutine shift_mdot3p


   !> Extrapolate a scalar field from pure to interfacial cells
   subroutine pure_interfacial_extp(this,phase,A)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(lgpc), intent(inout) :: this
      integer, intent(in) :: phase
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout)  :: A  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer  :: i,j,k,ihat,jhat,khat,index,ierr
      real(WP), dimension(:,:,:), allocatable :: w
      real(WP) :: wsum,denum
      integer :: stxm,stym,stzm
      integer :: stxp,styp,stzp
      integer :: stx,sty,stz
      real(WP), dimension(3) :: d
      
      ! Check the dimensions
      if (this%nCell(1).gt.1) then
         stxm=-extp_stc; stxp=+extp_stc
      else
         stxm= 0; stxp= 0
      end if
      if (this%nCell(2).gt.1) then
         stym=-extp_stc; styp=+extp_stc
      else
         stym= 0; styp= 0
      end if
      if (this%nCell(3).gt.1) then
         stzm=-extp_stc; stzp=+extp_stc
      else
         stzm= 0; stzp= 0
      end if

      ! Allocate the weight matrix
      allocate(w(stxm:stxp,stym:styp,stzm:stzp))

      ! Get the interface normal
      call this%get_normal()

      ! Loop over the interfacial cells
      do index=1,this%vf%band_count(0)
         ! Get the interfacial cell indices
         ihat=this%vf%band_map(1,index)
         jhat=this%vf%band_map(2,index)
         khat=this%vf%band_map(3,index)
         ! Initialize weights
         w=0.0_WP
         ! Loop over the stencil
         do stz=stzm,stzp
            k=khat+stz
            do sty=stym,styp
               j=jhat+sty
               do stx=stxm,stxp
                  i=ihat+stx
                  ! Calculate the weight
                  if (this%vf%VF(i,j,k).eq.(1.0_WP-real(phase,WP))) then
                     d=[this%cfg%xm(ihat)-this%cfg%xm(i),this%cfg%ym(jhat)-this%cfg%ym(j),this%cfg%zm(khat)-this%cfg%zm(k)]
                     denum=sum(d*d)
                     ! w(stx,sty,stz)=abs(sum(d*this%normal(ihat,jhat,khat,:)))*denum
                     denum=denum*denum
                     if (denum.gt.0.0_WP) w(stx,sty,stz)=abs(sum(d*this%normal(ihat,jhat,khat,:)))/denum
                  end if
               end do
            end do
         end do
         ! Normalize the weights
         wsum=sum(W)
         if (wsum.gt.epsilon(wsum)) w=w/wsum
         ! Initialize with zero
         A(ihat,jhat,khat)=0.0_WP
         ! Loop over the stencil
         do stz=stzm,stzp
            k=khat+stz
            do sty=stym,styp
               j=jhat+sty
               do stx=stxm,stxp
                  i=ihat+stx
                  ! Update the interfacial value
                  if (w(stx,sty,stz).gt.0.0_WP) A(ihat,jhat,khat)=A(ihat,jhat,khat)+w(stx,sty,stz)*A(i,j,k)
               end do
            end do
         end do
      end do

      ! Sync the scalar field
      call this%cfg%sync(A)

      ! Deallocate weights
      deallocate(w)

   end subroutine pure_interfacial_extp


   !> Get the gradient of a scalar field
   subroutine get_grad(this,phase,A,A_grd)
      implicit none
      class(lgpc), intent(in) :: this
      integer, intent(in) :: phase
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: A     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: A_grd !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: index,index_pure,i,j,k
      real(WP) :: grdX,grdY,grdZ
      real(WP) :: grdm,grdp
      real(WP) :: Am,Ap,lm,lp,vof
      ! Initialize with zeros
      A_grd=0.0_WP
      ! Loop over the pure cells within the stencil
      do index=1,sum(this%vf%band_count(1:extp_stc))
         ! Offset for the interfacial cells' index
         index_pure=this%vf%band_count(0)+index
         ! Get the cell indices
         i=this%vf%band_map(1,index_pure)
         j=this%vf%band_map(2,index_pure)
         k=this%vf%band_map(3,index_pure)
         ! Calculate the Gauss gradient
         if (this%vf%VF(i,j,k).eq.(1.0_WP-real(phase,WP))) then
            ! ! X grad
            ! vof=this%vf%VF(i-1,j,k)
            ! if ((vof.gt.0.0_WP).and.(vof.lt.1.0_WP)) then
            !    lm=0.5_WP*this%cfg%dx(i)
            !    lp=this%cfg%dx(i-1)*abs(real(phase,WP)-vof)
            !    Am=(lm*A(i-1,j,k)+lp*A(i,j,k))/(lm+lp)
            ! else
            !    Am=this%itp(1)%arr(-1,i,j,k)*A(i-1,j,k)+this%itp(1)%arr(0,i,j,k)*A(i,j,k)
            ! end if
            ! vof=this%vf%VF(i+1,j,k)
            ! if ((vof.gt.0.0_WP).and.(vof.lt.1.0_WP)) then
            !    lm=this%cfg%dx(i+1)*abs(real(phase,WP)-vof)
            !    lp=0.5_WP*this%cfg%dx(i)
            !    Ap=(lm*A(i,j,k)+lp*A(i+1,j,k))/(lm+lp)
            ! else
            !    Ap=this%itp(1)%arr(-1,i+1,j,k)*A(i,j,k)+this%itp(1)%arr(0,i+1,j,k)*A(i+1,j,k)
            ! end if
            ! grdX=this%div(1)%arr(0,i,j,k)*Am+this%div(1)%arr(1,i,j,k)*Ap
            ! ! Y grad
            ! vof=this%vf%VF(i,j-1,k)
            ! if ((vof.gt.0.0_WP).and.(vof.lt.1.0_WP)) then
            !    lm=0.5_WP*this%cfg%dy(j)
            !    lp=this%cfg%dy(j-1)*abs(real(phase,WP)-vof)
            !    Am=(lm*A(i,j-1,k)+lp*A(i,j,k))/(lm+lp)
            ! else
            !    Am=this%itp(2)%arr(-1,i,j,k)*A(i,j-1,k)+this%itp(2)%arr(0,i,j,k)*A(i,j,k)
            ! end if
            ! vof=this%vf%VF(i,j+1,k)
            ! if ((vof.gt.0.0_WP).and.(vof.lt.1.0_WP)) then
            !    lm=this%cfg%dy(j+1)*abs(real(phase,WP)-vof)
            !    lp=0.5_WP*this%cfg%dy(j)
            !    Ap=(lm*A(i,j,k)+lp*A(i,j+1,k))/(lm+lp)
            ! else
            !    Ap=this%itp(2)%arr(-1,i,j+1,k)*A(i,j,k)+this%itp(2)%arr(0,i,j+1,k)*A(i,j+1,k)
            ! end if
            ! grdY=this%div(2)%arr(0,i,j,k)*Am+this%div(2)%arr(1,i,j,k)*Ap
            ! ! Z grad
            ! vof=this%vf%VF(i,j,k-1)
            ! if ((vof.gt.0.0_WP).and.(vof.lt.1.0_WP)) then
            !    lm=0.5_WP*this%cfg%dz(k)
            !    lp=this%cfg%dz(k-1)*abs(real(phase,WP)-vof)
            !    Am=(lm*A(i,j,k-1)+lp*A(i,j,k))/(lm+lp)
            ! else
            !    Am=this%itp(3)%arr(-1,i,j,k)*A(i,j,k-1)+this%itp(3)%arr(0,i,j,k)*A(i,j,k)
            ! end if
            ! vof=this%vf%VF(i,j,k+1)
            ! if ((vof.gt.0.0_WP).and.(vof.lt.1.0_WP)) then
            !    lm=this%cfg%dz(k+1)*abs(real(phase,WP)-vof)
            !    lp=0.5_WP*this%cfg%dz(k)
            !    Ap=(lm*A(i,j,k)+lp*A(i,j,k+1))/(lm+lp)
            ! else
            !    Ap=this%itp(3)%arr(-1,i,j,k+1)*A(i,j,k)+this%itp(3)%arr(0,i,j,k+1)*A(i,j,k+1)
            ! end if
            ! grdZ=this%div(3)%arr(0,i,j,k)*Am+this%div(3)%arr(1,i,j,k)*Ap

            A_grd(i,j,k)=((this%div(1)%arr(0,i,j,k)*(this%itp(1)%arr(-1,i,j,k)*A(i-1,j,k)+this%itp(1)%arr(0,i,j,k)*A(i,j,k))+this%div(1)%arr(1,i,j,k)*(this%itp(1)%arr(-1,i+1,j,k)*A(i,j,k)+this%itp(1)%arr(0,i+1,j,k)*A(i+1,j,k)))**2.0_WP &
            &            +(this%div(2)%arr(0,i,j,k)*(this%itp(2)%arr(-1,i,j,k)*A(i,j-1,k)+this%itp(2)%arr(0,i,j,k)*A(i,j,k))+this%div(2)%arr(1,i,j,k)*(this%itp(1)%arr(-1,i,j+1,k)*A(i,j,k)+this%itp(2)%arr(0,i,j+1,k)*A(i,j+1,k)))**2.0_WP &
            &            +(this%div(3)%arr(0,i,j,k)*(this%itp(3)%arr(-1,i,j,k)*A(i,j,k-1)+this%itp(3)%arr(0,i,j,k)*A(i,j,k))+this%div(3)%arr(1,i,j,k)*(this%itp(1)%arr(-1,i,j,k+1)*A(i,j,k)+this%itp(3)%arr(0,i,j,k+1)*A(i,j,k+1)))**2.0_WP)**0.5_WP

            ! grdm=this%cfg%dxmi(i  )*(A(i  ,j,k)-A(i-1,j,k))
            ! grdp=this%cfg%dxmi(i+1)*(A(i+1,j,k)-A(i  ,j,k))
            ! grdX=minmod(grdm,grdp)
            ! grdm=this%cfg%dymi(j  )*(A(i,j  ,k)-A(i,j-1,k))
            ! grdp=this%cfg%dymi(j+1)*(A(i,j+1,k)-A(i,j ,k))
            ! grdY=minmod(grdm,grdp)
            ! grdm=this%cfg%dzmi(k  )*(A(i,j,k  )-A(i,j,k-1))
            ! grdp=this%cfg%dzmi(k+1)*(A(i,j,k+1)-A(i,j,k  ))
            ! grdZ=minmod(grdm,grdp)

            ! A_grd(i,j,k)=sqrt(grdX**2.0_WP+grdY**2.0_WP+grdZ**2.0_WP)
         end if
      end do
      ! Sync it
      call this%cfg%sync(A_grd)

      contains
      
         !> Minmod gradient
         function minmod(g1,g2) result(g)
            implicit none
            real(WP), intent(in) :: g1,g2
            real(WP) :: g
            if (g1*g2.le.0.0_WP) then
               g=0.0_WP
            else
               if (abs(g1).lt.abs(g2)) then
                  g=g1
               else
                  g=g2
               end if
            end if
         end function minmod

   end subroutine get_grad


   !> Calculate the divergence induced by phase change
   subroutine get_div(this)
      implicit none
      class(lgpc), intent(inout) :: this
      this%div_vel=this%mdot3pLG(:,:,:,Gphase)/this%rho_g+this%mdot3pLG(:,:,:,Lphase)/this%rho_l
   end subroutine get_div


   !> Calculate temperature gradients on the liquid and gas sides (Boyd and Ling 2023)
   subroutine get_temperature_grad(this)
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(lgpc), intent(inout) :: this
      ! Apply boundary conditions
      call this%apply_bcond()
      ! Get the temperature gradients
      call this%get_grad(Lphase,this%Tl,this%Tl_grd)
      call this%get_grad(Gphase,this%Tg,this%Tg_grd)
      ! Extrapolate the gradients to interface
      call this%pure_interfacial_extp(Lphase,this%Tl_grd)
      call this%pure_interfacial_extp(Gphase,this%Tg_grd)
   end subroutine get_temperature_grad


   ! !> Calculate temperature gradients on the liquid and gas sides (Bothe and Fleckenstein 2013)
   ! subroutine get_temperature_grad(this)
   !    use irl_fortran_interface, only: calculateCentroid
   !    implicit none
   !    class(lgpc), intent(inout) :: this
   !    real(WP), dimension(3) :: posI
   !    integer  :: i,j,k,index
   !    ! Get the interface normal
   !    call this%get_normal()
   !    ! Apply boundary conditions
   !    call this%apply_bcond()
   !    ! Zero out the gradients
   !    this%Tl_grd=0.0_WP
   !    this%Tg_grd=0.0_WP
   !    ! Loop over the interfacial cells
   !    do index=1,this%vf%band_count(0)
   !       ! Get the interfacial cell indices
   !       i=this%vf%band_map(1,index)
   !       j=this%vf%band_map(2,index)
   !       k=this%vf%band_map(3,index)
   !       ! Get the interface center
   !       posI=calculateCentroid(this%vf%interface_polygon(1,i,j,k))
   !       ! Get the liquid side gradient
   !       call this%get_one_sided_grad(phase=Lphase,F=this%Tl,posI=posI,i=i,j=j,k=k,normal=this%normal(i,j,k,:),Fgrd=this%Tl_grd(i,j,k))
   !       ! Get the gas side gradient
   !       call this%get_one_sided_grad(phase=Gphase,F=this%TG,posI=posI,i=i,j=j,k=k,normal=this%normal(i,j,k,:),Fgrd=this%Tg_grd(i,j,k))
   !    end do
   ! end subroutine get_temperature_grad


   !> Interpolates the input scalar field (F) at the intersection between the input line (defined by a point and a normal vector)
   !> and the closest plane that goes through the cell centers of the grid
   subroutine itp_ccplane(this,F,x0,y0,i,j,k,nx,ny,Fitp,xitp,yitp)
      implicit none
      class(lgpc), intent(in) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: F
      real(WP), intent(in)  :: x0,y0
      integer,  intent(in)  :: i,j,k
      real(WP), intent(inout)  :: nx,ny
      real(WP), intent(out) :: Fitp
      real(WP), optional, intent(out) :: xitp,yitp
      real(WP) :: dx,dy
      real(WP) :: Sx,Sy
      real(WP) :: w_m,w_p,wsum
      integer  :: ip,jp
      ! Get the distances to the coordinate aligned planes
      if (abs(nx).lt.1e-6) then
         nx=0.0_WP
         Sx=0.0_WP
         ip=i+Sx
         dx=huge(1.0_WP)
      else
         Sx=sign(1.0_WP,nx)
         ip=i+Sx
         dx=(this%cfg%xm(ip)-x0)/nx
      end if
      if (abs(ny).lt.1e-6) then
         ny=0.0_WP
         Sy=0.0_WP
         jp=j+Sy
         dy=huge(1.0_WP)
      else
         Sy=sign(1.0_WP,ny)
         jp=j+Sy
         dy=(this%cfg%ym(jp)-y0)/ny
      end if
      ! Identify the closest plane
      if (abs(dx).lt.abs(dy)) then
         ! Interpolate at the intersection between the interface normal and the x plane
         if (ny.eq.0.0_WP) then
            xitp=this%cfg%xm(ip)
            yitp=y0
            Fitp=F(ip,j,k)
         else
            xitp=this%cfg%xm(ip)
            yitp=ny/nx*(xitp-x0)+y0
            w_m=abs(this%cfg%ym(jp)-yitp)
            w_p=abs(this%cfg%ym(j )-yitp)
            wsum=w_p+w_m
            w_p=w_p/(wsum)
            w_m=w_m/(wsum)
            Fitp=w_p*F(ip,jp,k)+w_m*F(ip,j,k)
         end if
      else
         ! Interpolate at the intersection between the interface normal and the y plane
         if (nx.eq.0.0_WP) then
            yitp=this%cfg%ym(jp)
            xitp=x0
            Fitp=F(i,jp,k)
         else
            ! Interpolate at the intersection
            yitp=this%cfg%ym(jp)
            xitp=nx/ny*(yitp-y0)+x0
            w_m=abs(this%cfg%xm(ip)-xitp)
            w_p=abs(this%cfg%xm(i )-xitp)
            wsum=w_p+w_m
            w_p=w_p/(wsum)
            w_m=w_m/(wsum)
            Fitp=w_p*F(ip,jp,k)+w_m*F(i,jp,k)
         end if
      end if
   end subroutine itp_ccplane


   subroutine get_one_sided_grad(this,phase,F,posI,i,j,k,normal,Fgrd)
      use messager, only: die
      implicit none
      class(lgpc), intent(in) :: this
      integer, intent(in) :: phase
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: F
      real(WP), dimension(3), intent(in) :: posI
      integer, intent(in) :: i,j,k
      real(WP), dimension(3), intent(in) :: normal
      real(WP), intent(out) :: Fgrd
      integer, dimension(3) :: ind
      real(WP) :: nx,ny,FI,F1,F2,mult,p
      real(WP) :: x1,y1,x2,y2,d1,d2,vof
      ! Get the phase multiplier (Follows IRL convention)
      p=real(phase,WP)     ! ( 0 for liquid, 1 for gas)
      mult=2.0_WP*p-1.0_WP ! (-1 for liquid, 1 for gas)
      ! Adjust the normal vecor
      nx=mult*normal(1)
      ny=mult*normal(2)
      ! Interpolate at the first plane
      x1=posI(1)
      y1=posI(2)
      vof=this%vf%VF(i,j,k)
      ind=[i,j,k]
      ! Get the cell indices containing the interpolated point (skip the interfacial cells)
      do while ((vof.gt.0.0_WP).and.(vof.lt.1.0_WP))
         call this%itp_ccplane(F=F,x0=x1,y0=y1,i=ind(1),j=ind(2),k=ind(3),nx=nx,ny=ny,Fitp=F1,xitp=x2,yitp=y2)
         x1=x2
         y1=y2
         ind=this%cfg%get_ijk_global(pos=[x1,y1,this%cfg%zm(1)],ind_guess=[i,j,k])
         vof=this%vf%VF(ind(1),ind(2),ind(3))
      end do
      ! Interpolate at the second plane
      call this%itp_ccplane(F=F,x0=x1,y0=y1,i=ind(1),j=ind(2),k=ind(3),nx=nx,ny=ny,Fitp=F2,xitp=x2,yitp=y2)
      ! Calculate the one-sided gradient
      FI=F(i,j,k)
      d1=norm2([posI(1)-x1,posI(2)-y1,0.0_WP])
      d2=norm2([posI(1)-x2,posI(2)-y2,0.0_WP])
      Fgrd=(p-this%vf%VF(i,j,k))*(F1-FI)/d1+(this%vf%VF(i,j,k)+mult-p)*(F2-FI)/d2
      ! if ((i.eq.154).and.(j.eq.129).and.(k.eq.1)) then
      !    print*,'----------------------------------------------------'
      !    print*,'i,j,k = ',i,j,k
      !    print*,'VOF = ',this%vf%VF(i,j,k)
      !    print*,'normal = ',nx,ny
      !    print*,'FI = ',FI
      !    print*,'----------------------------------------------------'
      !    print*,'x1 = ',x1
      !    print*,'y1 = ',y1
      !    print*,'F1 = ',F1
      !    print*,'d1 = ',d1
      !    print*,'grd1 = ',(F1-FI)/d1
      !    print*,'----------------------------------------------------'
      !    print*,'x2 = ',x2
      !    print*,'y2 = ',y2
      !    print*,'F2 = ',F2
      !    print*,'d2 = ',d2
      !    print*,'grd2 = ',(F2-FI)/d2
      !    print*,'----------------------------------------------------'
      !    print*,'Fgrd = ',Fgrd
      !    call die('')
      ! end if
   end subroutine get_one_sided_grad


   !> Filter an interfacial field in a conservative way
   subroutine filter(this,F,lvl,stc)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(lgpc), intent(in) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: F
      integer, intent(in) :: lvl,stc
      real(WP), dimension(:,:,:), allocatable   :: F_c
      real(WP), dimension(:,:,:), allocatable   :: w
      integer  :: i,j,k,ihat,jhat,khat,index,ierr,ilvl
      integer :: stxm,stym,stzm
      integer :: stxp,styp,stzp
      integer :: stx,sty,stz
      real(WP) :: my_F_int,F_int,F_int_org
      
      if (lvl.lt.1) return

      ! Check the dimensions
      if (this%nCell(1).gt.1) then
         stxm=-stc; stxp=+stc
      else
         stxm= 0; stxp= 0
      end if
      if (this%nCell(2).gt.1) then
         stym=-stc; styp=+stc
      else
         stym= 0; styp= 0
      end if
      if (this%nCell(3).gt.1) then
         stzm=-stc; stzp=+stc
      else
         stzm= 0; stzp= 0
      end if

      ! Allocate the arrays
      allocate(w(stxm:stxp,stym:styp,stzm:stzp))
      allocate(F_c(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Loop over the filtering levels
      do ilvl=1,lvl

         ! Initialize 
         F_c=F
         F=0.0_WP

         ! Loop over the interfacial cells
         do index=1,this%vf%band_count(0)
            ! Get the interfacial cell indices
            ihat=this%vf%band_map(1,index)
            jhat=this%vf%band_map(2,index)
            khat=this%vf%band_map(3,index)
            ! Initialize weights
            w=0.0_WP
            ! Loop over the stencil
            do stz=stzm,stzp
               k=khat+stz
               do sty=stym,styp
                  j=jhat+sty
                  do stx=stxm,stxp
                     i=ihat+stx
                     ! Calculate the weight
                     if (this%vf%VF(i,j,k).gt.VFlo.and.this%vf%VF(i,j,k).lt.VFhi) w(stx,sty,stz)=this%cfg%vol(i,j,k)
                  end do
               end do
            end do
            ! Normalize the weights
            w=w/sum(w)
            ! Loop over the stencil
            do stz=stzm,stzp
               k=khat+stz
               do sty=stym,styp
                  j=jhat+sty
                  do stx=stxm,stxp
                     i=ihat+stx
                     ! Apply the filter
                     if (w(stx,sty,stz).gt.0.0_WP) F(i,j,k)=F(i,j,k)+w(stx,sty,stz)*F_c(ihat,jhat,khat)
                  end do
               end do
            end do
         end do

         ! Sync the field
         call this%cfg%syncsum(F)
         call this%cfg%sync(F)

      end do

      ! ! Calculate the integral of the field
      ! my_F_int=0.0_WP
      ! do k=this%cfg%kmin_,this%cfg%kmax_
      !    do j=this%cfg%jmin_,this%cfg%jmax_
      !       do i=this%cfg%imin_,this%cfg%imax_
      !          my_F_int=my_F_int+this%cfg%vol(i,j,k)*F(i,j,k)
      !       end do
      !    end do
      ! end do
      ! call MPI_ALLREDUCE(my_F_int,F_int,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)

      ! my_F_int=0.0_WP
      ! do k=this%cfg%kmin_,this%cfg%kmax_
      !    do j=this%cfg%jmin_,this%cfg%jmax_
      !       do i=this%cfg%imin_,this%cfg%imax_
      !          my_F_int=my_F_int+this%cfg%vol(i,j,k)*F_c(i,j,k)
      !       end do
      !    end do
      ! end do
      ! call MPI_ALLREDUCE(my_F_int,F_int_org,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)

      ! if (this%cfg%amRoot) print*,'Integral error = ',F_int-F_int_org

      ! Deallocate arrays
      deallocate(w,F_c)

   end subroutine filter


   !> Filter a gas field in a conservative way
   subroutine filterG(this,F,lvl,stc)
      implicit none
      class(lgpc), intent(in) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: F
      integer, intent(in) :: lvl,stc
      real(WP), dimension(:,:,:), allocatable   :: F_c
      real(WP), dimension(:,:,:), allocatable   :: w
      integer :: i,j,k,ihat,jhat,khat,index,index_pure,ilvl
      integer :: stxm,stym,stzm
      integer :: stxp,styp,stzp
      integer :: stx,sty,stz
      
      if (lvl.lt.1) return

      ! Check the dimensions
      if (this%nCell(1).gt.1) then
         stxm=-stc; stxp=+stc
      else
         stxm= 0; stxp= 0
      end if
      if (this%nCell(2).gt.1) then
         stym=-stc; styp=+stc
      else
         stym= 0; styp= 0
      end if
      if (this%nCell(3).gt.1) then
         stzm=-stc; stzp=+stc
      else
         stzm= 0; stzp= 0
      end if

      ! Allocate the arrays
      allocate(w(stxm:stxp,stym:styp,stzm:stzp))
      allocate(F_c(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Loop over the filtering levels
      do ilvl=1,lvl

         ! Initialize 
         F_c=F
         F=0.0_WP

         ! Loop over the gas cells
         do khat=this%cfg%kmin_,this%cfg%kmax_
            do jhat=this%cfg%jmin_,this%cfg%jmax_
               do ihat=this%cfg%imin_,this%cfg%imax_
                  if (F_c(ihat,jhat,khat).eq.0.0_WP) cycle
                  ! Initialize weights
                  w=0.0_WP
                  ! Loop over the stencil
                  do stz=stzm,stzp
                     k=khat+stz
                     do sty=stym,styp
                        j=jhat+sty
                        do stx=stxm,stxp
                           i=ihat+stx
                           ! Calculate the weight
                           if (this%vf%VF(i,j,k).le.VFlo) w(stx,sty,stz)=this%cfg%vol(i,j,k)
                        end do
                     end do
                  end do
                  ! Normalize the weights
                  w=w/sum(w)
                  ! Loop over the stencil
                  do stz=stzm,stzp
                     k=khat+stz
                     do sty=stym,styp
                        j=jhat+sty
                        do stx=stxm,stxp
                           i=ihat+stx
                           ! Apply the filter
                           if (w(stx,sty,stz).gt.0.0_WP) F(i,j,k)=F(i,j,k)+w(stx,sty,stz)*F_c(ihat,jhat,khat)
                        end do
                     end do
                  end do
               end do
            end do
         end do

         ! Sync the field
         call this%cfg%syncsum(F)
         call this%cfg%sync(F)

      end do

      ! Deallocate arrays
      deallocate(w,F_c)

   end subroutine filterG


   !> Filter a liquid field in a conservative way
   subroutine filterL(this,F,lvl,stc)
      implicit none
      class(lgpc), intent(in) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: F
      integer, intent(in) :: lvl,stc
      real(WP), dimension(:,:,:), allocatable   :: F_c
      real(WP), dimension(:,:,:), allocatable   :: w
      integer :: i,j,k,ihat,jhat,khat,index,index_pure,ilvl
      integer :: stxm,stym,stzm
      integer :: stxp,styp,stzp
      integer :: stx,sty,stz
      
      if (lvl.lt.1) return

      ! Check the dimensions
      if (this%nCell(1).gt.1) then
         stxm=-stc; stxp=+stc
      else
         stxm= 0; stxp= 0
      end if
      if (this%nCell(2).gt.1) then
         stym=-stc; styp=+stc
      else
         stym= 0; styp= 0
      end if
      if (this%nCell(3).gt.1) then
         stzm=-stc; stzp=+stc
      else
         stzm= 0; stzp= 0
      end if

      ! Allocate the arrays
      allocate(w(stxm:stxp,stym:styp,stzm:stzp))
      allocate(F_c(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Loop over the filtering levels
      do ilvl=1,lvl

         ! Initialize 
         F_c=F
         F=0.0_WP

         ! Loop over the liquid cells
         do khat=this%cfg%kmin_,this%cfg%kmax_
            do jhat=this%cfg%jmin_,this%cfg%jmax_
               do ihat=this%cfg%imin_,this%cfg%imax_
                  if (F_c(ihat,jhat,khat).eq.0.0_WP) cycle
                  ! Initialize weights
                  w=0.0_WP
                  ! Loop over the stencil
                  do stz=stzm,stzp
                     k=khat+stz
                     do sty=stym,styp
                        j=jhat+sty
                        do stx=stxm,stxp
                           i=ihat+stx
                           ! Calculate the weight
                           if (this%vf%VF(i,j,k).ge.VFhi) w(stx,sty,stz)=this%cfg%vol(i,j,k)
                        end do
                     end do
                  end do
                  ! Normalize the weights
                  w=w/sum(w)
                  ! Loop over the stencil
                  do stz=stzm,stzp
                     k=khat+stz
                     do sty=stym,styp
                        j=jhat+sty
                        do stx=stxm,stxp
                           i=ihat+stx
                           ! Apply the filter
                           if (w(stx,sty,stz).gt.0.0_WP) F(i,j,k)=F(i,j,k)+w(stx,sty,stz)*F_c(ihat,jhat,khat)
                        end do
                     end do
                  end do
               end do
            end do
         end do

         ! Sync the field
         call this%cfg%syncsum(F)
         call this%cfg%sync(F)

      end do

      ! Deallocate arrays
      deallocate(w,F_c)

   end subroutine filterL


   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,face,dir)
      use string,         only: lowercase
      use messager,       only: die
      use iterator_class, only: locator_ftype
      implicit none
      class(lgpc), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer,  intent(in) :: type
      procedure(locator_ftype) :: locator
      character(len=1), intent(in) :: face
      integer, intent(in) :: dir
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type
      select case (new_bc%type)
         case (symmetry)
         case default
            call die('[lgpc add_bcond] Unknown bcond type')
      end select
      select case (lowercase(face))
         case ('x'); new_bc%face='x'
         case ('y'); new_bc%face='y'
         case ('z'); new_bc%face='z'
         case default; call die('[lgpc add_bcond] Unknown bcond face - expecting x, y, or z')
      end select
      select case (dir) ! Outward-oriented
         case (+1); new_bc%dir=+1
         case (-1); new_bc%dir=-1
         case ( 0); new_bc%dir= 0
         case default; call die('[lgpc add_bcond] Unknown bcond dir - expecting -1, +1, or 0')
      end select
      new_bc%itr=iterator(this%cfg,new_bc%name,locator,'c')
      
      ! Insert it up front
      new_bc%next=>this%first_bc
      this%first_bc=>new_bc
      
      ! Increment bcond counter
      this%nbc=this%nbc+1
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(lgpc), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[lgpc get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary conditions
   subroutine apply_bcond(this)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(lgpc), intent(inout) :: this
      integer :: i,j,k,dir,n
      integer :: ii,jj,kk
      integer :: i_out,j_out,k_out
      integer :: i_in,j_in,k_in
      type(bcond), pointer :: my_bc
      
      ! Traverse bcond list
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside the bcond work here
         if (my_bc%itr%amIn) then
            
            ! Select appropriate action based on the bcond type
            select case (my_bc%type)
               
               case (symmetry)             ! Apply symmetry conditions
                  
                  select case (my_bc%face)
                     case ('x')
                        do n=1,my_bc%itr%n_
                           i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                           do ii=1,extp_stc
                              i_out=i+my_bc%dir*(ii-1)
                              i_in =i-my_bc%dir*ii
                              this%vf%VF  (i_out,j,k)        = this%vf%VF  (i_in,j,k)
                              this%normal (i_out,j,k,1)      =-this%normal (i_in,j,k,1)
                              this%normal (i_out,j,k,2)      = this%normal (i_in,j,k,2)
                              this%normal (i_out,j,k,3)      = this%normal (i_in,j,k,3)
                              this%mdot3pLG(i_out,j,k,Lphase)= this%mdot3pLG(i_in,j,k,Lphase)
                              this%mdot3pLG(i_out,j,k,Gphase)= this%mdot3pLG(i_in,j,k,Gphase)
                              this%Tl     (i_out,j,k)        = this%Tl     (i_in,j,k)
                              this%Tg     (i_out,j,k)        = this%Tg     (i_in,j,k)
                           end do
                        end do
                     case ('y')
                        do n=1,my_bc%itr%n_
                           i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                           do jj=1,extp_stc
                              j_out=j+my_bc%dir*(jj-1)
                              j_in =j-my_bc%dir*jj
                              this%vf%VF  (i,j_out,k)        = this%vf%VF  (i,j_in,k)
                              this%normal (i,j_out,k,1)      = this%normal (i,j_in,k,1)
                              this%normal (i,j_out,k,2)      =-this%normal (i,j_in,k,2)
                              this%normal (i,j_out,k,3)      = this%normal (i,j_in,k,3)
                              this%mdot3pLG(i,j_out,k,Lphase)= this%mdot3pLG(i,j_in,k,Lphase)
                              this%mdot3pLG(i,j_out,k,Gphase)= this%mdot3pLG(i,j_in,k,Gphase)
                              this%Tl     (i,j_out,k)        = this%Tl     (i,j_in,k)
                              this%Tg     (i,j_out,k)        = this%Tg     (i,j_in,k)
                           end do
                        end do
                     case ('z')
                        do n=1,my_bc%itr%n_
                           i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                           do kk=1,extp_stc
                              k_out=k+my_bc%dir*(kk-1)
                              k_in =k-my_bc%dir*kk
                              this%vf%VF  (i,j,k_out)        = this%vf%VF  (i,j,k_in)
                              this%normal (i,j,k_out,1)      = this%normal (i,j,k_in,1)
                              this%normal (i,j,k_out,2)      = this%normal (i,j,k_in,2)
                              this%normal (i,j,k_out,3)      =-this%normal (i,j,k_in,3)
                              this%mdot3pLG(i,j,k_out,Lphase)= this%mdot3pLG(i,j,k_in,Lphase)
                              this%mdot3pLG(i,j,k_out,Gphase)= this%mdot3pLG(i,j,k_in,Gphase)
                              this%Tl     (i,j,k_out)        = this%Tl     (i,j,k_in)
                              this%Tg     (i,j,k_out)        = this%Tg     (i,j,k_in)
                           end do
                        end do
                  end select
                  
               case default
                  call die('[lgpc apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
      ! Sync full fields after all bcond
      call this%cfg%sync(this%vf%VF)
      do dir=1,3
         call this%cfg%sync(this%normal(:,:,:,dir))
      end do
      call this%cfg%sync(this%mdot3pLG(:,:,:,Lphase))
      call this%cfg%sync(this%mdot3pLG(:,:,:,Gphase))
      call this%cfg%sync(this%Tl)
      call this%cfg%sync(this%Tg)
      
   end subroutine apply_bcond


   !> Calculate the CFL
   subroutine get_cfl(this,dt,U,V,W,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(lgpc), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(out) :: cfl
      integer :: i,j,k,ierr
      real(WP) :: my_CFL
      
      ! Set the CFL to zero
      my_CFL=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_CFL=max(my_CFL,abs(U(i,j,k))*this%cfg%dxmi(i))
               my_CFL=max(my_CFL,abs(V(i,j,k))*this%cfg%dymi(j))
               my_CFL=max(my_CFL,abs(W(i,j,k))*this%cfg%dzmi(k))
            end do
         end do
      end do
      my_CFL=my_CFL*dt
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_CFL,cfl,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_cfl


end module lgpc_class
