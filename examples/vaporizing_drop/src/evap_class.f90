!> Evaporation class:
module evap_class
   use precision,         only: WP
   use string,            only: str_medium
   use config_class,      only: config
   use timetracker_class, only: timetracker
   use vfs_class,         only: vfs
   use tpscalar_class,    only: Lphase,Gphase
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: evap
   
   ! Index shift
   integer, dimension(1:3,1:3) :: ind_shift=reshape([1,0,0,0,1,0,0,0,1], shape(ind_shift))

   ! Default parameters for the evaporation class
   integer, parameter :: ext_lvl=5                                     !< The extension level for the interface normal (needs to be smaller than the VOF solver nband for the shift_mflux to work properly)

   type :: arr_ptr_4d
      real(WP), pointer :: arr(:,:,:,:)
   end type arr_ptr_4d

   !> Evaporation object definition
   type :: evap

      ! This is our config
      class(config), pointer :: cfg                                    !< This is the config the object is build for
      
      ! This is the name of the object
      character(len=str_medium) :: name='UNNAMED_EVAP'                 !< Object name (default=UNNAMED_EVAP)

      ! This is the corresponding VOF solver
      class(vfs), pointer :: vf

      ! Liquid and gas densities
      real(WP) :: rho_l,rho_g

      ! Evaporation mass flux data
      real(WP), dimension(:,:,:),   allocatable :: mdotdp              !< Evaporation mass flux (m dot double prime of dimension M/(T*L^2))
      real(WP), dimension(:,:,:),   allocatable :: mflux               !< Evaporation mass flux scaled by the surface density
      real(WP), dimension(:,:,:,:), allocatable :: mfluxLG,mfluxLG_old !< Liquid/Gas side shifted evaporation mass flux scaled by the surface density
      real(WP), dimension(:,:,:),   allocatable :: div_src             !< Evaporatin source term (div(U) = div_src)
      real(WP), dimension(:,:,:),   allocatable :: div_src_old         !< The old evaporatin source term

      ! Pseudo time over which the mflux is being shifted
      type(timetracker), public :: pseudo_time

      ! Phase-change and interfacial velocity
      real(WP), dimension(:,:,:,:), allocatable :: normal              !< Interface normal vector
      real(WP), dimension(:,:,:,:), allocatable :: pseudo_vel          !< Pseudo velocity for shifting mflux

      ! Metrics (point to flow solver metrics)
      type(arr_ptr_4d), dimension(:), allocatable :: itp               !< Cell to face interpolation coefficients
      type(arr_ptr_4d), dimension(:), allocatable :: div               !< Divergence operator (cell-centered)

      ! Number of cells in each direction
      integer, dimension(3) :: nCell
      
      ! Monitoring quantities
      real(WP) :: mflux_int,mflux_tol                                  !< Integral and tolerance of the scaled evap mass flux
      real(WP) :: mfluxL_int,mfluxL_err,mfluxL_int_err                 !< Liquid side scaled evap mass flux maximum, integral, and error
      real(WP) :: mfluxG_int,mfluxG_err,mfluxG_int_err                 !< Gas side scaled evap mass flux maximum, integral, and error
      
   contains

      procedure :: initialize                                          !< Class initializer
      procedure :: init_mfluxLG                                        !< Initialize evaporation mass flux on the liquid and gas sides
      procedure :: get_normal                                          !< Get the interface normal vector
      procedure :: get_pseudo_vel                                      !< Get the face-centered normilized gradient of VOF
      procedure :: get_mflux                                           !< Get the volumetric evaporation mass fllux
      procedure :: get_mflux                                           !< Get the volumetric evaporation mass fllux
      procedure :: get_dmfluxdt                                        !< Get the time derivative of the evaporation mass fllux
      procedure :: shift_mflux                                         !< Shift the evaporation mass flux
      procedure :: get_div                                             !< Get the evaporation source term
      procedure :: get_cfl                                             !< Get the CFL
      procedure :: extend_normal                                       !< Extend the interface normal                                             

   end type evap

   
contains
   
   
   !> Object initializer
   subroutine initialize(this,cfg,vf,itp_x,itp_y,itp_z,div_x,div_y,div_z,name)
      use messager,  only: die
      implicit none
      class(evap), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      class(vfs), target, intent(in) :: vf
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
      allocate(this%mdotdp     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));    this%mdotdp     =0.0_WP
      allocate(this%mflux      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));    this%mflux      =0.0_WP
      allocate(this%mfluxLG    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,0:1));this%mfluxLG    =0.0_WP
      allocate(this%mfluxLG_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,0:1));this%mfluxLG_old=0.0_WP
      allocate(this%div_src    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));    this%div_src    =0.0_WP
      allocate(this%div_src_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));    this%div_src_old=0.0_WP
      allocate(this%normal     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3));this%normal     =0.0_WP
      allocate(this%pseudo_vel (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3));this%pseudo_vel =0.0_WP
      allocate(this%itp(3))
      allocate(this%div(3))

      ! Set metrics
      this%itp(1)%arr=>itp_x
      this%itp(2)%arr=>itp_y
      this%itp(3)%arr=>itp_z
      this%div(1)%arr=>div_x
      this%div(2)%arr=>div_y
      this%div(3)%arr=>div_z

      ! Number of cells in each direction
      this%nCell(1)=this%cfg%nx
      this%nCell(2)=this%cfg%ny
      this%nCell(3)=this%cfg%nz

      ! Create a pseudo time
      this%pseudo_time=timetracker(amRoot=this%cfg%amRoot,name='Pseudo',print_info=.false.)
      this%pseudo_time=timetracker(amRoot=this%cfg%amRoot,name='Pseudo',print_info=.false.)

   end subroutine initialize


   !> Initialize the liquid and gas side volumetric mass fluxes
   subroutine init_mfluxLG(this)
      implicit none
      class(evap), intent(inout) :: this
      ! Initialize the liquid and gas side mfluxes
      this%mfluxLG(:,:,:,Lphase)=-this%mflux
      this%mfluxLG(:,:,:,Gphase)= this%mflux
      ! Initialize errors to zero
      this%mfluxL_err    =0.0_WP
      this%mfluxG_err    =0.0_WP
      this%mfluxL_int_err=0.0_WP
      this%mfluxG_int_err=0.0_WP
      ! Get the integrals
      call this%cfg%integrate(this%mflux,this%mflux_int)
      call this%cfg%integrate(this%mfluxLG(:,:,:,Lphase),this%mfluxL_int)
      call this%cfg%integrate(this%mfluxLG(:,:,:,Gphase),this%mfluxG_int)
   end subroutine init_mfluxLG
   

   !> Get the interfacial normal vector
   subroutine get_normal(this)
      use irl_fortran_interface, only: calculateNormal,getNumberOfVertices
      implicit none
      class(evap), intent(inout) :: this
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
      ! Check dimensions
      do dir=1,3
         if (this%nCell(dir).eq.1) this%normal(:,:,:,dir)=0.0_WP
      end do
   end subroutine get_normal


   !> Extend the interface normal for a smoother transition to zero
   subroutine extend_normal(this)
      implicit none
      class(evap), intent(inout) :: this
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
                     do sty=-1,1
                        do stx=-1,1
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
   end subroutine extend_normal


   !> Calculate the pseudo velocity used to shift the evaporation mass flux
   subroutine get_pseudo_vel(this)
      implicit none
      class(evap), intent(inout) :: this
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
   subroutine get_mflux(this)
      implicit none
      class(evap), intent(inout) :: this
      this%mflux=this%mdotdp*this%vf%SD
   end subroutine get_mflux

   
   !> Calculate the explicit mflux time derivative
   subroutine get_dmfluxdt(this,vel,mflux_old,dmfluxdt)
      implicit none
      class(evap), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(in)  :: vel        !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:3)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:),    intent(in)  :: mflux_old  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:),    intent(out) :: dmfluxdt   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,dir,index,n
      integer :: ic,jc,kc
      integer :: im,jm,km
      integer :: ip,jp,kp
      real(WP), dimension(:,:,:,:), allocatable :: F
      ! Zero out dmflux/dt array
      dmfluxdt=0.0_WP
      ! Allocate flux arrays
      allocate(F(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:3)); F=0.0_WP
      ! Fluxes of mflux
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
                  F(i,j,k,dir)=-0.5_WP*(vel(i,j,k,dir)+abs(vel(i,j,k,dir)))*mflux_old(im,jm,km) &
                  &            -0.5_WP*(vel(i,j,k,dir)-abs(vel(i,j,k,dir)))*mflux_old(ip,jp,kp)
               end if
            end do
         end do
      end do
      ! Time derivative of mflux
      do dir=1,3
         ! Loop over the cells
         do index=1,sum(this%vf%band_count(0:ext_lvl+1))
            ! Get the cell indices
            i=this%vf%band_map(1,index)
            j=this%vf%band_map(2,index)
            k=this%vf%band_map(3,index)
            ! Get the residual
            dmfluxdt(i,j,k)=dmfluxdt(i,j,k)                                                                            &
            &              +this%div(dir)%arr(0,i,j,k)*F(i,j,k,dir)                                                    &
            &              +this%div(dir)%arr(1,i,j,k)*F(i+ind_shift(1,dir),j+ind_shift(2,dir),k+ind_shift(3,dir),dir)
         end do
      end do
      ! Deallocate flux arrays
      deallocate(F)
      ! Sync residual
      call this%cfg%sync(dmfluxdt)
   end subroutine get_dmfluxdt


   !> Shift mflux away from the interface
   subroutine shift_mflux(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(evap), intent(inout) :: this
      real(WP), dimension(:,:,:), allocatable :: resmfluxL,resmfluxG
      integer  :: i,j,k,ierr
      real(WP) :: my_mflux_int,mflux_err

      ! Allocate memory for mflux residuals
      allocate(resmfluxL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(resmfluxG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Get the interface normal
      call this%get_normal()

      ! Extend the interface normal
      call this%extend_normal()

      ! Get the normalized gradient of VOF
      call this%get_pseudo_vel()

      ! Get the CFL based on the gradient of the VOF
      call this%get_cfl(this%pseudo_time%dt,this%pseudo_vel(:,:,:,1),this%pseudo_vel(:,:,:,2),this%pseudo_vel(:,:,:,3),this%pseudo_time%cfl)
      
      ! Reset the pseudo time
      call this%pseudo_time%reset()
      
      ! Adjust the pseudo time step
      call this%pseudo_time%adjust_dt()
      
      ! Initialize the evaporation mass fluxes on the liquid and gas sides
      this%mfluxLG(:,:,:,Lphase)=-this%mflux
      this%mfluxLG(:,:,:,Gphase)= this%mflux

      ! Move the evaporation mass flux away from the interface
      do while (.not.this%pseudo_time%done())
         
         ! Remember old mflux
         this%mfluxLG_old=this%mfluxLG
         
         ! Increment pseudo time
         call this%pseudo_time%increment()
         
         ! Assemble explicit residual
         call this%get_dmfluxdt(vel= this%pseudo_vel,mflux_old=this%mfluxLG_old(:,:,:,Lphase),dmfluxdt=resmfluxL)
         call this%get_dmfluxdt(vel=-this%pseudo_vel,mflux_old=this%mfluxLG_old(:,:,:,Gphase),dmfluxdt=resmfluxG)
         
         ! Apply these residuals
         this%mfluxLG(:,:,:,Lphase)=this%mfluxLG_old(:,:,:,Lphase)+this%pseudo_time%dt*resmfluxL
         this%mfluxLG(:,:,:,Gphase)=this%mfluxLG_old(:,:,:,Gphase)+this%pseudo_time%dt*resmfluxG

         ! Calculate the integral of the residual error of mfluxL
         my_mflux_int=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%vf%VF(i,j,k).lt.1) my_mflux_int=my_mflux_int+this%cfg%vol(i,j,k)*this%cfg%VF(i,j,k)*this%vf%VF(i,j,k)*abs(this%mfluxLG(i,j,k,Lphase))
               end do
            end do
         end do
         call MPI_ALLREDUCE(my_mflux_int,this%mfluxL_err,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         
         ! Calculate the integral of the residual error of mfluxG
         my_mflux_int=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%vf%VF(i,j,k).gt.0) my_mflux_int=my_mflux_int+this%cfg%vol(i,j,k)*this%cfg%VF(i,j,k)*(1.0_WP-this%vf%VF(i,j,k))*abs(this%mfluxLG(i,j,k,Gphase))
               end do
            end do
         end do
         call MPI_ALLREDUCE(my_mflux_int,this%mfluxG_err,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)

         ! Check convergence
         mflux_err=max(this%mfluxL_err,this%mfluxG_err)
         if (mflux_err.lt.this%mflux_tol) exit

      end do

      ! Integral of mflux
      call this%cfg%integrate(this%mflux,this%mflux_int)
      call this%cfg%integrate(this%mfluxLG(:,:,:,Lphase),this%mfluxL_int)
      call this%cfg%integrate(this%mfluxLG(:,:,:,Gphase),this%mfluxG_int)
      this%mfluxL_int_err=abs(abs(this%mfluxL_int)-this%mflux_int)
      this%mfluxG_int_err=abs(abs(this%mfluxG_int)-this%mflux_int)

      ! Deallocate mflux residuals
      deallocate(resmfluxL,resmfluxG)

   end subroutine shift_mflux


   ! !> Shift mflux away from the interface (Boyd and Ling) serial
   ! subroutine shift_mflux(this)
   !    use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
   !    use parallel,  only: MPI_REAL_WP
   !    implicit none
   !    class(evap), intent(inout) :: this
   !    real(WP), dimension(:,:,:), allocatable :: resmfluxL,resmfluxG
   !    integer  :: i,j,k,ihat,jhat,khat,index,ierr
   !    real(WP) :: my_mflux_int,mflux_err,w
   !    real(WP), dimension(:,:,:), allocatable :: wG,wL
   !    integer :: stxm,stym,stzm
   !    integer :: stxp,styp,stzp
   !    integer :: stx,sty,stz
   !    real(WP), dimension(3) :: d
      
   !    ! Check the dimensions
   !    if (this%nCell(1).gt.1) then
   !       stxm=-2; stxp=+2
   !    else
   !       stxm= 0; stxp= 0
   !    end if
   !    if (this%nCell(2).gt.1) then
   !       stym=-2; styp=+2
   !    else
   !       stym= 0; styp= 0
   !    end if
   !    if (this%nCell(3).gt.1) then
   !       stzm=-2; stzp=+2
   !    else
   !       stzm= 0; stzp= 0
   !    end if

   !    ! Allocate the weight matrix
   !    allocate(wG(stxm:stxp,stym:styp,stzm:stzp))
   !    allocate(wL(stxm:stxp,stym:styp,stzm:stzp))

   !    ! Get the interface normal
   !    call this%get_normal()

   !    ! Initialize one sided mfluxes
   !    this%mfluxG=0.0_WP
   !    this%mfluxL=0.0_WP

   !    ! Loop over the interfacial cells
   !    do index=1,this%vf%band_count(0)
   !       ! Get the interfacial cell indices
   !       ihat=this%vf%band_map(1,index)
   !       jhat=this%vf%band_map(2,index)
   !       khat=this%vf%band_map(3,index)
   !       ! Initialize weights
   !       wG=0.0_WP
   !       wL=0.0_WP
   !       ! Loop over the stencil
   !       do stz=stzm,stzp
   !          k=khat+stz
   !          do sty=stym,styp
   !             j=jhat+sty
   !             do stx=stxm,stxp
   !                i=ihat+stx
   !                ! Calculate the weight
   !                d=[this%cfg%xm(ihat)-this%cfg%xm(i),this%cfg%ym(jhat)-this%cfg%ym(j),this%cfg%zm(khat)-this%cfg%zm(k)]
   !                w=abs(sum(d*this%normal(ihat,jhat,khat,:)))/norm2(d)
   !                ! Assign the weight to the pure cell
   !                if (this%vf%VF(i,j,k).eq.0.0_WP) then
   !                   wG(stx,sty,stz)=w
   !                else if (this%vf%VF(i,j,k).eq.1.0_WP) then
   !                   wL(stx,sty,stz)=w
   !                end if
   !             end do
   !          end do
   !       end do
   !       ! Normalize the weights
   !       wG=wG/(sum(wG)+tiny(1.0_WP))
   !       wL=wL/(sum(wL)+tiny(1.0_WP))
   !       ! Loop over the stencil
   !       do stz=stzm,stzp
   !          k=khat+stz
   !          do sty=stym,styp
   !             j=jhat+sty
   !             do stx=stxm,stxp
   !                i=ihat+stx
   !                ! Update the mflux of the pure cell
   !                if (this%vf%VF(i,j,k).eq.0.0_WP) then
   !                   this%mfluxG(i,j,k)=this%mfluxG(i,j,k)+wG(stx,sty,stz)*this%mflux(ihat,jhat,khat)
   !                else if (this%vf%VF(i,j,k).eq.1.0_WP) then
   !                   this%mfluxL(i,j,k)=this%mfluxL(i,j,k)+wL(stx,sty,stz)*this%mflux(ihat,jhat,khat)
   !                end if
   !             end do
   !          end do
   !       end do
   !    end do

   !    ! Sync
   !    ! call this%cfg%sync(this%mfluxL)
   !    ! call this%cfg%sync(this%mfluxG)

   !    ! Integral of mflux
   !    call this%cfg%integrate(this%mflux,this%mflux_int)
   !    call this%cfg%integrate(this%mfluxL,this%mfluxL_int)
   !    call this%cfg%integrate(this%mfluxG,this%mfluxG_int)
   !    this%mfluxL_int_err=abs(this%mfluxL_int-this%mflux_int)
   !    this%mfluxG_int_err=abs(this%mfluxG_int-this%mflux_int)

   !    ! Deallocate
   !    deallocate(wG,wL)

   ! end subroutine shift_mflux


   ! !> Shift mflux away from the interface (Boyd and Ling) parallel
   ! subroutine shift_mflux(this)
   !    use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
   !    use parallel,  only: MPI_REAL_WP
   !    use vfs_class, only: VFlo,VFhi
   !    implicit none
   !    class(evap), intent(inout) :: this
   !    real(WP), dimension(:,:,:), allocatable :: resmfluxL,resmfluxG
   !    integer  :: i,j,k,ihat,jhat,khat,index,ierr
   !    real(WP) :: my_mflux_int,mflux_err,w
   !    real(WP), dimension(:,:,:), allocatable :: wG,wL
   !    integer :: stxm,stym,stzm
   !    integer :: stxp,styp,stzp
   !    integer :: stx,sty,stz
   !    real(WP), dimension(3) :: d
      
   !    ! Check the dimensions
   !    if (this%nCell(1).gt.1) then
   !       stxm=-2; stxp=+2
   !    else
   !       stxm= 0; stxp= 0
   !    end if
   !    if (this%nCell(2).gt.1) then
   !       stym=-2; styp=+2
   !    else
   !       stym= 0; styp= 0
   !    end if
   !    if (this%nCell(3).gt.1) then
   !       stzm=-2; stzp=+2
   !    else
   !       stzm= 0; stzp= 0
   !    end if

   !    ! Allocate the weight matrix
   !    ! allocate(wG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
   !    ! allocate(wL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
   !    allocate(wG(stxm:stxp,stym:styp,stzm:stzp))
   !    allocate(wL(stxm:stxp,stym:styp,stzm:stzp))

   !    ! Get the interface normal
   !    call this%get_normal()

   !    ! Initialize one sided mfluxes
   !    this%mfluxG=0.0_WP
   !    this%mfluxL=0.0_WP

   !    ! Loop over the interfacial cells
   !    do khat=this%cfg%kmin_+stzm,this%cfg%kmax_+stzp
   !       do jhat=this%cfg%jmin_+stym,this%cfg%jmax_+styp
   !          do ihat=this%cfg%imin_+stxm,this%cfg%imax_+stxp
   !             if (this%vf%VF(ihat,khat,jhat).gt.VFlo.and.this%vf%VF(ihat,khat,jhat).lt.VFhi) then
   !                ! Initialize weights
   !                wG=0.0_WP
   !                wL=0.0_WP
   !                ! Loop over the stencil
   !                ! do k=khat+stzm,khat+stzp
   !                !    do j=jhat+stym,jhat+styp
   !                !       do i=ihat+stxm,ihat+stxp
   !                do stz=stzm,stzp
   !                   k=khat+stz
   !                   do sty=stym,styp
   !                      j=jhat+sty
   !                      do stx=stxm,stxp
   !                         i=ihat+stx
   !                         ! Calculate the weight
   !                         d=[this%cfg%xm(ihat)-this%cfg%xm(i),this%cfg%ym(jhat)-this%cfg%ym(j),this%cfg%zm(khat)-this%cfg%zm(k)]
   !                         w=abs(sum(d*this%normal(ihat,jhat,khat,:)))/norm2(d)
   !                         ! w=1.0_WP
   !                         ! Assign the weight to the pure cell
   !                         if (this%vf%VF(i,j,k).eq.0.0_WP) then
   !                            ! wG(i,j,k)=w
   !                            wG(stx,sty,stz)=w
   !                         else if (this%vf%VF(i,j,k).eq.1.0_WP) then
   !                            ! wL(i,j,k)=w
   !                            wL(stx,sty,stz)=w
   !                         end if
   !                      end do
   !                   end do
   !                end do
   !                ! Normalize the weights
   !                wG=wG/(sum(wG)+tiny(1.0_WP))
   !                wL=wL/(sum(wL)+tiny(1.0_WP))
   !                ! Loop over the stencil
   !                ! do k=khat+stzm,khat+stzp
   !                !    do j=jhat+stym,jhat+styp
   !                !       do i=ihat+stxm,ihat+stxp
   !                do stz=stzm,stzp
   !                   k=khat+stz
   !                   do sty=stym,styp
   !                      j=jhat+sty
   !                      do stx=stxm,stxp
   !                         i=ihat+stx
   !                         ! Update the mflux of the pure cell
   !                         if (this%vf%VF(i,j,k).eq.0.0_WP) then
   !                            ! this%mfluxG(i,j,k)=this%mfluxG(i,j,k)+wG(i,j,k)*this%mflux(ihat,jhat,khat)
   !                            this%mfluxG(i,j,k)=this%mfluxG(i,j,k)+wG(stx,sty,stz)*this%mflux(ihat,jhat,khat)
   !                         else if (this%vf%VF(i,j,k).eq.1.0_WP) then
   !                            ! this%mfluxL(i,j,k)=this%mfluxL(i,j,k)+wL(i,j,k)*this%mflux(ihat,jhat,khat)
   !                            this%mfluxL(i,j,k)=this%mfluxL(i,j,k)+wL(stx,sty,stz)*this%mflux(ihat,jhat,khat)
   !                         end if
   !                      end do
   !                   end do
   !                end do
   !             end if
   !          end do
   !       end do
   !    end do

   !    ! Sync
   !    call this%cfg%sync(this%mfluxL)
   !    call this%cfg%sync(this%mfluxG)

   !    ! Integral of mflux
   !    call this%cfg%integrate(this%mflux,this%mflux_int)
   !    call this%cfg%integrate(this%mfluxL,this%mfluxL_int)
   !    call this%cfg%integrate(this%mfluxG,this%mfluxG_int)
   !    this%mfluxL_int_err=abs(this%mfluxL_int-this%mflux_int)
   !    this%mfluxG_int_err=abs(this%mfluxG_int-this%mflux_int)

   !    ! Deallocate
   !    deallocate(wG,wL)

   ! end subroutine shift_mflux


   !> Calculate the divergence induced by phase change
   subroutine get_div(this)
      implicit none
      class(evap), intent(inout) :: this
      this%div_src=this%mfluxLG(:,:,:,Gphase)/this%rho_g+this%mfluxLG(:,:,:,Lphase)/this%rho_l
   end subroutine get_div


   !> Calculate the CFL
   subroutine get_cfl(this,dt,U,V,W,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(evap), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
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


   !> Calculate the maximum of the phase-change velocity
   subroutine get_max_vel_pc(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(evap), intent(inout) :: this
      real(WP) :: my_Upcmax,my_Vpcmax,my_Wpcmax
      integer  :: i,j,k,ierr
      ! Set all to zero
      my_Upcmax=0.0_WP; my_Vpcmax=0.0_WP; my_Wpcmax=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_Upcmax=max(my_Upcmax,abs(this%vel_pc(i,j,k,1)))
               my_Vpcmax=max(my_Vpcmax,abs(this%vel_pc(i,j,k,2)))
               my_Wpcmax=max(my_Wpcmax,abs(this%vel_pc(i,j,k,3)))
            end do
         end do
      end do
      ! Get the parallel max
      call MPI_ALLREDUCE(my_Upcmax,this%Upcmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Vpcmax,this%Vpcmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Wpcmax,this%Wpcmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
   end subroutine get_max_vel_pc


end module evap_class
