!> Multiphase compressible flow solver class:
!> Provides support for RHS calculation only
module mpcomp_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use irl_fortran_interface
   implicit none
   private
   
   ! Expose type
   public :: mpcomp,VFlo,VFhi
   
   ! Default parameters for volume fraction solver
   integer,  parameter :: nband=1                                    !< Number of cells around the interfacial cells on which localized work is performed
   integer,  parameter :: advect_band=1                              !< How far we do SL transport
   integer,  parameter :: half_edge=1                                !< Half-edge cutting (default)
   real(WP), parameter :: VFlo=1.0e-12_WP                            !< Minimum VF value considered
   real(WP), parameter :: VFhi=1.0_WP-VFlo                           !< Maximum VF value considered
   real(WP), parameter :: volume_epsilon_factor =1.0e-15_WP          !< Minimum volume  to consider for computational geometry (normalized by min_meshsize**3)
   real(WP), parameter :: surface_epsilon_factor=1.0e-15_WP          !< Minimum surface to consider for computational geometry (normalized by min_meshsize**2)
   real(WP), parameter :: iterative_distfind_tol=1.0e-12_WP          !< Tolerance for iterative plane distance finding
   
   !> Multiphase compressible solver object definition
   type :: mpcomp
      
      ! This is the config around which solver is built
      class(config), pointer :: cfg
      
      ! Solver name
      character(len=str_medium) :: name='UNNAMED_MPCOMP'
      
      ! Gravitational acceleration
      real(WP), dimension(3) :: gravity=0.0_WP
      
      ! Volume moments
      real(WP), dimension(:,:,:)  , allocatable :: VF,VFold
      real(WP), dimension(:,:,:,:), allocatable :: Lbary,Gbary
      
      ! IRL-native data
      type(ByteBuffer_type) :: send_byte_buffer
      type(ByteBuffer_type) :: recv_byte_buffer
      type(ObjServer_PlanarSep_type)  :: planar_separator_allocation
      type(ObjServer_PlanarLoc_type)  :: planar_localizer_allocation
      type(ObjServer_LocSepLink_type) :: localized_separator_link_allocation
      type(ObjServer_LocLink_type)    :: localizer_link_allocation
      type(PlanarLoc_type),   dimension(:,:,:),   allocatable :: localizer
      type(PlanarSep_type),   dimension(:,:,:),   allocatable :: liquid_gas_interface
      type(LocSepLink_type),  dimension(:,:,:),   allocatable :: localized_separator_link
      type(LocLink_type),     dimension(:,:,:),   allocatable :: localizer_link
      type(Poly_type),        dimension(:,:,:),   allocatable :: interface_polygon
      type(SepVM_type),       dimension(:,:,:,:), allocatable :: face_flux
      type(TagAccVM_SepVM_type), dimension(:,:,:,:), allocatable :: detailed_face_flux    !< Cell-decomposed face flux geometric data
      real(WP) :: dt_detailed_face_flux

      ! Conserved variables
      integer :: nVAR
      real(WP), dimension(:,:,:,:), allocatable :: VAR,VARold
      
      ! Flow variables
      real(WP), dimension(:,:,:), allocatable :: U,V,W
      real(WP), dimension(:,:,:), allocatable :: RHOl
      real(WP), dimension(:,:,:), allocatable :: T
      
      ! Viscosity, diffusivity, and specific heat
      real(WP), dimension(:,:,:), allocatable :: visc,beta,diff
      real(WP) :: Cv=1.0_WP
      
      ! Store mesh size
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      
      ! CFL numbers
      real(WP) :: CFLc_x,CFLc_y,CFLc_z                    !< Convective CFL numbers
      real(WP) :: CFLa_x,CFLa_y,CFLa_z                    !< Acoustic CFL numbers
      real(WP) :: CFLv_x,CFLv_y,CFLv_z                    !< Viscous CFL numbers
      
      ! Monitoring quantities
      real(WP) :: VFmin,VFmax,VFint                       !< Volume fraction stats
      real(WP) :: Umax,Vmax,Wmax                          !< Velocity stats
      real(WP) :: RHOlmin,RHOlmax,RHOlint                 !< Liquid density stats
      !real(WP) :: RHOUint,RHOVint,RHOWint                 !< Momentum stats
      !real(WP) :: RHOmin,RHOmax,RHOint                    !< Density stats
      !real(WP) :: Emin,Emax,RHOEint                       !< Internal energy stats
      !real(WP) :: Pmin,Pmax                               !< Pressure stats
      
   contains
      procedure :: print=>mpcomp_print                    !< Output solver to the screen
      procedure :: initialize                             !< Initialize the flow solver
      procedure :: initialize_irl                         !< Initialize interface with interface reconstruction library
      procedure :: advance_volume                         !< Advance volume equation
      procedure :: build_interface                        !< Build interface from volume moments
      procedure :: get_RHS                                !< Calculate dQ/dt
      procedure :: get_velocity                           !< Calculate velocity from momentum
      procedure :: get_momentum                           !< Calculate momentum from velocity
      procedure :: interp_vel                             !< Calculate interpolated velocity
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: get_info                               !< Calculate maximum field values
      procedure :: update_surfmesh                        !< Update a surfmesh object using current polygons
   end type mpcomp
   
   
contains
   
   
   !> Initialization for compressible flow solver
   subroutine initialize(this,cfg,name)
      use messager, only: die
      implicit none
      class(mpcomp) :: this
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to config object
      this%cfg=>cfg
      
      ! Check that config is uniform with at least 3 cells of overlap
      if (this%cfg%no.lt.3) call die('[mpcomp initialize] mpcomp solver requires at least 3 cells of overlap')
      if (.not.all([this%cfg%uniform_x,this%cfg%uniform_y,this%cfg%uniform_z])) call die('[mpcomp initialize] mpcomp solver requires a uniform mesh')
      
      ! Store constant cell size and its inverse
      this%dx=this%cfg%dx(this%cfg%imin_); this%dxi=1.0_WP/this%dx
      this%dy=this%cfg%dy(this%cfg%jmin_); this%dyi=1.0_WP/this%dy
      this%dz=this%cfg%dz(this%cfg%kmin_); this%dzi=1.0_WP/this%dz
      
      ! Allocate volume moments
      allocate(this%VF   (  this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VF   =0.0_WP
      allocate(this%VFold(  this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VFold=0.0_WP
      allocate(this%Lbary(3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Lbary=0.0_WP
      allocate(this%Gbary(3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Gbary=0.0_WP
      
      ! Initialize Interface Reconstruction Library
      call this%initialize_irl()
      
      ! Allocate and zero out conserved variables storage
      this%nVAR=1
      allocate(this%VAR   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nVAR)); this%VAR   =0.0_WP
      allocate(this%VARold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nVAR)); this%VARold=0.0_WP
      
      ! Flow variables
      allocate(this%U(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U=0.0_WP
      allocate(this%V(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%V=0.0_WP
      allocate(this%W(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%W=0.0_WP
      allocate(this%RHOl(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOl=0.0_WP
      
      
      !allocate(this%E(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%E=0.0_WP
      !allocate(this%P(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%P=0.0_WP
      !allocate(this%T(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%T=0.0_WP
      !allocate(this%visc(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc=0.0_WP
      !allocate(this%beta(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%beta=0.0_WP
      !allocate(this%diff(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%diff=0.0_WP
      
   end subroutine initialize
   
   
   !> Initialize IRL interface
   subroutine initialize_irl(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k,n,tag
      real(WP), dimension(3,4) :: vert
      integer(IRL_LargeOffsetIndex_t) :: total_cells
      
      ! Transfer small constants to IRL
      call setVFBounds(VFlo)
      call setVFTolerance_IterativeDistanceFinding(iterative_distfind_tol)
      call setMinimumVolToTrack(volume_epsilon_factor*this%cfg%min_meshsize**3)
      call setMinimumSAToTrack(surface_epsilon_factor*this%cfg%min_meshsize**2)
      
      ! Set IRL's moment calculation method
      call getMoments_setMethod(half_edge)
      
      ! Allocate IRL arrays
      allocate(this%localizer               (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%liquid_gas_interface    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%localized_separator_link(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%localizer_link          (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%interface_polygon       (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Arrays for storing face fluxes and detailed flux geometry
      allocate(this%face_flux         (1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%detailed_face_flux(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               do n=1,3
                  call new(this%face_flux(n,i,j,k))
                  call new(this%detailed_face_flux(n,i,j,k))
               end do
            end do
         end do
      end do
      
      ! Initialize size for IRL
      total_cells=int(this%cfg%nxo_,8)*int(this%cfg%nyo_,8)*int(this%cfg%nzo_,8)
      call new(this%planar_localizer_allocation,total_cells)
      call new(this%planar_separator_allocation,total_cells)
      call new(this%localized_separator_link_allocation,total_cells)
      call new(this%localizer_link_allocation,total_cells)
      
      ! Initialize arrays and setup linking
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Transfer cell to IRL
               call new(this%localizer(i,j,k),this%planar_localizer_allocation)
               call setFromRectangularCuboid(this%localizer(i,j,k),[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
               ! PLIC interface(s)
               call new(this%liquid_gas_interface(i,j,k),this%planar_separator_allocation)
               ! PLIC+mesh with connectivity (i.e., link)
               call new(this%localized_separator_link(i,j,k),this%localized_separator_link_allocation,this%localizer(i,j,k),this%liquid_gas_interface(i,j,k))
               ! Mesh with connectivity
               call new(this%localizer_link(i,j,k),this%localizer_link_allocation,this%localizer(i,j,k))
            end do
         end do
      end do
      
      ! Polygonal representation of the surface
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
                  call new(this%interface_polygon(i,j,k))
            end do
         end do
      end do
      
      ! Give each link a unique lexicographic tag (per processor)
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               tag=this%cfg%get_lexico_from_ijk([i,j,k])
               call setId(this%localized_separator_link(i,j,k),tag)
               call setId(this%localizer_link(i,j,k),tag)
            end do
         end do
      end do
      
      ! Set the connectivity
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! In the x- direction
               if (i.gt.this%cfg%imino_) then
                  call setEdgeConnectivity(this%localized_separator_link(i,j,k),0,this%localized_separator_link(i-1,j,k))
                  call setEdgeConnectivity(this%localizer_link(i,j,k),0,this%localizer_link(i-1,j,k))
               end if
               ! In the x+ direction
               if (i.lt.this%cfg%imaxo_) then
                  call setEdgeConnectivity(this%localized_separator_link(i,j,k),1,this%localized_separator_link(i+1,j,k))
                  call setEdgeConnectivity(this%localizer_link(i,j,k),1,this%localizer_link(i+1,j,k))
               end if
               ! In the y- direction
               if (j.gt.this%cfg%jmino_) then
                  call setEdgeConnectivity(this%localized_separator_link(i,j,k),2,this%localized_separator_link(i,j-1,k))
                  call setEdgeConnectivity(this%localizer_link(i,j,k),2,this%localizer_link(i,j-1,k))
               end if
               ! In the y+ direction
               if (j.lt.this%cfg%jmaxo_) then
                  call setEdgeConnectivity(this%localized_separator_link(i,j,k),3,this%localized_separator_link(i,j+1,k))
                  call setEdgeConnectivity(this%localizer_link(i,j,k),3,this%localizer_link(i,j+1,k))
               end if
               ! In the z- direction
               if (k.gt.this%cfg%kmino_) then
                  call setEdgeConnectivity(this%localized_separator_link(i,j,k),4,this%localized_separator_link(i,j,k-1))
                  call setEdgeConnectivity(this%localizer_link(i,j,k),4,this%localizer_link(i,j,k-1))
               end if
               ! In the z+ direction
               if (k.lt.this%cfg%kmaxo_) then
                  call setEdgeConnectivity(this%localized_separator_link(i,j,k),5,this%localized_separator_link(i,j,k+1))
                  call setEdgeConnectivity(this%localizer_link(i,j,k),5,this%localizer_link(i,j,k+1))
               end if
            end do
         end do
      end do
      
      ! Prepare byte storage for synchronization
      call new(this%send_byte_buffer)
      call new(this%recv_byte_buffer)
      
   end subroutine initialize_irl
   
   
   !> Advance volume equation using passed U/V/W and current plic
   subroutine advance_volume(this,dt,U,V,W)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt              !< Time step by which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer, dimension(:,:,:), allocatable :: band    !< Band to localize workload around the interface
      integer :: i,j,k
      
      ! Allocate the band
      allocate(band(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! First prepare a band
      prepare_band: block
         integer :: ii,jj,kk,dir,n
         integer, dimension(3) :: ind
         integer :: ibmin_,ibmax_,jbmin_,jbmax_,kbmin_,kbmax_
         ! Loop extents
         ibmin_=this%cfg%imin_; if (this%cfg%iproc.eq.1           .and..not.this%cfg%xper) ibmin_=this%cfg%imin-1
         ibmax_=this%cfg%imax_; if (this%cfg%iproc.eq.this%cfg%npx.and..not.this%cfg%xper) ibmax_=this%cfg%imax+1
         jbmin_=this%cfg%jmin_; if (this%cfg%jproc.eq.1           .and..not.this%cfg%yper) jbmin_=this%cfg%jmin-1
         jbmax_=this%cfg%jmax_; if (this%cfg%jproc.eq.this%cfg%npy.and..not.this%cfg%yper) jbmax_=this%cfg%jmax+1
         kbmin_=this%cfg%kmin_; if (this%cfg%kproc.eq.1           .and..not.this%cfg%zper) kbmin_=this%cfg%kmin-1
         kbmax_=this%cfg%kmax_; if (this%cfg%kproc.eq.this%cfg%npz.and..not.this%cfg%zper) kbmax_=this%cfg%kmax+1
         ! Initialize band
         band=(nband+1)*int(sign(1.0_WP,this%VF-0.5_WP))
         ! First sweep to identify cells with interface
         do k=kbmin_,kbmax_; do j=jbmin_,jbmax_; do i=ibmin_,ibmax_
            ! Identify cells with interface
            if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).le.VFhi) band(i,j,k)=0
            ! Check if cell-face is an interface
            do dir=1,3; do n=-1,+1,2
               ind=[i,j,k]; ind(dir)=ind(dir)+n
               if (this%VF(i,j,k).lt.VFlo.and.this%VF(ind(1),ind(2),ind(3)).gt.VFhi.or.&
               &   this%VF(i,j,k).gt.VFhi.and.this%VF(ind(1),ind(2),ind(3)).lt.VFlo) band(i,j,k)=0
            end do; end do
         end do; end do; end do
         call this%cfg%sync(band)
         ! Sweep to identify the bands up to nband
         do n=1,nband
            ! For each band
            do k=kbmin_,kbmax_; do j=jbmin_,jbmax_; do i=ibmin_,ibmax_
               ! Work on one band at a time
               if (abs(band(i,j,k)).gt.n) then
                  ! Loop over 26 neighbors and extend the band
                  do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                     if (abs(band(ii,jj,kk)).eq.n-1) band(i,j,k)=int(sign(real(n,WP),this%VF(i,j,k)-0.5_WP))
                  end do; end do; end do
               end if
            end do; end do; end do
            call this%cfg%sync(band)
         end do
      end block prepare_band
      
      ! Reset face fluxes to crude estimate (just needs to be valid for volume away from interface)
      reset_face_fluxes: block
         ! Loop through all cell faces
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            if (band(i,j,k).lt.0) then
               call construct(this%face_flux(1,i,j,k),[dt*U(i,j,k)*this%dy*this%dz,[0.0_WP,0.0_WP,0.0_WP],0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
               call construct(this%face_flux(2,i,j,k),[dt*V(i,j,k)*this%dz*this%dx,[0.0_WP,0.0_WP,0.0_WP],0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
               call construct(this%face_flux(3,i,j,k),[dt*W(i,j,k)*this%dx*this%dy,[0.0_WP,0.0_WP,0.0_WP],0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
            else
               call construct(this%face_flux(1,i,j,k),[0.0_WP,[0.0_WP,0.0_WP,0.0_WP],dt*U(i,j,k)*this%dy*this%dz,0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
               call construct(this%face_flux(2,i,j,k),[0.0_WP,[0.0_WP,0.0_WP,0.0_WP],dt*V(i,j,k)*this%dz*this%dx,0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
               call construct(this%face_flux(3,i,j,k),[0.0_WP,[0.0_WP,0.0_WP,0.0_WP],dt*W(i,j,k)*this%dx*this%dy,0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
            end if
            call clear(this%detailed_face_flux(1,i,j,k))
            call clear(this%detailed_face_flux(2,i,j,k))
            call clear(this%detailed_face_flux(3,i,j,k))
         end do; end do; end do
         ! Remember dt for which the fluxes are calculated
         this%dt_detailed_face_flux=dt
      end block reset_face_fluxes
      
      ! Loop over the domain and compute fluxes using semi-Lagrangian algorithm
      calculate_fluxes: block
         real(IRL_double), dimension(3,9) :: face
         type(CapDod_type) :: flux_polyhedron
         type(SepVM_type) :: my_SepVM
         real(WP), dimension(3) :: lbar,gbar
         real(WP) :: lvol,gvol
         integer :: n
         ! Allocate flux polyhedron
         call new(flux_polyhedron)
         ! Loop through all cell faces - extra cell on the left because mass flux will be needed by momentum later
         do k=this%cfg%kmin_-1,this%cfg%kmax_+1; do j=this%cfg%jmin_-1,this%cfg%jmax_+1; do i=this%cfg%imin_-1,this%cfg%imax_+1
            ! X flux
            if (minval(abs(band(i-1:i,j,k))).le.advect_band) then
               ! Construct and project face
               face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=back_project(face(:,1))
               face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=back_project(face(:,2))
               face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=back_project(face(:,3))
               face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=back_project(face(:,4))
               face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]; face(:,9)=back_project(face(:,9))
               ! Form flux polyhedron
               call construct(flux_polyhedron,face)
               ! Add solenoidal correction
               call adjustCapToMatchVolume(flux_polyhedron,dt*U(i,j,k)*this%dy*this%dz)
               ! Need full geometric flux
               call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%detailed_face_flux(1,i,j,k))
               ! Rebuild face flux from detailed face flux
               lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP
               do n=0,getSize(this%detailed_face_flux(1,i,j,k))-1
                  call getSepVMAtIndex(this%detailed_face_flux(1,i,j,k),n,my_SepVM)
                  lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0)
                  gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1)
               end do
               call construct(this%face_flux(1,i,j,k),[lvol,lbar,gvol,gbar])
            end if
            ! Y flux
            if (minval(abs(band(i,j-1:j,k))).le.advect_band) then
               ! Construct and project face
               face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,5)=back_project(face(:,1))
               face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=back_project(face(:,2))
               face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=back_project(face(:,3))
               face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,8)=back_project(face(:,4))
               face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]; face(:,9)=back_project(face(:,9))
               ! Form flux polyhedron
               call construct(flux_polyhedron,face)
               ! Add solenoidal correction
               call adjustCapToMatchVolume(flux_polyhedron,dt*V(i,j,k)*this%dx*this%dz)
               ! Need full geometric flux
               call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%detailed_face_flux(2,i,j,k))
               ! Rebuild face flux from detailed face flux
               lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP
               do n=0,getSize(this%detailed_face_flux(2,i,j,k))-1
                  call getSepVMAtIndex(this%detailed_face_flux(2,i,j,k),n,my_SepVM)
                  lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0)
                  gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1)
               end do
               call construct(this%face_flux(2,i,j,k),[lvol,lbar,gvol,gbar])
            end if
            ! Z flux
            if (minval(abs(band(i,j,k-1:k))).le.advect_band) then
               ! Construct and project face
               face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=back_project(face(:,1))
               face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=back_project(face(:,2))
               face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,7)=back_project(face(:,3))
               face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,8)=back_project(face(:,4))
               face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]; face(:,9)=back_project(face(:,9))
               ! Form flux polyhedron
               call construct(flux_polyhedron,face)
               ! Add solenoidal correction
               call adjustCapToMatchVolume(flux_polyhedron,dt*W(i,j,k)*this%dx*this%dy)
               ! Need full geometric flux
               call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%detailed_face_flux(3,i,j,k))
               ! Rebuild face flux from detailed face flux
               lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP
               do n=0,getSize(this%detailed_face_flux(3,i,j,k))-1
                  call getSepVMAtIndex(this%detailed_face_flux(3,i,j,k),n,my_SepVM)
                  lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0)
                  gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1)
               end do
               call construct(this%face_flux(3,i,j,k),[lvol,lbar,gvol,gbar])
            end if
         end do; end do; end do
      end block calculate_fluxes
      
      ! Compute transported moments
      update_moments: block
         real(WP) :: Lvolold,Gvolold,Lvolinc,Gvolinc,Lvolnew,Gvolnew
         ! Loop through each cell
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            ! Only work in band
            if (abs(band(i,j,k)).gt.advect_band) cycle
            ! Old liquid and gas volumes
            Lvolold=        this%VFold(i,j,k) *this%cfg%vol(i,j,k)
            Gvolold=(1.0_WP-this%VFold(i,j,k))*this%cfg%vol(i,j,k)
            ! Compute incoming liquid and gas volumes
            Lvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),0)+getVolumePtr(this%face_flux(1,i,j,k),0) &
            &       -getVolumePtr(this%face_flux(2,i,j+1,k),0)+getVolumePtr(this%face_flux(2,i,j,k),0) &
            &       -getVolumePtr(this%face_flux(3,i,j,k+1),0)+getVolumePtr(this%face_flux(3,i,j,k),0)
            Gvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),1)+getVolumePtr(this%face_flux(1,i,j,k),1) &
            &       -getVolumePtr(this%face_flux(2,i,j+1,k),1)+getVolumePtr(this%face_flux(2,i,j,k),1) &
            &       -getVolumePtr(this%face_flux(3,i,j,k+1),1)+getVolumePtr(this%face_flux(3,i,j,k),1)
            ! Compute new liquid and gas volumes
            Lvolnew=Lvolold+Lvolinc
            Gvolnew=Gvolold+Gvolinc
            ! Compute new liquid volume fraction
            this%VF(i,j,k)=Lvolnew/(Lvolnew+Gvolnew)
            ! Only work on higher order moments if VF is in [VFlo,VFhi]
            if (this%VF(i,j,k).lt.VFlo) then
               this%VF(i,j,k)=0.0_WP
            else if (this%VF(i,j,k).gt.VFhi) then
               this%VF(i,j,k)=1.0_WP
            else
               ! Compute old phase barycenters
               this%Lbary(:,i,j,k)=(this%Lbary(:,i,j,k)*Lvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),0)+getCentroidPtr(this%face_flux(1,i,j,k),0) &
               &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),0)+getCentroidPtr(this%face_flux(2,i,j,k),0) &
               &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),0)+getCentroidPtr(this%face_flux(3,i,j,k),0))/Lvolnew
               this%Gbary(:,i,j,k)=(this%Gbary(:,i,j,k)*Gvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),1)+getCentroidPtr(this%face_flux(1,i,j,k),1) &
               &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),1)+getCentroidPtr(this%face_flux(2,i,j,k),1) &
               &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),1)+getCentroidPtr(this%face_flux(3,i,j,k),1))/Gvolnew
               ! Project forward in time
               this%Lbary(:,i,j,k)=project(this%Lbary(:,i,j,k))
               this%Gbary(:,i,j,k)=project(this%Gbary(:,i,j,k))
            end if
         end do; end do; end do
      end block update_moments
      
      ! Synchronize and clean up barycenter
      sync_and_cleanup: block
         ! Synchronize VF
         call this%cfg%sync(this%VF)
         ! Clean up barycenters everywhere
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
               this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end if
         end do; end do; end do
         ! Synchronize barycenters
         call this%cfg%sync(this%Lbary)
         call this%cfg%sync(this%Gbary)
         ! Fix barycenter synchronization across periodic boundaries
         if (this%cfg%xper.and.this%cfg%iproc.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino,this%cfg%imin-1
               this%Lbary(1,i,j,k)=this%Lbary(1,i,j,k)-this%cfg%xL
               this%Gbary(1,i,j,k)=this%Gbary(1,i,j,k)-this%cfg%xL
            end do; end do; end do
         end if
         if (this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imax+1,this%cfg%imaxo
               this%Lbary(1,i,j,k)=this%Lbary(1,i,j,k)+this%cfg%xL
               this%Gbary(1,i,j,k)=this%Gbary(1,i,j,k)+this%cfg%xL
            end do; end do; end do
         end if
         if (this%cfg%yper.and.this%cfg%jproc.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino,this%cfg%jmin-1; do i=this%cfg%imino_,this%cfg%imaxo_
               this%Lbary(2,i,j,k)=this%Lbary(2,i,j,k)-this%cfg%yL
               this%Gbary(2,i,j,k)=this%Gbary(2,i,j,k)-this%cfg%yL
            end do; end do; end do
         end if
         if (this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmax+1,this%cfg%jmaxo; do i=this%cfg%imino_,this%cfg%imaxo_
               this%Lbary(2,i,j,k)=this%Lbary(2,i,j,k)+this%cfg%yL
               this%Gbary(2,i,j,k)=this%Gbary(2,i,j,k)+this%cfg%yL
            end do; end do; end do
         end if
         if (this%cfg%zper.and.this%cfg%kproc.eq.1) then
            do k=this%cfg%kmino,this%cfg%kmin-1; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               this%Lbary(3,i,j,k)=this%Lbary(3,i,j,k)-this%cfg%zL
               this%Gbary(3,i,j,k)=this%Gbary(3,i,j,k)-this%cfg%zL
            end do; end do; end do
         end if
         if (this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
            do k=this%cfg%kmax+1,this%cfg%kmaxo; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               this%Lbary(3,i,j,k)=this%Lbary(3,i,j,k)+this%cfg%zL
               this%Gbary(3,i,j,k)=this%Gbary(3,i,j,k)+this%cfg%zL
            end do; end do; end do
         end if
         ! Handle 2D barycenters
         if (this%cfg%nx.eq.1) then
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%Lbary(1,i,:,:)=this%cfg%xm(i)
               this%Gbary(1,i,:,:)=this%cfg%xm(i)
            end do
         end if
         if (this%cfg%ny.eq.1) then
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               this%Lbary(2,:,j,:)=this%cfg%ym(j)
               this%Gbary(2,:,j,:)=this%cfg%ym(j)
            end do
         end if
         if (this%cfg%nz.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_
               this%Lbary(3,:,:,k)=this%cfg%zm(k)
               this%Gbary(3,:,:,k)=this%cfg%zm(k)
            end do
         end if
      end block sync_and_cleanup
      
      ! Reconstruct the interface
      call this%build_interface()
      
      ! Deallocate the band
      deallocate(band)
      
   contains
      !> Project function that moves a point p1 in vicinity of (i,j,k) at (U,V,W) for -dt
      function back_project(p1) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), dimension(3) :: p2
         p2=p1-dt*this%cfg%get_velocity(        p1    ,i,j,k,U,V,W)
         p2=p1-dt*this%cfg%get_velocity(0.5_WP*(p1+p2),i,j,k,U,V,W)
      end function back_project
      !> Project function that moves a point p1 in vicinity of (i,j,k) at (U,V,W) for dt
      function project(p1) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), dimension(3) :: p2
         p2=p1+dt*this%cfg%get_velocity(        p1    ,i,j,k,U,V,W)
         p2=p1+dt*this%cfg%get_velocity(0.5_WP*(p1+p2),i,j,k,U,V,W)
      end function project
   end subroutine advance_volume
   
   
   !> Build the interface from volume moments using PLICnet
   subroutine build_interface(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      
      ! Perfrom PLICnet reconstruction of the interface
      plicnet_reconstruct: block
         use mathtools, only: normalize
         use plicnet,   only: get_normal,reflect_moments
         integer(IRL_SignedIndex_t) :: i,j,k
         integer :: ind,ii,jj,kk
         real(IRL_double), dimension(0:2) :: normal
         real(IRL_double), dimension(0:188) :: moments
         integer :: direction
         logical :: flip
         real(IRL_double) :: m000,m100,m010,m001
         real(IRL_double), dimension(0:2) :: center
         real(IRL_double) :: initial_dist
         type(RectCub_type) :: cell
         ! Get a cell
         call new(cell)
         ! Traverse domain and reconstruct interface
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            ! Handle full cells differently
            if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
               call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
               call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
               cycle
            end if
            ! Liquid-gas symmetry
            flip=.false.
            if (this%VF(i,j,k).ge.0.5_WP) flip=.true.
            m000=0; m100=0; m010=0; m001=0
            ! Construct neighborhood of volume moments
            if (flip) then
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=1.0_WP-this%VF(ii,jj,kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%Gbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%Gbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%Gbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%Lbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%Lbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%Lbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  ! Calculate geometric moments of neighborhood
                  m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
               end do; end do; end do
            else
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=this%VF(ii,jj,kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%Lbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%Lbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%Lbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%Gbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%Gbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%Gbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  ! Calculate geometric moments of neighborhood
                  m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
               end do; end do; end do
            end if
            ! Calculate geometric center of neighborhood
            center=[m100,m010,m001]/m000
            ! Symmetry about Cartesian planes
            call reflect_moments(moments,center,direction)
            ! Get PLIC normal vector from neural network
            call get_normal(moments,normal)
            normal=normalize(normal)
            ! Rotate normal vector to original octant
            if (direction.eq.1) then
               normal(0)=-normal(0)
            else if (direction.eq.2) then
               normal(1)=-normal(1)
            else if (direction.eq.3) then
               normal(2)=-normal(2)
            else if (direction.eq.4) then
               normal(0)=-normal(0)
               normal(1)=-normal(1)
            else if (direction.eq.5) then
               normal(0)=-normal(0)
               normal(2)=-normal(2)
            else if (direction.eq.6) then
               normal(1)=-normal(1)
               normal(2)=-normal(2)
            else if (direction.eq.7) then
               normal(0)=-normal(0)
               normal(1)=-normal(1)
               normal(2)=-normal(2)
            end if
            if (.not.flip) then
               normal(0)=-normal(0)
               normal(1)=-normal(1)
               normal(2)=-normal(2)
            end if
            ! Locate PLIC plane in cell
            call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
            initial_dist=dot_product(normal,[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)])
            call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
            call setPlane(this%liquid_gas_interface(i,j,k),0,normal,initial_dist)
            call matchVolumeFraction(cell,this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
         end do; end do; end do
      end block plicnet_reconstruct
      
      ! Synchronize PLIC interface across boundaries
      synchronize_interface: block
         integer :: i,j,k,ni
         real(WP), dimension(1:4) :: plane
         integer , dimension(2,3) :: send_range,recv_range
         ! Synchronize in x
         if (this%cfg%nx.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               call copy(this%liquid_gas_interface(i,j,k),this%liquid_gas_interface(this%cfg%imin,j,k))
            end do; end do; end do
         else
            ! Send minus
            send_range(1:2,1)=[this%cfg%imin_   ,this%cfg%imin_ +this%cfg%no-1]
            send_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
            send_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
            recv_range(1:2,1)=[this%cfg%imax_ +1,this%cfg%imaxo_              ]
            recv_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
            recv_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
            call sync_side(send_range,recv_range,0,-1)
            ! Send plus
            send_range(1:2,1)=[this%cfg%imax_ -this%cfg%no+1,this%cfg%imax_   ]
            send_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
            send_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
            recv_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imin_ -1]
            recv_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
            recv_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
            call sync_side(send_range,recv_range,0,+1)
         end if
         ! Synchronize in y
         if (this%cfg%ny.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               call copy(this%liquid_gas_interface(i,j,k),this%liquid_gas_interface(i,this%cfg%jmin,k))
            end do; end do; end do
         else
            ! Send minus side
            send_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
            send_range(1:2,2)=[this%cfg%jmin_   ,this%cfg%jmin_ +this%cfg%no-1]
            send_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
            recv_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
            recv_range(1:2,2)=[this%cfg%jmax_ +1,this%cfg%jmaxo_              ]
            recv_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
            call sync_side(send_range,recv_range,1,-1)
            ! Send plus side
            send_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
            send_range(1:2,2)=[this%cfg%jmax_ -this%cfg%no+1,this%cfg%jmax_   ]
            send_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
            recv_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
            recv_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmin_ -1]
            recv_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
            call sync_side(send_range,recv_range,1,+1)
         end if
         ! Synchronize in z
         if (this%cfg%nz.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               call copy(this%liquid_gas_interface(i,j,k),this%liquid_gas_interface(i,j,this%cfg%kmin))
            end do; end do; end do
         else
            ! Send minus side
            send_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
            send_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
            send_range(1:2,3)=[this%cfg%kmin_   ,this%cfg%kmin_ +this%cfg%no-1]
            recv_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
            recv_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
            recv_range(1:2,3)=[this%cfg%kmax_ +1,this%cfg%kmaxo_              ]
            call sync_side(send_range,recv_range,2,-1)
            ! Send plus side
            send_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
            send_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
            send_range(1:2,3)=[this%cfg%kmax_ -this%cfg%no+1,this%cfg%kmax_   ]
            recv_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
            recv_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
            recv_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmin_ -1]
            call sync_side(send_range,recv_range,2,+1)
         end if
         ! Fix plane position if we are periodic in x
         if (this%cfg%xper.and.this%cfg%iproc.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino,this%cfg%imin-1
               plane=getPlane(this%liquid_gas_interface(i,j,k),0)
               plane(4)=plane(4)-plane(1)*this%cfg%xL
               call setPlane(this%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         if (this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imax+1,this%cfg%imaxo
               plane=getPlane(this%liquid_gas_interface(i,j,k),0)
               plane(4)=plane(4)+plane(1)*this%cfg%xL
               call setPlane(this%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         ! Fix plane position if we are periodic in y
         if (this%cfg%yper.and.this%cfg%jproc.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino,this%cfg%jmin-1; do i=this%cfg%imino_,this%cfg%imaxo_
               plane=getPlane(this%liquid_gas_interface(i,j,k),0)
               plane(4)=plane(4)-plane(2)*this%cfg%yL
               call setPlane(this%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         if (this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmax+1,this%cfg%jmaxo; do i=this%cfg%imino_,this%cfg%imaxo_
               plane=getPlane(this%liquid_gas_interface(i,j,k),0)
               plane(4)=plane(4)+plane(2)*this%cfg%yL
               call setPlane(this%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         ! Fix plane position if we are periodic in z
         if (this%cfg%zper.and.this%cfg%kproc.eq.1) then
            do k=this%cfg%kmino,this%cfg%kmin-1; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               plane=getPlane(this%liquid_gas_interface(i,j,k),0)
               plane(4)=plane(4)-plane(3)*this%cfg%zL
               call setPlane(this%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         if (this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
            do k=this%cfg%kmax+1,this%cfg%kmaxo; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               plane=getPlane(this%liquid_gas_interface(i,j,k),0)
               plane(4)=plane(4)+plane(3)*this%cfg%zL
               call setPlane(this%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
      end block synchronize_interface
      
      ! Reset volume moments to match with reconstructed interface
      reset_moments: block
         integer :: i,j,k
         type(RectCub_type) :: cell
         type(SepVM_type) :: separated_volume_moments
         call new(cell)
         call new(separated_volume_moments)
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            ! Form the grid cell
            call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
            ! Cut it by the current interface
            call getNormMoments(cell,this%liquid_gas_interface(i,j,k),separated_volume_moments)
            ! Recover relevant moments
            this%VF(i,j,k)     =getVolumePtr(separated_volume_moments,0)/this%cfg%vol(i,j,k)
            this%Lbary(:,i,j,k)=getCentroid(separated_volume_moments,0)
            this%Gbary(:,i,j,k)=getCentroid(separated_volume_moments,1)
            ! Clean up
            if (this%VF(i,j,k).lt.VFlo) then
               this%VF(i,j,k)     =0.0_WP
               this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end if
            if (this%VF(i,j,k).gt.VFhi) then
               this%VF(i,j,k)     =1.0_WP
               this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end if
         end do; end do; end do
      end block reset_moments
      
      ! Polygonize interface
      polygonize_interface: block
         integer :: i,j,k,n
         type(RectCub_type) :: cell
         real(WP), dimension(1:3,1:4) :: vert
         real(WP), dimension(1:3) :: norm
         ! Create a cell object
         call new(cell)
         ! Loop over full domain and form polygon
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            ! Zero out the polygon
            call zeroPolygon(this%interface_polygon(i,j,k))
            ! Create polygons for cells with interfaces, zero for those without
            if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).le.VFhi) then
               call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
               call getPoly(cell,this%liquid_gas_interface(i,j,k),0,this%interface_polygon(i,j,k))
            end if
         end do; end do; end do
         ! Find inferface between filled and empty cells on x-face
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_+1,this%cfg%imaxo_
            if (this%VF(i,j,k).lt.VFlo.and.this%VF(i-1,j,k).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(i-1,j,k).lt.VFlo) then
               norm=[sign(1.0_WP,0.5_WP-this%VF(i,j,k)),0.0_WP,0.0_WP]
               call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
               call setPlane(this%liquid_gas_interface(i,j,k),0,norm,sign(1.0_WP,0.5_WP-this%VF(i,j,k))*this%cfg%x(i))
               vert(:,1)=[this%cfg%x(i),this%cfg%y(j  ),this%cfg%z(k  )]
               vert(:,2)=[this%cfg%x(i),this%cfg%y(j+1),this%cfg%z(k  )]
               vert(:,3)=[this%cfg%x(i),this%cfg%y(j+1),this%cfg%z(k+1)]
               vert(:,4)=[this%cfg%x(i),this%cfg%y(j  ),this%cfg%z(k+1)]
               call construct(this%interface_polygon(i,j,k),4,vert)
               call setPlaneOfExistence(this%interface_polygon(i,j,k),getPlane(this%liquid_gas_interface(i,j,k),0))
               if (this%VF(i,j,k).gt.VFhi) call reversePtOrdering(this%interface_polygon(i,j,k))
            end if
         end do; end do; end do
         ! Find inferface between filled and empty cells on y-face
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_+1,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            if (this%VF(i,j,k).lt.VFlo.and.this%VF(i,j-1,k).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(i,j-1,k).lt.VFlo) then
               norm=[0.0_WP,sign(1.0_WP,0.5_WP-this%VF(i,j,k)),0.0_WP]
               call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
               call setPlane(this%liquid_gas_interface(i,j,k),0,norm,sign(1.0_WP,0.5_WP-this%VF(i,j,k))*this%cfg%y(j))
               vert(:,1)=[this%cfg%x(i  ),this%cfg%y(j),this%cfg%z(k  )]
               vert(:,2)=[this%cfg%x(i  ),this%cfg%y(j),this%cfg%z(k+1)]
               vert(:,3)=[this%cfg%x(i+1),this%cfg%y(j),this%cfg%z(k+1)]
               vert(:,4)=[this%cfg%x(i+1),this%cfg%y(j),this%cfg%z(k  )]
               call construct(this%interface_polygon(i,j,k),4,vert)
               call setPlaneOfExistence(this%interface_polygon(i,j,k),getPlane(this%liquid_gas_interface(i,j,k),0))
               if (this%VF(i,j,k).gt.VFhi) call reversePtOrdering(this%interface_polygon(i,j,k))
            end if
         end do; end do; end do
         ! Find inferface between filled and empty cells on z-face
         do k=this%cfg%kmino_+1,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            if (this%VF(i,j,k).lt.VFlo.and.this%VF(i,j,k-1).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(i,j,k-1).lt.VFlo) then
               norm=[0.0_WP,0.0_WP,sign(1.0_WP,0.5_WP-this%VF(i,j,k))]
               call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
               call setPlane(this%liquid_gas_interface(i,j,k),0,norm,sign(1.0_WP,0.5_WP-this%VF(i,j,k))*this%cfg%z(k))
               vert(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k)]
               vert(:,2)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k)]
               vert(:,3)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k)]
               vert(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k)]
               call construct(this%interface_polygon(i,j,k),4,vert)
               call setPlaneOfExistence(this%interface_polygon(i,j,k),getPlane(this%liquid_gas_interface(i,j,k),0))
               if (this%VF(i,j,k).gt.VFhi) call reversePtOrdering(this%interface_polygon(i,j,k))
            end if
         end do; end do; end do
      end block polygonize_interface
      
   contains
      
      !> Private procedure to perform communication across one boundary
      subroutine sync_side(a_send_range,a_recv_range,a_dimension,a_direction)
         implicit none
         integer, dimension(2,3), intent(in) :: a_send_range
         integer, dimension(2,3), intent(in) :: a_recv_range
         integer, intent(in) :: a_dimension
         integer, intent(in) :: a_direction
         integer :: i,j,k
         logical :: something_received
         ! Pack the buffer
         call resetBufferPointer(this%send_byte_buffer)
         call setSize(this%send_byte_buffer,int(0,8))
         do k=a_send_range(1,3),a_send_range(2,3); do j=a_send_range(1,2),a_send_range(2,2); do i=a_send_range(1,1),a_send_range(2,1)
            call serializeAndPack(this%liquid_gas_interface(i,j,k),this%send_byte_buffer)
         end do; end do; end do
         ! Communicate
         call sync_ByteBuffer(this%send_byte_buffer,a_dimension,a_direction,this%recv_byte_buffer,something_received)
         ! If something was received, unpack it: traversal order is important and must be aligned with how the sent data was packed
         if (something_received) then
            call resetBufferPointer(this%recv_byte_buffer)
            do k=a_recv_range(1,3),a_recv_range(2,3); do j=a_recv_range(1,2),a_recv_range(2,2); do i=a_recv_range(1,1),a_recv_range(2,1)
               call unpackAndStore(this%liquid_gas_interface(i,j,k),this%recv_byte_buffer)
            end do; end do; end do
         end if
      end subroutine sync_side
      
      !> Private procedure to communicate a package of bytes across one boundary
      subroutine sync_ByteBuffer(a_send_buffer,a_dimension,a_direction,a_receive_buffer,a_received_something)
         use mpi_f08
         implicit none
         type(ByteBuffer_type), intent(inout)  :: a_send_buffer !< Inout needed because it is preallocated
         integer, intent(in) :: a_dimension  !< Should be 0/1/2 for x/y/z
         integer, intent(in) :: a_direction  !< Should be -1 for left or +1 for right
         type(ByteBuffer_type), intent(inout) :: a_receive_buffer !< Inout needed because it is preallocated
         logical, intent(out) :: a_received_something
         type(MPI_Status) :: status
         integer :: isrc,idst,ierr
         integer(IRL_LargeOffsetIndex_t) :: my_size
         integer(IRL_LargeOffsetIndex_t) :: incoming_size
         integer :: my_size_small,incoming_size_small
         ! Figure out source and destination
         call MPI_CART_SHIFT(this%cfg%comm,a_dimension,a_direction,isrc,idst,ierr)
         ! Communicate sizes so that each processor knows what to expect in main communication
         my_size=getSize(a_send_buffer)
         call MPI_SENDRECV(my_size,1,MPI_INTEGER8,idst,0,incoming_size,1,MPI_INTEGER8,isrc,0,this%cfg%comm,status,ierr)
         ! Set size of recv buffer to appropriate size and perform send-receive
         if (isrc.ne.MPI_PROC_NULL) then
            a_received_something=.true.
            call setSize(a_receive_buffer,incoming_size)
         else
            a_received_something=.false.
            incoming_size=0
            call setSize(a_receive_buffer,int(1,8))
         end if
         ! Convert integers
         my_size_small=int(my_size,4)
         incoming_size_small=int(incoming_size,4)
         call MPI_SENDRECV(dataPtr(a_send_buffer),my_size_small,MPI_BYTE,idst,0,dataPtr(a_receive_buffer),incoming_size_small,MPI_BYTE,isrc,0,this%cfg%comm,status,ierr)
      end subroutine sync_ByteBuffer
      
   end subroutine build_interface
   
   
   !> Calculate the explicit time derivative of conserved variables given primitive variables
   subroutine get_RHS(this,dVARdt)
      use irl_fortran_interface
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(out) :: dVARdt  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nVAR)
      real(WP), dimension(:,:,:,:), allocatable :: FX,FY,FZ
      real(WP), dimension(:,:,:), pointer :: buf
      integer :: i,j,k,n
      real(WP) :: w,div
      real(WP), parameter :: eps=1.0e-15_WP
      real(WP), dimension(-2: 0) :: wxp,wyp,wzp
      real(WP), dimension(-1:+1) :: wxm,wym,wzm
      type(SepVM_type) :: my_SepVM
      integer, dimension(3) :: ind
      
      ! Zero out array
      dVARdt=0.0_WP
      
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nVAR))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nVAR))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nVAR))
      
      ! ================================================================ !
      ! ======================== INVISID FLUXES ======================== !
      ! ================================================================ !
      
      ! Calculate face-centered mass fluxes
      do k=this%cfg%kmin_-1,this%cfg%kmax_+1
         do j=this%cfg%jmin_-1,this%cfg%jmax_+1
            do i=this%cfg%imin_-1,this%cfg%imax_+1
               ! Mass flux in X
               if (getSize(this%detailed_face_flux(1,i,j,k)).gt.0) then
                  ! Detailed geometric flux is available, use geometric fluxing to accumulate mass flux
                  FX(i,j,k,1)=0.0_WP
                  do n=0,getSize(this%detailed_face_flux(1,i,j,k))-1
                     ind=this%cfg%get_ijk_from_lexico(getTagForIndex(this%detailed_face_flux(1,i,j,k),n))
                     call getSepVMAtIndex(this%detailed_face_flux(1,i,j,k),n,my_SepVM)
                     FX(i,j,k,1)=FX(i,j,k,1)-getVolume(my_SepVM,0)*this%RHOl(ind(1),ind(2),ind(3))
                  end do
                  FX(i,j,k,1)=FX(i,j,k,1)/(this%dt_detailed_face_flux*this%dy*this%dz)
               else
                  ! Only calculate in liquid
                  if (this%VF(i,j,k).gt.0.0_WP) then
                     ! No detailed geometric flux is available, prepare WENO scheme for mass transport
                     w=weno_weight((abs(this%RHOl(i-1,j,k)-this%RHOl(i-2,j,k))+eps)/(abs(this%RHOl(i,j,k)-this%RHOl(i-1,j,k))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOl(i+1,j,k)-this%RHOl(i  ,j,k))+eps)/(abs(this%RHOl(i,j,k)-this%RHOl(i-1,j,k))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     ! WENO-based mass flux
                     FX(i,j,k,1)=-0.5_WP*(this%U(i,j,k)+abs(this%U(i,j,k)))*sum(wxp*this%RHOl(i-2:i  ,j,k)) &
                     &           -0.5_WP*(this%U(i,j,k)-abs(this%U(i,j,k)))*sum(wxm*this%RHOl(i-1:i+1,j,k))
                  else
                     FX(i,j,k,1)=0.0_WP
                  end if
               end if
               ! Mass flux in Y
               if (getSize(this%detailed_face_flux(2,i,j,k)).gt.0) then
                  ! Detailed geometric flux is available, use geometric fluxing to accumulate mass flux
                  FY(i,j,k,1)=0.0_WP
                  do n=0,getSize(this%detailed_face_flux(2,i,j,k))-1
                     ind=this%cfg%get_ijk_from_lexico(getTagForIndex(this%detailed_face_flux(2,i,j,k),n))
                     call getSepVMAtIndex(this%detailed_face_flux(2,i,j,k),n,my_SepVM)
                     FY(i,j,k,1)=FY(i,j,k,1)-getVolume(my_SepVM,0)*this%RHOl(ind(1),ind(2),ind(3))
                  end do
                  FY(i,j,k,1)=FY(i,j,k,1)/(this%dt_detailed_face_flux*this%dx*this%dz)
               else
                  ! Only calculate in liquid
                  if (this%VF(i,j,k).gt.0.0_WP) then
                     ! No detailed geometric flux is available, prepare WENO scheme for mass transport
                     w=weno_weight((abs(this%RHOl(i,j-1,k)-this%RHOl(i,j-2,k))+eps)/(abs(this%RHOl(i,j,k)-this%RHOl(i,j-1,k))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOl(i,j+1,k)-this%RHOl(i,j  ,k))+eps)/(abs(this%RHOl(i,j,k)-this%RHOl(i,j-1,k))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     ! WENO-based mass flux
                     FY(i,j,k,1)=-0.5_WP*(this%V(i,j,k)+abs(this%V(i,j,k)))*sum(wyp*this%RHOl(i,j-2:j  ,k)) &
                     &           -0.5_WP*(this%V(i,j,k)-abs(this%V(i,j,k)))*sum(wym*this%RHOl(i,j-1:j+1,k))
                  else
                     FY(i,j,k,1)=0.0_WP
                  end if
               end if
               ! Mass flux in Z
               if (getSize(this%detailed_face_flux(3,i,j,k)).gt.0) then
                  ! Detailed geometric flux is available, use geometric fluxing to accumulate mass flux
                  FZ(i,j,k,1)=0.0_WP
                  do n=0,getSize(this%detailed_face_flux(3,i,j,k))-1
                     ind=this%cfg%get_ijk_from_lexico(getTagForIndex(this%detailed_face_flux(3,i,j,k),n))
                     call getSepVMAtIndex(this%detailed_face_flux(2,i,j,k),n,my_SepVM)
                     FZ(i,j,k,1)=FZ(i,j,k,1)-getVolume(my_SepVM,0)*this%RHOl(ind(1),ind(2),ind(3))
                  end do
                  FZ(i,j,k,1)=FZ(i,j,k,1)/(this%dt_detailed_face_flux*this%dx*this%dy)
               else
                  ! Only calculate in liquid
                  if (this%VF(i,j,k).gt.0.0_WP) then
                     ! No detailed geometric flux is available, prepare WENO scheme for mass transport
                     w=weno_weight((abs(this%RHOl(i,j,k-1)-this%RHOl(i,j,k-2))+eps)/(abs(this%RHOl(i,j,k)-this%RHOl(i,j,k-1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOl(i,j,k+1)-this%RHOl(i,j,k  ))+eps)/(abs(this%RHOl(i,j,k)-this%RHOl(i,j,k-1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     ! WENO-based mass fluxes
                     FZ(i,j,k,1)=-0.5_WP*(this%W(i,j,k)+abs(this%W(i,j,k)))*sum(wzp*this%RHOl(i,j,k-2:k  )) &
                     &           -0.5_WP*(this%W(i,j,k)-abs(this%W(i,j,k)))*sum(wzm*this%RHOl(i,j,k-1:k+1))
                  else
                     FZ(i,j,k,1)=0.0_WP
                  end if
               end if
            end do
         end do
      end do
      
      ! Calculate cell-centered momentum fluxes
      !do k=this%cfg%kmin_-1,this%cfg%kmax_
      !   do j=this%cfg%jmin_-1,this%cfg%jmax_
      !      do i=this%cfg%imin_-1,this%cfg%imax_
      !         FUX(i,j,k)=0.25_WP*sum(FRX(i:i+1,j,k))*sum(this%U(i:i+1,j,k))-this%P(i,j,k)
      !         FVY(i,j,k)=0.25_WP*sum(FRY(i,j:j+1,k))*sum(this%V(i,j:j+1,k))-this%P(i,j,k)
      !         FWZ(i,j,k)=0.25_WP*sum(FRZ(i,j,k:k+1))*sum(this%W(i,j,k:k+1))-this%P(i,j,k)
      !      end do
      !   end do
      !end do
      
      ! Calculate edge-centered momentum fluxes and face-centered internal energy fluxes
      !do k=this%cfg%kmin_,this%cfg%kmax_+1
      !   do j=this%cfg%jmin_,this%cfg%jmax_+1
      !      do i=this%cfg%imin_,this%cfg%imax_+1
      !         ! Momentum fluxes
      !         FUY(i,j,k)=0.25_WP*sum(FRY(i-1:i,j,k))*sum(this%U(i,j-1:j,k))
      !         FUZ(i,j,k)=0.25_WP*sum(FRZ(i-1:i,j,k))*sum(this%U(i,j,k-1:k))
      !         FVX(i,j,k)=0.25_WP*sum(FRX(i,j-1:j,k))*sum(this%V(i-1:i,j,k))
      !         FVZ(i,j,k)=0.25_WP*sum(FRZ(i,j-1:j,k))*sum(this%V(i,j,k-1:k))
      !         FWX(i,j,k)=0.25_WP*sum(FRX(i,j,k-1:k))*sum(this%W(i-1:i,j,k))
      !         FWY(i,j,k)=0.25_WP*sum(FRY(i,j,k-1:k))*sum(this%W(i,j-1:j,k))
               ! Prepare WENO scheme for internal energy transport
               !w=weno_weight((abs(this%E(i-1,j,k)-this%E(i-2,j,k))+eps)/(abs(this%E(i,j,k)-this%E(i-1,j,k))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               !w=weno_weight((abs(this%E(i+1,j,k)-this%E(i  ,j,k))+eps)/(abs(this%E(i,j,k)-this%E(i-1,j,k))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               !w=weno_weight((abs(this%E(i,j-1,k)-this%E(i,j-2,k))+eps)/(abs(this%E(i,j,k)-this%E(i,j-1,k))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               !w=weno_weight((abs(this%E(i,j+1,k)-this%E(i,j  ,k))+eps)/(abs(this%E(i,j,k)-this%E(i,j-1,k))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               !w=weno_weight((abs(this%E(i,j,k-1)-this%E(i,j,k-2))+eps)/(abs(this%E(i,j,k)-this%E(i,j,k-1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               !w=weno_weight((abs(this%E(i,j,k+1)-this%E(i,j,k  ))+eps)/(abs(this%E(i,j,k)-this%E(i,j,k-1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               ! WENO-based internal energy fluxes
               !FEX(i,j,k)=0.5_WP*(FRX(i,j,k)-abs(-FRX(i,j,k)))*sum(wxp*this%E(i-2:i  ,j,k)) &
               !&         +0.5_WP*(FRX(i,j,k)+abs(-FRX(i,j,k)))*sum(wxm*this%E(i-1:i+1,j,k))
               !FEY(i,j,k)=0.5_WP*(FRY(i,j,k)-abs(-FRY(i,j,k)))*sum(wyp*this%E(i,j-2:j  ,k)) &
               !&         +0.5_WP*(FRY(i,j,k)+abs(-FRY(i,j,k)))*sum(wym*this%E(i,j-1:j+1,k))
               !FEZ(i,j,k)=0.5_WP*(FRZ(i,j,k)-abs(-FRZ(i,j,k)))*sum(wzp*this%E(i,j,k-2:k  )) &
               !&         +0.5_WP*(FRZ(i,j,k)+abs(-FRZ(i,j,k)))*sum(wzm*this%E(i,j,k-1:k+1))
      !         ! Centered internal energy fluxes
      !         FEX(i,j,k)=0.5_WP*FRX(i,j,k)*sum(this%E(i-1:i,j,k))
      !         FEY(i,j,k)=0.5_WP*FRY(i,j,k)*sum(this%E(i,j-1:j,k))
      !         FEZ(i,j,k)=0.5_WP*FRZ(i,j,k)*sum(this%E(i,j,k-1:k))
      !      end do
      !   end do
      !end do
      
      ! Assemble time derivative
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Divergence of fluxes
               dVARdt(i,j,k,1)=this%dxi*(FX(i+1,j,k,1)-FX(i  ,j,k,1))+this%dyi*(FY(i,j+1,k,1)-FY(i,j  ,k,1))+this%dzi*(FZ(i,j,k+1,1)-FZ(i,j,k  ,1))
               !dRHOUdt(i,j,k)=this%dxi*(FUX(i  ,j,k)-FUX(i-1,j,k))+this%dyi*(FUY(i,j+1,k)-FUY(i,j  ,k))+this%dzi*(FUZ(i,j,k+1)-FUZ(i,j,k  ))
               !dRHOVdt(i,j,k)=this%dxi*(FVX(i+1,j,k)-FVX(i  ,j,k))+this%dyi*(FVY(i,j  ,k)-FVY(i,j-1,k))+this%dzi*(FVZ(i,j,k+1)-FVZ(i,j,k  ))
               !dRHOWdt(i,j,k)=this%dxi*(FWX(i+1,j,k)-FWX(i  ,j,k))+this%dyi*(FWY(i,j+1,k)-FWY(i,j  ,k))+this%dzi*(FWZ(i,j,k  )-FWZ(i,j,k-1))
               !dRHOEdt(i,j,k)=this%dxi*(FEX(i+1,j,k)-FEX(i  ,j,k))+this%dyi*(FEY(i,j+1,k)-FEY(i,j  ,k))+this%dzi*(FEZ(i,j,k+1)-FEZ(i,j,k  ))
               ! Pressure dilatation term
               !dRHOEdt(i,j,k)=dRHOEdt(i,j,k)-this%P(i,j,k)*(this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k)))
            end do
         end do
      end do
      
      ! ================================================================ !
      ! ======================== VISCOUS FLUXES ======================== !
      ! ================================================================ !
      
      ! Compute cell-centered fluxes
      !do k=this%cfg%kmin_-1,this%cfg%kmax_
      !   do j=this%cfg%jmin_-1,this%cfg%jmax_
      !      do i=this%cfg%imin_-1,this%cfg%imax_
      !         ! Divergence of velocity
      !         div=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
      !         ! Viscous flux of momentum
      !         FUX(i,j,k)=2.0_WP*this%visc(i,j,k)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+(this%beta(i,j,k)-2.0_WP*this%visc(i,j,k)/3.0_WP)*div
      !         FVY(i,j,k)=2.0_WP*this%visc(i,j,k)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+(this%beta(i,j,k)-2.0_WP*this%visc(i,j,k)/3.0_WP)*div
      !         FWZ(i,j,k)=2.0_WP*this%visc(i,j,k)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+(this%beta(i,j,k)-2.0_WP*this%visc(i,j,k)/3.0_WP)*div
      !      end do
      !   end do
      !end do
      
      ! Calculate edge-centered momentum fluxes and face-centered internal energy fluxes
      !do k=this%cfg%kmin_,this%cfg%kmax_+1
      !   do j=this%cfg%jmin_,this%cfg%jmax_+1
      !      do i=this%cfg%imin_,this%cfg%imax_+1
      !         ! Momentum fluxes (symmetric)
      !         FUY(i,j,k)=0.25_WP*sum(this%visc(i-1:i,j-1:j,k))*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k)))
      !         FVZ(i,j,k)=0.25_WP*sum(this%visc(i,j-1:j,k-1:k))*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k)))
      !         FUZ(i,j,k)=0.25_WP*sum(this%visc(i-1:i,j,k-1:k))*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1)))
      !         FVX(i,j,k)=FUY(i,j,k)
      !         FWY(i,j,k)=FVZ(i,j,k)
      !         FWX(i,j,k)=FUZ(i,j,k)
      !         ! Internal energy fluxes
      !         FEX(i,j,k)=0.5_WP*sum(this%diff(i-1:i,j,k))*this%dxi*(this%T(i,j,k)-this%T(i-1,j,k))
      !         FEY(i,j,k)=0.5_WP*sum(this%diff(i,j-1:j,k))*this%dyi*(this%T(i,j,k)-this%T(i,j-1,k))
      !         FEZ(i,j,k)=0.5_WP*sum(this%diff(i,j,k-1:k))*this%dzi*(this%T(i,j,k)-this%T(i,j,k-1))
      !         ! Also prepare viscous heating term
      !         FRZ(i,j,k)=FUY(i,j,k)*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k)))
      !         FRX(i,j,k)=FVZ(i,j,k)*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k)))
      !         FRY(i,j,k)=FUZ(i,j,k)*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1)))
      !      end do
      !   end do
      !end do
      
      ! Increment time derivative
      !do k=this%cfg%kmin_,this%cfg%kmax_
      !   do j=this%cfg%jmin_,this%cfg%jmax_
      !      do i=this%cfg%imin_,this%cfg%imax_
      !         ! Divergence of fluxes
      !         dRHOUdt(i,j,k)=dRHOUdt(i,j,k)+this%dxi*(FUX(i  ,j,k)-FUX(i-1,j,k))+this%dyi*(FUY(i,j+1,k)-FUY(i,j  ,k))+this%dzi*(FUZ(i,j,k+1)-FUZ(i,j,k  ))
      !         dRHOVdt(i,j,k)=dRHOVdt(i,j,k)+this%dxi*(FVX(i+1,j,k)-FVX(i  ,j,k))+this%dyi*(FVY(i,j  ,k)-FVY(i,j-1,k))+this%dzi*(FVZ(i,j,k+1)-FVZ(i,j,k  ))
      !         dRHOWdt(i,j,k)=dRHOWdt(i,j,k)+this%dxi*(FWX(i+1,j,k)-FWX(i  ,j,k))+this%dyi*(FWY(i,j+1,k)-FWY(i,j  ,k))+this%dzi*(FWZ(i,j,k  )-FWZ(i,j,k-1))
      !         dRHOEdt(i,j,k)=dRHOEdt(i,j,k)+this%dxi*(FEX(i+1,j,k)-FEX(i  ,j,k))+this%dyi*(FEY(i,j+1,k)-FEY(i,j  ,k))+this%dzi*(FEZ(i,j,k+1)-FEZ(i,j,k  ))
      !         ! Viscous heating term
      !         dRHOEdt(i,j,k)=dRHOEdt(i,j,k)+FUX(i,j,k)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+&
      !         &                             FVY(i,j,k)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+&
      !         &                             FWZ(i,j,k)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+&
      !         &                             0.25_WP*sum(FRZ(i:i+1,j:j+1,k))+0.25_WP*sum(FRX(i,j:j+1,k:k+1))+0.25_WP*sum(FRY(i:i+1,j,k:k+1))
      !      end do
      !   end do
      !end do
      
      ! Synchronize
      do n=1,this%nVAR
         call this%cfg%sync(dVARdt(:,:,:,n))
      end do
      
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      
   contains
      real(WP) function weno_weight(ratio)
         implicit none
         real(WP), intent(in) :: ratio
         real(WP), parameter :: lambda=0.13_WP ! Switching parameter
         real(WP), parameter :: delta=0.01_WP  ! Switching thickness
         weno_weight=(1.0_WP-tanh((ratio-lambda)/delta))/3.0_WP+(1.0_WP-tanh((ratio-1.0_WP/lambda)/delta))/6.0_WP
      end function weno_weight
   end subroutine get_RHS
   
   
   !> Calculate velocity from momentum
   subroutine get_velocity(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k
      ! Calculate velocity as far as possible
      !do k=this%cfg%kmino_+1,this%cfg%kmaxo_
      !   do j=this%cfg%jmino_+1,this%cfg%jmaxo_
      !      do i=this%cfg%imino_+1,this%cfg%imaxo_
      !         this%U(i,j,k)=2.0_WP*this%RHOU(i,j,k)/(this%RHO(i-1,j,k)+this%RHO(i,j,k))
      !         this%V(i,j,k)=2.0_WP*this%RHOV(i,j,k)/(this%RHO(i,j-1,k)+this%RHO(i,j,k))
      !         this%W(i,j,k)=2.0_WP*this%RHOW(i,j,k)/(this%RHO(i,j,k-1)+this%RHO(i,j,k))
      !      end do
      !   end do
      !end do
      ! Sync velocity
      !call this%cfg%sync(this%U)
      !call this%cfg%sync(this%V)
      !call this%cfg%sync(this%W)
      ! Add last layer in each direction
      !if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then
      !   this%U(this%cfg%imino,:,:)=this%RHOU(this%cfg%imino,:,:)/this%RHO(this%cfg%imino,:,:)
      !   this%V(this%cfg%imino,:,:)=this%RHOV(this%cfg%imino,:,:)/this%RHO(this%cfg%imino,:,:)
      !   this%W(this%cfg%imino,:,:)=this%RHOW(this%cfg%imino,:,:)/this%RHO(this%cfg%imino,:,:)
      !end if
      !if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
      !   this%U(:,this%cfg%jmino,:)=this%RHOU(:,this%cfg%jmino,:)/this%RHO(:,this%cfg%jmino,:)
      !   this%V(:,this%cfg%jmino,:)=this%RHOV(:,this%cfg%jmino,:)/this%RHO(:,this%cfg%jmino,:)
      !   this%W(:,this%cfg%jmino,:)=this%RHOW(:,this%cfg%jmino,:)/this%RHO(:,this%cfg%jmino,:)
      !end if
      !if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
      !   this%U(:,:,this%cfg%kmino)=this%RHOU(:,:,this%cfg%kmino)/this%RHO(:,:,this%cfg%kmino)
      !   this%V(:,:,this%cfg%kmino)=this%RHOV(:,:,this%cfg%kmino)/this%RHO(:,:,this%cfg%kmino)
      !   this%W(:,:,this%cfg%kmino)=this%RHOW(:,:,this%cfg%kmino)/this%RHO(:,:,this%cfg%kmino)
      !end if
   end subroutine get_velocity
   
   
   !> Calculate momentum from velocity
   subroutine get_momentum(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k
      ! Calculate momentum as far as possible
      !do k=this%cfg%kmino_+1,this%cfg%kmaxo_
      !   do j=this%cfg%jmino_+1,this%cfg%jmaxo_
      !      do i=this%cfg%imino_+1,this%cfg%imaxo_
      !         this%RHOU(i,j,k)=0.5_WP*(this%RHO(i-1,j,k)+this%RHO(i,j,k))*this%U(i,j,k)
      !         this%RHOV(i,j,k)=0.5_WP*(this%RHO(i,j-1,k)+this%RHO(i,j,k))*this%V(i,j,k)
      !         this%RHOW(i,j,k)=0.5_WP*(this%RHO(i,j,k-1)+this%RHO(i,j,k))*this%W(i,j,k)
      !      end do
      !   end do
      !end do
      ! Sync momentum
      !call this%cfg%sync(this%RHOU)
      !call this%cfg%sync(this%RHOV)
      !call this%cfg%sync(this%RHOW)
      ! Add last layer in each direction
      !if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then
      !   this%RHOU(this%cfg%imino,:,:)=this%RHO(this%cfg%imino,:,:)*this%U(this%cfg%imino,:,:)
      !   this%RHOV(this%cfg%imino,:,:)=this%RHO(this%cfg%imino,:,:)*this%V(this%cfg%imino,:,:)
      !   this%RHOW(this%cfg%imino,:,:)=this%RHO(this%cfg%imino,:,:)*this%W(this%cfg%imino,:,:)
      !end if
      !if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
      !   this%RHOU(:,this%cfg%jmino,:)=this%RHO(:,this%cfg%jmino,:)*this%U(:,this%cfg%jmino,:)
      !   this%RHOV(:,this%cfg%jmino,:)=this%RHO(:,this%cfg%jmino,:)*this%V(:,this%cfg%jmino,:)
      !   this%RHOW(:,this%cfg%jmino,:)=this%RHO(:,this%cfg%jmino,:)*this%W(:,this%cfg%jmino,:)
      !end if
      !if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
      !   this%RHOU(:,:,this%cfg%kmino)=this%RHO(:,:,this%cfg%kmino)*this%U(:,:,this%cfg%kmino)
      !   this%RHOV(:,:,this%cfg%kmino)=this%RHO(:,:,this%cfg%kmino)*this%V(:,:,this%cfg%kmino)
      !   this%RHOW(:,:,this%cfg%kmino)=this%RHO(:,:,this%cfg%kmino)*this%W(:,:,this%cfg%kmino)
      !end if
   end subroutine get_momentum
   
   
   !> Calculate the interpolated velocity, including overlap and ghosts
   subroutine interp_vel(this,Ui,Vi,Wi)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Ui !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Vi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Wi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      ! Calculate interpolated velocity as far as possible
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_,this%cfg%jmaxo_-1
            do i=this%cfg%imino_,this%cfg%imaxo_-1
               Ui(i,j,k)=0.5_WP*(this%U(i,j,k)+this%U(i+1,j,k))
               Vi(i,j,k)=0.5_WP*(this%V(i,j,k)+this%V(i,j+1,k))
               Wi(i,j,k)=0.5_WP*(this%W(i,j,k)+this%W(i,j,k+1))
            end do
         end do
      end do
      ! Sync interpolated velocity
      call this%cfg%sync(Ui)
      call this%cfg%sync(Vi)
      call this%cfg%sync(Wi)
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
         Ui(this%cfg%imaxo,:,:)=this%U(this%cfg%imaxo,:,:)
         Vi(this%cfg%imaxo,:,:)=this%V(this%cfg%imaxo,:,:)
         Wi(this%cfg%imaxo,:,:)=this%W(this%cfg%imaxo,:,:)
      end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
         Ui(:,this%cfg%jmaxo,:)=this%U(:,this%cfg%jmaxo,:)
         Vi(:,this%cfg%jmaxo,:)=this%V(:,this%cfg%jmaxo,:)
         Wi(:,this%cfg%jmaxo,:)=this%W(:,this%cfg%jmaxo,:)
      end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
         Ui(:,:,this%cfg%kmaxo)=this%U(:,:,this%cfg%kmaxo)
         Vi(:,:,this%cfg%kmaxo)=this%V(:,:,this%cfg%kmaxo)
         Wi(:,:,this%cfg%kmaxo)=this%W(:,:,this%cfg%kmaxo)
      end if
   end subroutine interp_vel
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,C,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: C !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(out) :: cfl
      integer :: ierr
      real(WP) :: maxvisc,maxC
      ! Compute convective+acoustic CFLs
      this%CFLc_x=maxval(abs(this%U)+abs(C))*dt*this%dxi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLc_y=maxval(abs(this%V)+abs(C))*dt*this%dyi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLc_z=maxval(abs(this%W)+abs(C))*dt*this%dzi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      ! Compute acoustic CFLs
      maxC=maxval(C); call MPI_ALLREDUCE(MPI_IN_PLACE,maxC,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLa_x=maxC*dt*this%dxi
      this%CFLa_y=maxC*dt*this%dyi
      this%CFLa_z=maxC*dt*this%dzi
      ! Compute viscous CFLs
      !maxvisc=maxval(this%visc/this%RHO); call MPI_ALLREDUCE(MPI_IN_PLACE,maxvisc,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLv_x=0.0_WP!4.0_WP*maxvisc*dt*this%dxi**2
      this%CFLv_y=0.0_WP!4.0_WP*maxvisc*dt*this%dyi**2
      this%CFLv_z=0.0_WP!4.0_WP*maxvisc*dt*this%dzi**2
      ! Return the maximum overall CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,&
      &       this%CFLa_x,this%CFLa_y,this%CFLa_z,&
      &       this%CFLv_x,this%CFLv_y,this%CFLv_z)
   end subroutine get_cfl
   
   
   !> Calculate info about our fields
   subroutine get_info(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k,ierr
      ! Compute integrals
      call this%cfg%integrate(this%VF  ,integral=this%VFint  )
      call this%cfg%integrate(this%RHOl,integral=this%RHOlint)
      !call this%cfg%integrate(this%RHOU,integral=this%RHOUint)
      !call this%cfg%integrate(this%RHOV,integral=this%RHOVint)
      !call this%cfg%integrate(this%RHOW,integral=this%RHOWint)
      !call this%cfg%integrate(this%RHOE,integral=this%RHOEint)
      ! Calculate extrema
      this%VFmin=+huge(1.0_WP)
      this%VFmax=-huge(1.0_WP)
      this%RHOlmin=+huge(1.0_WP)
      this%RHOlmax=-huge(1.0_WP)
      this%Umax=0.0_WP
      this%Vmax=0.0_WP
      this%Wmax=0.0_WP
      !this%Emin=+huge(1.0_WP)
      !this%Emax=-huge(1.0_WP)
      !this%Pmin=+huge(1.0_WP)
      !this%Pmax=-huge(1.0_WP)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%VFmin  =min(this%VFmin  ,    this%VF  (i,j,k) )
               this%VFmax  =max(this%VFmax  ,    this%VF  (i,j,k) )
               this%RHOlmin=min(this%RHOlmin,    this%RHOl(i,j,k) )
               this%RHOlmax=max(this%RHOlmax,    this%RHOl(i,j,k) )
               this%Umax  =max(this%Umax  ,abs(this%U  (i,j,k)))
               this%Vmax  =max(this%Vmax  ,abs(this%V  (i,j,k)))
               this%Wmax  =max(this%Wmax  ,abs(this%W  (i,j,k)))
               !this%Emin  =min(this%Emin  ,    this%E  (i,j,k) )
               !this%Emax  =max(this%Emax  ,    this%E  (i,j,k) )
               !this%Pmin  =min(this%Pmin  ,    this%P  (i,j,k) )
               !this%Pmax  =max(this%Pmax  ,    this%P  (i,j,k) )
            end do
         end do
      end do
      ! Get the parallel min/max
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%VFmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%VFmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOlmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOlmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax   ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmax   ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmax   ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      !call MPI_ALLREDUCE(MPI_IN_PLACE,this%Emin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      !call MPI_ALLREDUCE(MPI_IN_PLACE,this%Emax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      !call MPI_ALLREDUCE(MPI_IN_PLACE,this%Pmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      !call MPI_ALLREDUCE(MPI_IN_PLACE,this%Pmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
   end subroutine get_info
   
   
   !> Update a surfmesh object from our current polygons
   subroutine update_surfmesh(this,smesh)
      use surfmesh_class, only: surfmesh
      implicit none
      class(mpcomp),   intent(inout) :: this
      class(surfmesh), intent(inout) :: smesh
      integer :: i,j,k,n,shape,nv,np
      real(WP), dimension(3) :: tmp_vert
      ! Reset surface mesh storage
      call smesh%reset()
      ! First pass to count how many vertices and polygons are inside our processor
      nv=0; np=0
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         shape=getNumberOfVertices(this%interface_polygon(i,j,k))
         if (shape.gt.0) then
            nv=nv+shape
            np=np+1
         end if
      end do; end do; end do
      ! Reallocate storage and fill out arrays
      if (np.gt.0) then
         call smesh%set_size(nvert=nv,npoly=np)
         allocate(smesh%polyConn(smesh%nVert)) ! Also allocate naive connectivity
         nv=0; np=0
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            shape=getNumberOfVertices(this%interface_polygon(i,j,k))
            if (shape.gt.0) then
               ! Increment polygon counter
               np=np+1
               smesh%polySize(np)=shape
               ! Loop over its vertices and add them
               do n=1,shape
                  tmp_vert=getPt(this%interface_polygon(i,j,k),n-1)
                  ! Increment node counter
                  nv=nv+1
                  smesh%xVert(nv)=tmp_vert(1)
                  smesh%yVert(nv)=tmp_vert(2)
                  smesh%zVert(nv)=tmp_vert(3)
                  smesh%polyConn(nv)=nv
               end do
            end if
         end do; end do; end do
      else
         ! Add a zero-area triangle if this proc doesn't have one
         np=1; nv=3
         call smesh%set_size(nvert=nv,npoly=np)
         allocate(smesh%polyConn(smesh%nVert)) ! Also allocate naive connectivity
         smesh%xVert(1:3)=this%cfg%x(this%cfg%imin)
         smesh%yVert(1:3)=this%cfg%y(this%cfg%jmin)
         smesh%zVert(1:3)=this%cfg%z(this%cfg%kmin)
         smesh%polySize(1)=3
         smesh%polyConn(1:3)=[1,2,3]
      end if
   end subroutine update_surfmesh
   
   
   !> Print out info for mpcomp flow solver
   subroutine mpcomp_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(mpcomp), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("mpcomp solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine mpcomp_print
   
   
end module mpcomp_class
