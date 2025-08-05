!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs,VFlo,VFhi
   use tpscalar_class,    only: tpscalar,Lphase,Gphase
   use evap_class,        only: evap
   use string,            only: str_short,str_medium
   use YAMLRead,          only: YAMLElement
   use chem_sys_class,    only: chem_sys
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use mathtools,         only: Pi
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver, a volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: ss,vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(tpscalar),    public :: sc
   type(evap),        public :: evp
   type(timetracker), public :: time,timeSC

   !> The array of the species. Eeach stored as a YAMLElement object
   type(YAMLElement), dimension(:), allocatable :: species

   !> Species names
   character(len=str_medium), dimension(:), allocatable :: sp_names

   !> Chemical system
   type(chem_sys)   :: sys
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,scfile,evpfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resSC
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: T
   real(WP), dimension(:),       allocatable :: MM
   
   !> Problem definition
   real(WP) :: R0,R,R_ext,V_b
   real(WP), dimension(3) :: center
   integer  :: iTl,iTg,iO2=1,iWg=2,iWl=3,iN2=4
   real(WP) :: rho_l,rho_g,k_l,k_g,Cp_l,Cp_g,alpha_l,alpha_g,h_lg,T_sat
   real(WP) :: T_inf,beta,f_b_cnst,t0
   real(WP) :: pressure,air2vw_rat,N2O_rat
   integer :: ns
   ! Debug
   real(WP) :: prhs_int
   real(WP) :: mfr_err
   
contains


   !> Function that defines a level set function
   function levelset_bubble(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the bubble
      G=norm2(xyz-center)-R0
   end function levelset_bubble


   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x- boundary for scalar fields
   function xm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1) isIn=.true.
   end function xm_locator_sc


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   

   !> Function that localizes y- boundary
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator


   !> Function that localizes y- boundary for scalar fields
   function ym_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin-1) isIn=.true.
   end function ym_locator_sc
   
   
   !> Function that localizes y+ boundary
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator


   !> Function that localizes z- boundary
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator


   !> Function that localizes z- boundary for scalar fields
   function zm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin-1) isIn=.true.
   end function zm_locator_sc
   
   
   !> Function that localizes z+ boundary
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator
   

   ! Integrand function for beta equation
   function integrand(z,b)
      real(WP), intent(in) :: z,b
      real(WP) :: integrand
      integrand=exp(-b**2*((1.0_WP-z)**(-2)-2.0_WP*(1.0_WP-rho_g/rho_l)*z-1.0_WP))
   end function integrand


   ! Simpson's rule integration
   function Simpson(integrand,b,zmin,zmax,n)
      interface
      function integrand(x1,x2)
            use precision, only: WP
            real(WP), intent(in) :: x1,x2
            real(WP) :: integrand
         end function integrand
      end interface
      real(WP), intent(in) :: b
      real(WP), intent(in) :: zmin,zmax
      integer,  intent(in) :: n
      real(WP) :: Simpson
      integer  :: i
      real(WP) :: h
      h=(zmax-zmin)/real(n,WP)
      Simpson=integrand(zmin,b)-integrand(zmax,b)
      do i=0,n-2,2
         Simpson=Simpson+4.0_WP*integrand(zmin+(i+1)*h,b)+2.0_WP*integrand(zmin+(i+2)*h,b)
      end do
      Simpson=Simpson*h/3.0_WP
   end function Simpson


   !> Function that governs beta
   function f_beta(b)
      real(WP), intent(in) :: b
      real(WP) :: f_beta
      f_beta=2.0_WP*b**2*Simpson(integrand,b,0.0_WP,1.0_WP,20)-f_b_cnst
   end function f_beta


   !> Function that returns the analytical bubble radius
   function get_Rext(tim)
      real(WP), intent(in) :: tim
      real(WP) :: get_Rext
      get_Rext=2.0_WP*beta*sqrt(alpha_l*tim)
   end function get_Rext


   !> Function that returns the numerical bubble radius
   function get_R()
      use mathtools, only: Pi
      real(WP) :: get_R
      call sc%cfg%integrate(sc%PVF(:,:,:,Gphase),V_b)
      get_R=(3.0_WP*V_b/(4.0_WP*Pi))**(1.0_WP/3.0_WP)
   end function get_R

   
   !> Initialization of problem solver
   subroutine simulation_init
      use param,    only: param_read,param_getsize
      use messager, only: die
      implicit none
      integer :: np=2,ne,ncs
      character(len=str_short), dimension(:), allocatable :: e_names
      real(WP), dimension(:,:), allocatable :: elem_mat
      real(WP), dimension(:,:), allocatable :: phse_mat
      real(WP), allocatable :: nasa_coef(:,:)
      ! real(WP), dimension(:), allocatable :: N_init
      character(len=str_medium), dimension(:), allocatable :: const_sp
      integer,  dimension(:), allocatable :: CS
      

      ! Parse the mechanism file
      parse_mech: block
         use chem_sys_class, only: ncof
         use YAMLRead,       only: YAMLHandler,YAMLSequence,YAMLMap,yaml_open_file,yaml_start_from_sequence,yaml_close_file
         character(len=str_medium) :: mch_file
         character(len=str_short), dimension(:), allocatable :: sp_names_copy,const_sp_copy
         ! real(WP), dimension(:), allocatable :: N_init_copy
         type(YAMLHandler)  :: domain
         type(YAMLSequence) :: sp_list,phases,elements
         type(YAMLElement)  :: sp,gas
         type(YAMLMap)      :: thermo,comp
         integer :: isc,nn,i,j,k,e,code
         character(len=:), allocatable :: name_arr(:)
         character(len=:), allocatable :: name
         real(WP), allocatable :: T_range(:)
         real(WP), dimension(:,:), allocatable :: a
         logical :: new_elem
         ! Get the target species from input
         ns=param_getsize('Species')
         ncs=param_getsize('Constrained species')
         allocate(sp_names(1:ns))
         allocate(sp_names_copy(1:ns))
         ! allocate(N_init(1:ns))
         ! allocate(N_init_copy(1:ns))
         allocate(const_sp(1:ncs))
         allocate(const_sp_copy(1:ncs))
         allocate(CS(ncs))
         call param_read('Species',sp_names)
         ! call param_read('Initial moles',N_init)
         call param_read('Constrained species',const_sp)
         sp_names_copy=sp_names
         ! N_init_copy=N_init
         const_sp_copy=const_sp
         ! Read the mechanism file path
         call param_read('Mechanism file',mch_file)
         ! Open the mechanism
         domain=yaml_open_file(trim(mch_file))
         ! Get the list of all species
         sp_list=yaml_start_from_sequence(domain,'species')
         ! Extract the target species from the mechanism
         allocate(species(1:ns))
         nn=0
         k=0
         do isc=0,sp_list%size-1 ! Index in YAMLSequence starts from 0
            sp=sp_list%element(isc)
            name_arr=sp%value_str('name',code)
            name=''
            do i=1,size(name_arr)
               name=trim(name//name_arr(i))
            end do
            do i=1,ns
               if (sp_names_copy(i).eq.name) then
                  nn=nn+1
                  species(nn)=sp
                  sp_names(nn)=name
                  ! N_init(nn)=N_init_copy(i)
                  do j=1,ncs
                     if (const_sp_copy(j).eq.name) then
                        k=k+1
                        const_sp(k)=name
                        CS(k)=nn
                     end if
                  end do
               end if
            end do
            call sp%destroy()
         end do
         if(nn.ne.ns) call die('Some species are missing in the mechanism file.')
         ! Get the elements that exist in the target species
         phases=yaml_start_from_sequence(domain,'phases')
         gas=phases%element(0)
         elements=gas%value_sequence('elements',code) ! For some reason I couldn't directly get the element names using yaml-fortran
         allocate(e_names(elements%size)); e_names=''
         ne=0
         do isc=1,ns
            sp=species(isc)
            comp=sp%value_map('composition')
            do e=1,size(comp%labels)
               name=''
               do i=1,size(comp%labels(e)%str)
                  name=name//trim(comp%labels(e)%str(i))
               end do
               new_elem=.true.
               do i=1,ne
                  if (trim(e_names(i)).eq.trim(name)) then
                     new_elem=.false.
                     exit
                  end if
               end do
               if (new_elem) then
                  ne=ne+1
                  e_names(ne)=trim(name)
               end if
            end do
         end do
         e_names=e_names(1:ne)
         ! Form the element matrix
         allocate(elem_mat(ns,ne)); elem_mat=0.0_WP
         do isc=1,ns
            sp=species(isc)
            comp=sp%value_map('composition')
            do i=1,size(comp%labels)
               name=''
               do nn=1,size(comp%labels(i)%str)
                  name=name//trim(comp%labels(i)%str(nn))
               end do
               do e=1,ne
                  if (trim(e_names(e)).eq.trim(name)) elem_mat(isc,e)=real(comp%value_int(name,code),WP)
               end do
            end do
         end do
         ! Read the NASA-7 polynomials
         allocate(nasa_coef(1:ns,2*ncof+1)); nasa_coef=0.0_WP
         allocate(a(1:2,1:ncof)); a=0.0_WP
         do isc=1,ns
            sp=species(isc)
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
            nasa_coef(isc,1)=T_range(2)
            nasa_coef(isc,2:  ncof+1)=a(1,:)
            nasa_coef(isc,9:2*ncof+1)=a(2,:)
         end do
         ! Form the phase summation matrix
         allocate(phse_mat(ns,Lphase+1:Gphase+1)); phse_mat(:,Lphase+1)=0.0_WP; phse_mat(:,Gphase+1)=1.0_WP
         do isc=1,ns
            if (len_trim(sp_names(isc)).ge.3) then
               if (sp_names(isc)(len_trim(sp_names(isc))-2:len_trim(sp_names(isc))).eq.'(L)') then
                  phse_mat(isc,Lphase+1)=1.0_WP
                  phse_mat(isc,Gphase+1)=0.0_WP
               end if
            end if
         end do
         ! Close the mechanism file and clean up
         call yaml_close_file(domain)
         call sp_list%destroy()
         call sp%destroy()
         call comp%destroy()
         call thermo%destroy()
         ! deallocate(sp_names_copy,N_init_copy,const_sp_copy)
         deallocate(sp_names_copy,const_sp_copy)
      end block parse_mech

      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resSC (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:2))
         allocate(resU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(T     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(MM(ns)); MM=[32.0_WP,18.0_WP,18.0_WP,28.0_WP]; MM=0.001_WP*MM
      end block allocate_work_arrays


      ! Initialize the chemical system
      sys_init: block
         use messager, only: die
         integer :: ng=1
         real(WP), dimension(:,:), allocatable :: Bg
         character(len=2) :: eq_cond
         integer :: isc
         ! Allocate arrays
         allocate(Bg(ns,ng));  Bg=0.0_WP
         ! Create the general constraints
         do isc=1,ns
            if (sp_names(isc).eq.'H2O')    Bg(isc,1)=1.0_WP
            if (sp_names(isc).eq.'H2O(L)') Bg(isc,1)=1.0_WP
         end do
         ! Chemical system object
         call sys%initialize(np=np,ns=ns,ne=ne,ncs=ncs,ng=ng,P=phse_mat,Ein=elem_mat,CS=CS,Bg=Bg,thermo_in=nasa_coef,diag=5)
         ! Deallocate arrays
         deallocate(Bg)
      end block sys_init
      
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='Main')
         call param_read('Initial time',t0)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Sub-iterations',time%itmax)
         ! time%t=t0
         time%dt=time%dtmax
         timeSC=timetracker(amRoot=cfg%amRoot,name='SC Time')
         call param_read('Scalar time step',timeSC%dtmax)
         timeSC%tmax=time%tmax
         timeSC%dt=timeSC%dtmax
      end block initialize_timetracker


      ! Read problem inputs
      read_inputs: block
         call param_read('Liquid density',rho_l)
         call param_read('Gas density',rho_g)
         call param_read('Saturation temperature',T_sat)
         call param_read('Latent heat',h_lg)
         call param_read('Liquid thermal conductivity',k_l)
         call param_read('Liquid specific heat capacity',Cp_l)
         call param_read('Gas thermal conductivity',k_g)
         call param_read('Gas specific heat capacity',Cp_g)
         call param_read('Far-field temperature',T_inf)
         call param_read('Pressure',Pressure)
         alpha_l=k_l/(rho_l*Cp_l)
         alpha_g=k_g/(rho_g*Cp_g)
         f_b_cnst=(rho_l*Cp_l*(T_inf-T_sat))/(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))
      end block read_inputs


      ! Analytical solution
      analytical_solution: block
         use mpi_f08,  only: MPI_BCAST
         use parallel, only: MPI_REAL_WP
         real(WP) :: err,tol
         real(WP) :: beta0,betaL,betaR,fL,fR,fM
         integer :: it,itmax,ierr
         logical :: convergence
         if (cfg%amRoot) then
            convergence=.false.
            itmax=200
            tol=1e-9
            betaL=0.0_WP
            betaR=5.0_WP
            do it=1,itmax
               beta=0.5_WP*(betaL+betaR)
               fM=f_beta(beta)
               err=abs(fM)
               if (err.lt.tol) then
                  convergence=.true.
                  exit
               end if
               fL=f_beta(betaL)
               fR=f_beta(betaR)
               if (fL*fM.lt.0.0_WP) then
                  betaR=beta
               else
                  betaL=beta
               endif
            end do
            print*,'Convergence = ',convergence
            print*,'beta =',beta
            print*,'Absolute error =',err
            print*,'Bi-section iterations =',it
            ! Debug
            print*,'Analytical mass flux = ',2.0_WP*k_l*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)/(get_Rext(t0)*h_lg)*integrand(0.0_WP,beta)
         end if
         call MPI_BCAST(beta,1,MPI_REAL_WP,0,cfg%comm,ierr)
      end block analytical_solution
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,flux_storage,neumann
         integer :: i,j,k,n,si,sj,sk,ierr
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux_storage,nband=6,name='VOF')
         ! Boundary conditinos
         call vf%add_bcond(name='xm',type=neumann,locator=xm_locator_sc,dir='-x')
         call vf%add_bcond(name='xp',type=neumann,locator=xp_locator   ,dir='+x')
         call vf%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         call vf%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         call vf%add_bcond(name='zm',type=neumann,locator=zm_locator_sc,dir='-z')
         call vf%add_bcond(name='zp',type=neumann,locator=zp_locator   ,dir='+z')
         ! Initialize the VOF field
         call param_read('Bubble center',center)
         R0=get_Rext(t0)
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_bubble,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Apply boundary conditions
         call vf%apply_bcond(time%t,time%dt)
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Apply symmetry condition
         ! symmetry_irl: block
         !    use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
         !    &                                setNumberOfPlanes,setPlane,matchVolumeFraction
         !    real(WP), dimension(1:4) :: plane
         !    type(RectCub_type) :: cell
         !    integer :: ii,jj,kk
         !    call new(cell)
         !    if (vf%cfg%iproc.eq.1) then
         !       do k=vf%cfg%kmino_,vf%cfg%kmaxo_
         !          do j=vf%cfg%jmino_,vf%cfg%jmaxo_
         !             do i=vf%cfg%imino_,vf%cfg%imin_-1
         !                ii=(vf%cfg%imin_-i)+vf%cfg%imin_-1
         !                plane=getPlane(vf%liquid_gas_interface(ii,j,k),0)
         !                call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
         !                &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
         !                plane(4)=dot_product([-plane(1),plane(2),plane(3)],[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
         !                call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
         !                call setPlane(vf%liquid_gas_interface(i,j,k),0,[-plane(1),plane(2),plane(3)],plane(4))
         !                call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
         !             end do
         !          end do
         !       end do
         !    end if
         !    if (vf%cfg%jproc.eq.1) then
         !       do k=vf%cfg%kmino_,vf%cfg%kmaxo_
         !          do j=vf%cfg%jmino_,vf%cfg%jmin_-1
         !             jj=(vf%cfg%jmin_-j)+vf%cfg%jmin_-1
         !             do i=vf%cfg%imino_,vf%cfg%imaxo_
         !                plane=getPlane(vf%liquid_gas_interface(i,jj,k),0)
         !                call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
         !                &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
         !                plane(4)=dot_product([plane(1),-plane(2),plane(3)],[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
         !                call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
         !                call setPlane(vf%liquid_gas_interface(i,j,k),0,[plane(1),-plane(2),plane(3)],plane(4))
         !                call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
         !             end do
         !          end do
         !       end do
         !    end if
         !    if (vf%cfg%kproc.eq.1) then
         !       do k=vf%cfg%kmino_,vf%cfg%kmin_-1
         !          kk=(vf%cfg%kmin_-k)+vf%cfg%kmin_-1
         !          do j=vf%cfg%jmino_,vf%cfg%jmaxo_
         !             do i=vf%cfg%imino_,vf%cfg%imaxo_
         !                plane=getPlane(vf%liquid_gas_interface(i,j,kk),0)
         !                call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
         !                &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
         !                plane(4)=dot_product([plane(1),plane(2),-plane(3)],[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
         !                call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
         !                call setPlane(vf%liquid_gas_interface(i,j,k),0,[plane(1),plane(2),-plane(3)],plane(4))
         !                call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
         !             end do
         !          end do
         !       end do
         !    end if
         ! end block symmetry_irl
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof


      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: gmres_pfmg2
         use mathtools,       only: Pi
         use tpns_class,      only: clipped_neumann,bcond
         type(bcond), pointer :: mybc
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         fs%rho_l=rho_l
         fs%rho_g=rho_g
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         fs%contact_angle=fs%contact_angle*Pi/180.0_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Boundary conditions
         call fs%add_bcond(name='xm',type=clipped_neumann,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
         call fs%add_bcond(name='xp',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         call fs%add_bcond(name='ym',type=clipped_neumann,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call fs%add_bcond(name='yp',type=clipped_neumann,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call fs%add_bcond(name='zm',type=clipped_neumann,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         call fs%add_bcond(name='zp',type=clipped_neumann,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=gmres_pfmg2,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         call param_read('Max coarsening levels',ps%maxlevel)
         ! Implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Apply boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      
      ! Create a one-sided scalar solver
      create_scalar: block
         use param, only: param_read
         use tpscalar_class,  only: bcond,neumann,dirichlet
         type(bcond), pointer :: mybc
         real(WP) :: radius,mp(2),N_init(ns)
         integer  :: i,j,k,isc,p
         integer  :: pos_open,pos_close
         ! Read-in inputs
         call param_read('Air to water vapor mole ratio',air2vw_rat)
         call param_read('Nitrogen to oxygen mole ratio',N2O_rat)
         ! Create scalar solver
         call sc%initialize(cfg=cfg,nscalar=ns+2,name='tpscalar')
         ! Boundary conditinos
         call sc%add_bcond(name='xm',type=neumann,locator=xm_locator_sc,dir='-x')
         call sc%add_bcond(name='xp',type=neumann,locator=xp_locator   ,dir='+x')
         call sc%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         call sc%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         call sc%add_bcond(name='zm',type=neumann,locator=zm_locator_sc,dir='-z')
         call sc%add_bcond(name='zp',type=neumann,locator=zp_locator   ,dir='+z')
         sc%SCname=[sp_names,'Tl','Tg']; iTl=ns+1; iTg=ns+2
         do isc=1,ns
            pos_open =index(sc%SCname(isc),'(')
            pos_close=index(sc%SCname(isc),')')
            if (pos_open > 0 .and. pos_close > pos_open) then
               sc%SCname(isc)(pos_open:pos_open)   = '_'
               sc%SCname(isc)(pos_close:pos_close) = ' '
               sc%SCname(isc)=adjustl(trim(sc%SCname(isc)))
            end if
            sc%phase(isc)=sys%get_pind(isc)
         end do
         sc%phase(iTl)=Lphase
         sc%phase(iTg)=Gphase
         sc%diff(:,:,:,iTl)=alpha_l
         sc%diff(:,:,:,iTg)=alpha_g
         ! Initialize the linear solver
         ss=ddadi(cfg=cfg,name='Scalar',nst=7)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
         ! Initialize the phasic density and VOF
         sc%Prho(Lphase)=fs%rho_l
         sc%Prho(Gphase)=fs%rho_g
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         ! Initialize scalar fields
         N_init(iWl)=1.0_WP
         N_init(iWg)=1.0_WP
         N_init(iO2)=air2vw_rat*N_init(iWg)
         N_init(iN2)=N2O_rat*N_init(iO2)
         mp=0.0_WP
         do isc=1,ns
            p=sc%phase(isc)
            mp(p)=mp(p)+N_init(isc)*MM(isc)
         end do
         do isc=1,ns
            p=sc%phase(isc)
            where (sc%PVF(:,:,:,p).gt.VFlo)
               sc%SC(:,:,:,isc)=MM(isc)*N_init(isc)/mp(p)
            else where
               sc%SC(:,:,:,isc)=0.0_WP
            end where
         end do
         do i=sc%cfg%imino_,sc%cfg%imaxo_
            do j=sc%cfg%jmino_,sc%cfg%jmaxo_
               do k=sc%cfg%kmino_,sc%cfg%kmaxo_
                  radius=norm2([sc%cfg%xm(i),sc%cfg%ym(j),sc%cfg%zm(k)])
                  if (vf%VF(i,j,k).gt.VFlo) then
                     sc%SC(i,j,k,iTl)=T_inf-2.0_WP*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)*Simpson(integrand,beta,1.0_WP-R0/radius,1.0_WP,20)
                  end if
                  if (vf%VF(i,j,k).lt.VFhi) then
                     sc%SC(i,j,k,iTg)=T_sat
                  end if
               end do
            end do
         end do
         where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
            sc%SC(:,:,:,iTl)=T_sat
            sc%SC(:,:,:,iTg)=T_sat
         end where
         call sc%apply_bcond(time%t,time%dt)
         ! Post process
         T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)
      end block create_scalar
      

      ! Create and initialize an evp object
      create_evp: block
         ! use evap_class, only: symmetry
         integer :: i,j,k
         ! Create the object
         call evp%initialize(cfg=cfg,vf=vf,sc=sc%SC,iTl=iTl,iTg=iTg,itp_x=fs%itpr_x,itp_y=fs%itpr_y,itp_z=fs%itpr_z,div_x=fs%divp_x,div_y=fs%divp_y,div_z=fs%divp_z,name='liquid gas pc')
         call param_read('Mass flux tolerence',     evp%mflux_tol)
         call param_read('Max pseudo timestep size',evp%pseudo_time%dtmax)
         call param_read('Max pseudo cfl number',   evp%pseudo_time%cflmax)
         call param_read('Max pseudo time steps',   evp%pseudo_time%nmax)
         evp%pseudo_time%dt=evp%pseudo_time%dtmax
         ! Boundary conditions
         ! call evp%add_bcond(name='xm',type=symmetry,face='x',dir=-1,locator=xm_locator_sc)
         ! call evp%add_bcond(name='ym',type=symmetry,face='y',dir=-1,locator=ym_locator_sc)
         ! call evp%add_bcond(name='zm',type=symmetry,face='z',dir=-1,locator=zm_locator_sc)
         ! Get densities from the flow solver
         evp%rho_l=fs%rho_l
         evp%rho_g=fs%rho_g
      end block create_evp


      ! Interface jump conditions
      interface_jump: block
         ! Debug
         use mpi_f08,          only: MPI_MAX,MPI_ALLREDUCE
         use parallel,         only: MPI_REAL_WP
         use chem_state_class, only: chem_state,fixed_PT,fixed_PH
         type(chem_state) :: state
         real(WP), dimension(:), allocatable :: vol,mp,N
         real(WP) :: V
         integer :: i,j,k,index,isc,p
         ! Debug
         real(WP) :: my_mdotdp_max,mdotdp_max,m_old
         real(WP) :: T_g
         integer :: ierr
         ! Allocate arrays
         allocate(vol(Lphase:Gphase))
         allocate(mp(Lphase:Gphase))
         allocate(N(ns))
         ! Initialize the chemical state
         ! call state%initialize(sys=sys,cond=fixed_PT,p=pressure)
         call state%initialize(sys=sys,cond=fixed_PH,p=pressure)
         call param_read('Newton tolerance',state%tol_N)
         call param_read('Newton max iterations',state%iter_N_max)
         if (state%cond.eq.fixed_PH) then
            call param_read('T tolerance',state%tol_T)
            call param_read('T max iterations',state%iter_T_max)
            call param_read('Temperature initial guess',T_g)
         end if
         ! Loop over the interfacial cells
         do index=1,vf%band_count(0)
            print*,''
            print*,''
            print*,'------------------------------------------------------------------------------------------------------------------------'
            ! Get the interfacial cell indices
            i=vf%band_map(1,index)
            j=vf%band_map(2,index)
            k=vf%band_map(3,index)
            print*,'VOF(',i,',',j,',',k,') = ',vf%VF(i,j,k)
            ! Calculate the phase masses
            mp=sc%Prho*sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
            m_old=sum(mp)
            print*,'Phase mass in the cell prior to equilibrium = ',mp
            print*,'Total mass in the cell prior to equilibrium = ',m_old
            ! Calculate the moles numbers
            do isc=1,ns
               p=sc%phase(isc)
               N(isc)=sc%SC(i,j,k,isc)*mp(p)/MM(isc)
            end do
            print*,'N prior to equilibrium: '
            do isc=1,ns
               print*,sc%SCname(isc),': ',N(isc)
            end do
            print*,'Total moles prior to equilibrium = ',sum(N)
            ! Get the chemical equilibrium
            ! call state%N_init(T=T_sat,N=N)
            call state%N_init(N=N,N_h=N,T_h=T_sat,T_g=T_g)
            print*,'N reinitialized = ',state%N
            print*,'Total moles prior to equilibrium after reinitialization = ',sum(N)
            print*,''
            call state%equilibrate()
            print*,'Relative change in total moles of the system = ',(sum(state%N)-sum(N))/sum(N)
            N=state%N
            print*,'Equilibrium N = ',N
            print*,'Total equilibrium moles = ',sum(N)
            ! Update the phase masses
            mp=0.0_WP
            do isc=1,ns
               mp(sc%phase(isc))=mp(sc%phase(isc))+N(isc)*MM(isc)
            end do
            print*,'Equilibrium mass of phases = ',mp
            print*,'Equilibrium mass of system = ',sum(mp)
            print*,'Relative change in the total mass of system = ',(sum(mp)-m_old)/m_old
            ! Get the phase volumes
            vol=0.0_WP
            vol=mp/sc%Prho
            V=sum(vol)
            print*,'cell vol = ',cfg%vol(i,j,k)
            print*,'Equilibrium phase volumes = ',vol
            print*,'Equilibrium system volume = ',V
            ! Get the phase change mass flux
            evp%mdotdp(i,j,k)=(V-cfg%vol(i,j,k))/(time%dt*(1.0_WP/sc%Prho(Gphase)-1.0_WP/sc%Prho(Lphase))*cfg%vol(i,j,k)*vf%SD(i,j,k))
            print*,'evp%mdotdp(i,j,k) = ',evp%mdotdp(i,j,k)
            ! Relax the interface temperature
            sc%SC(i,j,k,iTl)=state%T
            sc%SC(i,j,k,iTg)=state%T
            T(i,j,k)=state%T
            print*,'Temperature = ',state%T
            ! Relax the interface VOF (Need to update the corresponding IRL quantities in vf?)
            vf%VF(i,j,k)=vol(Lphase)/V
            sc%PVF(i,j,k,Lphase)=vf%VF(i,j,k)
            sc%PVF(i,j,k,Gphase)=1.0_WP-vf%VF(i,j,k)
            ! Relax the interface composition
            do isc=1,ns
               sc%SC(i,j,k,isc)=MM(isc)*N(isc)/mp(sc%phase(isc))
            end do
         end do
         ! Correct Y_H2O(G) inside the bubble (Assuming homogenous gas mixture at initial time)
         ! where (vf%VF.le.VFlo) sc%SC(:,:,:,iWg)=MM(iWg)*N(iWg)/mp(sc%phase(iWg))
         ! Sync fields
         do isc=1,sc%nscalar
            call cfg%sync(sc%SC(:,:,:,isc))
         end do
         call cfg%sync(vf%VF)
         call cfg%sync(sc%PVF(:,:,:,Lphase))
         call cfg%sync(sc%PVF(:,:,:,Gphase))
         call cfg%sync(evp%mdotdp)
         ! Apply boundary conditions
         call sc%apply_bcond(time%t,time%dt)
         call vf%apply_bcond(time%t,time%dt)
         ! Get the volumetric evaporation mass flux
         call evp%get_mflux()
         ! Initialize the liquid and gas side mass fluxes
         call evp%init_mfluxLG()
         ! Get the interface normal
         call evp%get_normal()
         ! Debug
         my_mdotdp_max=maxval(evp%mdotdp)
         call MPI_ALLREDUCE(my_mdotdp_max,mdotdp_max,1,MPI_REAL_WP,MPI_MAX,evp%cfg%comm,ierr)
         if (evp%cfg%amRoot) print*,'Numerical mass flux = ',mdotdp_max
         ! Deallocate arrays
         deallocate(vol,mp,N)
      end block interface_jump


      ! Initialize bubble radius
      R=get_R()
      R_ext=get_Rext(time%t+t0)


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         integer :: isc
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='Bubble_vapor')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_surface('plic',smesh)
         do isc=1,ns
           call ens_out%add_scalar('Y_'//trim(sc%SCname(isc)),sc%SC(:,:,:,isc))
         end do
         do isc=ns+1,sc%nscalar
           call ens_out%add_scalar(trim(sc%SCname(isc)),sc%SC(:,:,:,isc))
         end do
         call ens_out%add_scalar('mflux',evp%mflux)
         call ens_out%add_scalar('evp_div',evp%div_src)
         call ens_out%add_scalar('mdotdp',evp%mdotdp)
         call ens_out%add_scalar('mfluxL',evp%mfluxLG(:,:,:,Lphase))
         call ens_out%add_scalar('mfluxG',evp%mfluxLG(:,:,:,Gphase))
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('Temperature',T)
         call ens_out%add_vector('normal',evp%normal(:,:,:,1),evp%normal(:,:,:,2),evp%normal(:,:,:,3))
         call ens_out%add_scalar('Tl_grd',evp%Tl_grd)
         call ens_out%add_scalar('Tg_grd',evp%Tg_grd)
         ! Debug
         call ens_out%add_scalar('PVFG',sc%PVF(:,:,:,Gphase))
         call ens_out%add_scalar('PVFL',sc%PVF(:,:,:,Lphase))
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         integer :: isc
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call sc%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(R,'R')
         call mfile%add_column(R_ext,'R_ext')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%add_column(prhs_int,'prhs_int')
         ! Debug
         call fs%get_mfr()
         call evp%cfg%integrate(evp%div_src,mfr_err)
         mfr_err=abs(mfr_err-sum(fs%mfr))
         call mfile%add_column(mfr_err,'mfr_err')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create scalar monitor
         scfile=monitor(sc%cfg%amRoot,'scalar')
         call scfile%add_column(time%n,'Timestep number')
         call scfile%add_column(time%t,'Time')
         do isc=1,sc%nscalar
           call scfile%add_column(sc%SCmin(isc),trim(sc%SCname(isc))//'_min')
           call scfile%add_column(sc%SCmax(isc),trim(sc%SCname(isc))//'_max')
           call scfile%add_column(sc%SCint(isc),trim(sc%SCname(isc))//'_int')
         end do
         call scfile%write()
         ! Create evaporation monitor
         evpfile=monitor(evp%cfg%amRoot,'evaporation')
         call evpfile%add_column(time%n,'Timestep number')
         call evpfile%add_column(time%t,'Time')
         call evpfile%add_column(evp%pseudo_time%dt,'Pseudo time step')
         call evpfile%add_column(evp%pseudo_time%cfl,'Maximum pseudo CFL')
         call evpfile%add_column(evp%pseudo_time%n,'No. pseudo steps')
         call evpfile%add_column(evp%mflux_int,'mflux int')
         call evpfile%add_column(evp%mfluxL_int,'shifted mfluxL int')
         call evpfile%add_column(evp%mfluxG_int,'shifted mfluxG int')
         call evpfile%add_column(evp%mfluxL_int_err,'mfluxL int err')
         call evpfile%add_column(evp%mfluxG_int_err,'mfluxG int err')
         call evpfile%add_column(evp%mfluxL_err,'max mfluxL err')
         call evpfile%add_column(evp%mfluxG_err,'max mfluxG err')
         call evpfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation-this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use tpns_class, only: static_contact,harmonic_visc
      use mathtools,  only: Pi
      implicit none
      integer  :: i,j,k

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old VOF
         vf%VFold=vf%VF

         ! Remember old evaporation divergence
         evp%div_src_old=evp%div_src
         
         ! Remember old SC
         sc%SCold =sc%SC
         sc%PVFold=sc%PVF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)

         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         call vf%apply_bcond(time%t,time%dt)

         ! Apply symmetry condition
         ! symmetry_irl: block
         !    use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
         !    &                                setNumberOfPlanes,setPlane,matchVolumeFraction
         !    real(WP), dimension(1:4) :: plane
         !    type(RectCub_type) :: cell
         !    integer :: ii,jj,kk
         !    call new(cell)
         !    if (vf%cfg%iproc.eq.1) then
         !       do k=vf%cfg%kmino_,vf%cfg%kmaxo_
         !          do j=vf%cfg%jmino_,vf%cfg%jmaxo_
         !             do i=vf%cfg%imino_,vf%cfg%imin_-1
         !                ii=(vf%cfg%imin_-i)+vf%cfg%imin_-1
         !                plane=getPlane(vf%liquid_gas_interface(ii,j,k),0)
         !                call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
         !                &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
         !                plane(4)=dot_product([-plane(1),plane(2),plane(3)],[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
         !                call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
         !                call setPlane(vf%liquid_gas_interface(i,j,k),0,[-plane(1),plane(2),plane(3)],plane(4))
         !                call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
         !             end do
         !          end do
         !       end do
         !    end if
         !    if (vf%cfg%jproc.eq.1) then
         !       do k=vf%cfg%kmino_,vf%cfg%kmaxo_
         !          do j=vf%cfg%jmino_,vf%cfg%jmin_-1
         !             jj=(vf%cfg%jmin_-j)+vf%cfg%jmin_-1
         !             do i=vf%cfg%imino_,vf%cfg%imaxo_
         !                plane=getPlane(vf%liquid_gas_interface(i,jj,k),0)
         !                call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
         !                &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
         !                plane(4)=dot_product([plane(1),-plane(2),plane(3)],[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
         !                call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
         !                call setPlane(vf%liquid_gas_interface(i,j,k),0,[plane(1),-plane(2),plane(3)],plane(4))
         !                call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
         !             end do
         !          end do
         !       end do
         !    end if
         !    if (vf%cfg%kproc.eq.1) then
         !       do k=vf%cfg%kmino_,vf%cfg%kmin_-1
         !          kk=(vf%cfg%kmin_-k)+vf%cfg%kmin_-1
         !          do j=vf%cfg%jmino_,vf%cfg%jmaxo_
         !             do i=vf%cfg%imino_,vf%cfg%imaxo_
         !                plane=getPlane(vf%liquid_gas_interface(i,j,kk),0)
         !                call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
         !                &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
         !                plane(4)=dot_product([plane(1),plane(2),-plane(3)],[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
         !                call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
         !                call setPlane(vf%liquid_gas_interface(i,j,k),0,[plane(1),plane(2),-plane(3)],plane(4))
         !                call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
         !             end do
         !          end do
         !       end do
         !    end if
         ! end block symmetry_irl

         ! Transport scalars
         advance_scalar: block
            integer  :: isc
            real(WP) :: dt_sc

            ! Update the phas-specific VOF
            sc%PVF(:,:,:,Lphase)=vf%VF
            sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF

            ! Update the phas-specific face apertures
            call sc%get_face_apt()

            ! Explicit calculation of dVOFSC/dt from scalar advection
            call sc%get_dSCdt_adv(dSCdt=resSC,U=fs%U,V=fs%V,W=fs%W,detailed_face_flux=vf%detailed_face_flux,dt=time%dt)

            ! Advance scalar advection
            do isc=1,sc%nscalar
               where (sc%mask.eq.0.and.sc%PVF(:,:,:,sc%phase(isc)).gt.VFlo) sc%SC(:,:,:,isc)=(sc%PVFold(:,:,:,sc%phase(isc))*sc%SCold(:,:,:,isc)+time%dt*(resSC(:,:,:,isc)+evp%div_src_old(:,:,:)*sc%SCold(:,:,:,isc)))/sc%PVF(:,:,:,sc%phase(isc))
               where (sc%PVF(:,:,:,sc%phase(isc)).eq.0.0_WP) sc%SC(:,:,:,isc)=0.0_WP
            end do
            ! where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
            !    sc%SC(:,:,:,iTl)=T_sat
            !    sc%SC(:,:,:,iTg)=T_sat
            ! end where

            ! Advance scalar diffusion
            do while (timeSC%t.lt.time%t)

               if (timeSC%t+timeSC%dt.gt.time%t) then
                  dt_sc=timeSC%dt
                  timeSC%dt=time%t-timeSC%t
                  call timeSC%increment()
                  timeSC%dt=dt_sc
               else
                  call timeSC%increment()
               end if
               sc%SCold=sc%SC

               ! Explicit calculation of dVOFSC/dt from scalar diffusion
               call sc%get_dSCdt_dff(dSCdt=resSC)
               do isc=1,sc%nscalar
                  where (sc%mask.eq.0.and.sc%PVF(:,:,:,sc%phase(isc)).gt.VFlo) resSC(:,:,:,isc)=timeSC%dt*resSC(:,:,:,isc)/sc%PVF(:,:,:,sc%phase(isc))
                  where (sc%PVF(:,:,:,sc%phase(isc)).eq.0.0_WP) resSC(:,:,:,isc)=0.0_WP
               end do

               ! Form implicit diffusive residuals
               call sc%solve_implicit_diff(timeSC%dt,resSC)

               ! Apply the residuals
               sc%SC=sc%SC+resSC

               ! Apply boundary conditions
               ! where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
               !    sc%SC(:,:,:,iTl)=T_sat
               !    sc%SC(:,:,:,iTg)=T_sat
               ! end where
               call sc%apply_bcond(timeSC%t,timeSC%dt)

            end do

            ! Update the one-field temperature
            T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)

         end block advance_scalar

         ! Interface jump conditions
         interface_jump: block
            use chem_state_class, only: chem_state,fixed_PH
            type(chem_state) :: state
            real(WP), dimension(:), allocatable :: N_h,vol
            real(WP) :: T_h,V
            integer :: index,isc,p
            ! Allocate arrays
            allocate(N_h(ns))
            allocate(vol(Lphase:Gphase))
            ! Initialize the chemical state
            call state%initialize(sys=sys,cond=fixed_PH,p=pressure)
            ! Loop over the interfacial cells
            do index=1,vf%band_count(0)
               ! Get the interfacial cell indices
               i=vf%band_map(1,index)
               j=vf%band_map(2,index)
               k=vf%band_map(3,index)
               ! Calculate the mole numbers
               do isc=1,ns
                  p=sc%phase(isc)
                  N_h(isc)=sc%Prho(p)*cfg%vol(i,j,k)/MM(isc)*sc%PVF(i,j,k,p)*sc%SC(i,j,k,isc)
               end do
               ! Initialize the chemical state (what should be T_h?)
               T_h=T_sat
               ! Get the chemical equilibrium
               call state%N_init(N=N_h,N_h=N_h,T_h=T_h)
               call state%equilibrate()
               ! Get the new volume of the system
               vol=0.0_WP
               do isc=1,ns
                  vol(sc%phase(isc))=vol(sc%phase(isc))+state%N(isc)*MM(isc)
               end do
               vol=vol/sc%Prho
               V=sum(vol)
               ! Relax the interface temperature
               sc%SC(i,j,k,iTl)=state%T
               sc%SC(i,j,k,iTg)=state%T
               T(i,j,k)=state%T
               ! Relax the interface VOF (Need to update the corresponding IRL quantities in vf?)
               vf%VF(i,j,k)=vol(Lphase)/V
               sc%PVF(i,j,k,Lphase)=vf%VF(i,j,k)
               sc%PVF(i,j,k,Gphase)=1.0_WP-vf%VF(i,j,k)
               ! Relax the interface composition
               do isc=1,ns
                  p=sc%phase(isc)
                  sc%SC(i,j,k,isc)=MM(isc)*state%N(isc)/(sc%Prho(p)*cfg%vol(i,j,k)*sc%PVF(i,j,k,p))
               end do
               ! Get the phase change mass flux
               evp%mdotdp(i,j,k)=(V-cfg%vol(i,j,k))/(time%dt*(1.0_WP/sc%Prho(Gphase)-1.0_WP/sc%Prho(Lphase))*cfg%vol(i,j,k)*vf%SD(i,j,k))
            end do
            ! Sync fields
            do isc=1,ns+2
               call cfg%sync(sc%SC(:,:,:,isc))
            end do
            call cfg%sync(vf%VF)
            call cfg%sync(sc%PVF(:,:,:,Lphase))
            call cfg%sync(sc%PVF(:,:,:,Gphase))
            call cfg%sync(evp%mdotdp)
            ! Deallocate arrays
            deallocate(N_h,vol)
         end block interface_jump

         ! Get the volumetric evaporation mass flux
         call evp%get_mflux()

         ! Shift the evaporation mass flux
         call evp%shift_mflux()
         
         ! Get the phase-change induced divergence
         call evp%get_div()

         ! Advance flow
         advance_flow: block
            ! Prepare new staggered viscosity (at n+1)
            call fs%get_viscosity(vf=vf,strat=harmonic_visc)
            
            ! Perform sub-iterations
            do while (time%it.le.time%itmax)
               
               ! Build mid-time velocity
               fs%U=0.5_WP*(fs%U+fs%Uold)
               fs%V=0.5_WP*(fs%V+fs%Vold)
               fs%W=0.5_WP*(fs%W+fs%Wold)
               
               ! Preliminary mass and momentum transport step at the interface
               call fs%prepare_advection_upwind(dt=time%dt)
               
               ! Explicit calculation of drho*u/dt from NS
               call fs%get_dmomdt(resU,resV,resW)
               
               ! Add momentum mass fluxes
               call fs%addsrc_gravity(resU,resV,resW)
               
               ! Assemble explicit residual
               resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
               resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
               resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
               
               ! Form implicit residuals
               call fs%solve_implicit(time%dt,resU,resV,resW)

               ! Apply these residuals
               fs%U=2.0_WP*fs%U-fs%Uold+resU!/fs%rho_U
               fs%V=2.0_WP*fs%V-fs%Vold+resV!/fs%rho_V
               fs%W=2.0_WP*fs%W-fs%Wold+resW!/fs%rho_W
               
               ! Apply other boundary conditions
               call fs%apply_bcond(time%t,time%dt)
               
               ! Solve Poisson equation
               call fs%update_laplacian()
               call fs%correct_mfr(src=evp%div_src)
               call fs%get_div(src=evp%div_src)
               ! call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
               call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
               fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
               fs%psolv%sol=0.0_WP
               call fs%psolv%solve()
               call fs%shift_p(fs%psolv%sol)
               
               ! Correct velocity
               call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
               call cfg%integrate(fs%psolv%rhs,prhs_int)
               fs%P=fs%P+fs%psolv%sol
               fs%U=fs%U-time%dt*resU/fs%rho_U
               fs%V=fs%V-time%dt*resV/fs%rho_V
               fs%W=fs%W-time%dt*resW/fs%rho_W
               
               ! Increment sub-iteration counter
               time%it=time%it+1
               
            end do
            
            ! Recompute interpolated velocity and divergence
            call fs%interp_vel(Ui,Vi,Wi)
            call fs%get_div(src=evp%div_src)

            ! Debug
            call fs%get_mfr()
            call evp%cfg%integrate(evp%div_src,mfr_err)
            mfr_err=abs(mfr_err-sum(fs%mfr))

         end block advance_flow
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Update bubble radius
         R=get_R()
         R_ext=get_Rext(time%t+t0)

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call sc%get_max()
         call mfile%write()
         call cflfile%write()
         call scfile%write()
         call evpfile%write()
         
      end do

   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects-need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,resSC,T)

   end subroutine simulation_final
   
   
end module simulation