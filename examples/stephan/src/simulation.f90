!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs,VFlo,VFhi
   use tpscalar_class,    only: tpscalar,Lphase,Gphase
   use lgpc_class,        only: lgpc
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: ss,vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(tpscalar),    public :: sc
   type(lgpc),        public :: lg
   type(timetracker), public :: time,timeSC
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,scfile,lgpcfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resSC
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: T
   
   !> Problem definition
   real(WP) :: H0
   integer  :: iTl,iTg,iYv
   real(WP) :: rho_l,rho_g,k_l,k_g,Cp_l,Cp_g,alpha_l,alpha_g,h_lg,T_sat
   real(WP) :: T_w,beta,t0
   real(WP) :: x_itf,u_itf,x_itf_ext,u_itf_ext
   real(WP) :: mdot2p,prhs_int
   
contains


   !> Function that defines a level set function
   function levelset_vapor(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the drop
      G=xyz(1)-H0
   end function levelset_vapor


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
   

   !> Function that governs the 1D Stephan beta
   function f_beta(b)
      use mathtools, only: Pi
      real(WP) :: b
      real(WP) :: f_beta
      f_beta=b*exp(b**2)*erf(b)-Cp_g*(T_w-T_sat)/(h_lg*sqrt(Pi))
   end function
   

   !> Derivative of the function that governs the 1D Stephan beta
   function fp_beta(b)
      use mathtools, only: Pi
      real(WP) :: b
      real(WP) :: fp_beta
      fp_beta=(1.0_WP+2.0_WP*b**2)*exp(b**2)*erf(b)+2.0_WP*b/sqrt(Pi)
   end function


   subroutine get_interface()
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_MAX,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      integer :: i,j,ierr
      real(WP) :: my_x_itf,my_u_itf,vol,my_vol
      ! Get the interface location
      my_x_itf=0.0_WP
      my_u_itf=0.0_WP
      do i=vf%cfg%imin_,vf%cfg%imax_
         if (vf%VF(i,1,1).gt.VFlo.and.vf%VF(i,1,1).lt.VFhi) then
            my_x_itf=vf%cfg%xm(i)
            exit
         end if
      end do
      call MPI_ALLREDUCE(my_x_itf,x_itf,1,MPI_REAL_WP,MPI_MAX,vf%cfg%comm,ierr)
      my_vol=0.0_WP
      ! Get the interface velocity
      do j=vf%cfg%jmin_,vf%cfg%jmax_
         do i=vf%cfg%imin_,vf%cfg%imax_
            if (vf%VF(i,j,1).gt.VFlo.and.vf%VF(i,j,1).lt.VFhi) then
               my_vol=my_vol+fs%cfg%vol(i,j,1)
               my_u_itf=my_u_itf+fs%cfg%vol(i,j,1)*fs%U(i+1,j,1)
               exit
            end if
         end do
      end do
      call MPI_ALLREDUCE(my_u_itf,u_itf,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_vol,vol,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      u_itf=u_itf/vol
      ! Get the analytical solution
      x_itf_ext=2.0_WP*beta*sqrt(alpha_g*(time%t))
      u_itf_ext=beta*sqrt(alpha_g/(time%t))
   end subroutine get_interface


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
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
      end block allocate_work_arrays
      
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='Main')
         call param_read('Initial time',t0)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Sub-iterations',time%itmax)
         time%t=t0
         time%dt=time%dtmax
         timeSC=timetracker(amRoot=cfg%amRoot,name='SC Time')
         call param_read('Scalar time step',timeSC%dtmax)
         call param_read('Scalar sub-iterations',timeSC%itmax)
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
         call param_read('Wall temperature',T_w)
         alpha_l=k_l/(rho_l*Cp_l)
         alpha_g=k_g/(rho_g*Cp_g)
      end block read_inputs


      ! Analytical solution
      analytical_solution: block
         use mpi_f08,  only: MPI_BCAST
         use parallel, only: MPI_REAL_WP
         real(WP) :: err,tol,beta0
         integer :: it,itmax,ierr
         if (cfg%amRoot) then
            itmax=200
            tol=1e-7
            err=10*tol
            it=0
            beta=1
            do while (err.ge.tol)
               it=it+1
               if (it.ge.itmax) exit
               beta0=beta
               beta=beta0-f_beta(beta0)/fp_beta(beta0)
               err=abs((beta-beta0)/beta0)
            end do
            print*,'1D Stephan beta = ',beta
            print*,'Relative error = ',err
            print*,'Newton-Raphson iterations = ',it
         end if
         call MPI_BCAST(beta,1,MPI_REAL_WP,0,cfg%comm,ierr)
      end block analytical_solution
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo,flux_storage,neumann
         use mpi_f08,   only: MPI_ALLREDUCE,MPI_MAX
         use parallel,  only: MPI_REAL_WP
         real(WP) :: my_x_itf
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
         ! Initialize the VOF field
         H0=2.0_WP*beta*sqrt(alpha_g*t0)
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_vapor,0.0_WP,amr_ref_lvl)
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
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
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
         ! Initialize the interface location
         do i=vf%cfg%imin_,vf%cfg%imax_
            if (vf%VF(i,1,1).gt.VFlo.and.vf%VF(i,1,1).lt.VFhi) then
               my_x_itf=vf%cfg%xm(i)
               exit
            end if
         end do
         call MPI_ALLREDUCE(my_x_itf,x_itf,1,MPI_REAL_WP,MPI_MAX,vf%cfg%comm,ierr)
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
         call fs%add_bcond(name='xp',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
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
         integer :: i,j,k,n
         ! Create scalar solver
         call sc%initialize(cfg=cfg,nscalar=2,name='tpscalar')
         ! Boundary conditinos
         call sc%add_bcond(name='xm',type=dirichlet,locator=xm_locator_sc,dir='-x')
         call sc%add_bcond(name='xp',type=neumann,  locator=xp_locator   ,dir='+x')
         sc%SCname=[  'Tl',  'Tg']; iTl=1; iTg=2
         sc%phase =[Lphase,Gphase]
         sc%diff(:,:,:,iTl)=alpha_l
         sc%diff(:,:,:,iTg)=alpha_g
         ! Initialize the linear solver
         ss=ddadi(cfg=cfg,name='Scalar',nst=7)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
         ! Initialize scalar fields
         do i=sc%cfg%imin_,sc%cfg%imax_
            do j=sc%cfg%jmin_,sc%cfg%jmax_
               do k=sc%cfg%kmin_,sc%cfg%kmax_
                  if (vf%VF(i,j,k).gt.VFlo) sc%SC(i,j,k,iTl)=T_sat
                  ! if (vf%VF(i,j,k).lt.VFhi) sc%SC(i,j,k,iTg)=T_w+(T_sat-T_w)/erf(beta)*erf(sc%cfg%xm(i)/(2.0_WP*sqrt(k_g/(fs%rho_g*Cp_g)*t0)))
                  if (vf%VF(i,j,k).lt.VFhi) sc%SC(i,j,k,iTg)=T_w+(T_sat-T_w)/x_itf*sc%cfg%xm(i)
               end do
            end do
         end do
         where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
            sc%SC(:,:,:,iTl)=T_sat
            sc%SC(:,:,:,iTg)=T_sat
         end where
         call sc%cfg%sync(sc%SC(:,:,:,iTl))
         call sc%cfg%sync(sc%SC(:,:,:,iTg))
         ! Apply boundary conditions
         call sc%apply_bcond(time%t,time%dt)
         call sc%get_bcond('xm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            sc%SC(i,j,k,iTg)=2.0_WP*T_w-sc%SC(i+1,j,k,iTg)
         end do
         ! Initialize the phase specific density and VOF
         sc%Prho(Lphase)=fs%rho_l
         sc%Prho(Gphase)=fs%rho_g
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         ! Post process
         T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)
         ! Get the phasic face apertures
         call sc%get_face_apt()
      end block create_scalar
      

      ! Create and initialize an lgpc object
      create_lgpc: block
         integer :: index,index_pure,i,j,k
         ! Create the object
         call lg%initialize(cfg=cfg,vf=vf,sc=sc%SC,iTl=iTl,iTg=iTg,itp_x=fs%itpr_x,itp_y=fs%itpr_y,itp_z=fs%itpr_z,div_x=fs%divp_x,div_y=fs%divp_y,div_z=fs%divp_z,name='liquid gas pc')
         call param_read('Mass flux tolerence',     lg%mdot3p_tol)
         call param_read('Max pseudo timestep size',lg%pseudo_time%dtmax)
         call param_read('Max pseudo cfl number',   lg%pseudo_time%cflmax)
         call param_read('Max pseudo time steps',   lg%pseudo_time%nmax)
         lg%pseudo_time%dt=lg%pseudo_time%dtmax
         ! Get densities from the flow solver
         lg%rho_l=fs%rho_l
         lg%rho_g=fs%rho_g
         ! Get temperature gradient
         call lg%get_temperature_grad()
         ! Phase change mass flux
         where ((vf%VF.gt.VFlo).and.(vf%VF.lt.VFhi))
            lg%mdot2p=(k_g*lg%Tg_grd-k_l*lg%Tl_grd)/h_lg
         else where
            lg%mdot2p=0.0_WP
         end where
         ! Get the volumetric pahse change mass flux
         call lg%get_mdot3p()
         ! Initialize the liquid and gas side mass fluxes
         call lg%init_mdot3pLG()
      end block create_lgpc


      ! Initialize the interface location and velocity
      call get_interface()


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         integer :: nsc
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='stephan')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_surface('plic',smesh)
         do nsc=1,sc%nscalar
           call ens_out%add_scalar(trim(sc%SCname(nsc)),sc%SC(:,:,:,nsc))
         end do
         call ens_out%add_scalar('mdot2p',lg%mdot2p)
         call ens_out%add_scalar('mdot3p',lg%mdot3p)
         call ens_out%add_scalar('mdot3pL',lg%mdot3pLG(:,:,:,Lphase))
         call ens_out%add_scalar('mdot3pG',lg%mdot3pLG(:,:,:,Gphase))
         call ens_out%add_scalar('lgpc_div',lg%div_vel)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('Temperature',T)
         call ens_out%add_vector('normal',lg%normal(:,:,:,1),lg%normal(:,:,:,2),lg%normal(:,:,:,3))
         call ens_out%add_scalar('Tl_grd',lg%Tl_grd)
         call ens_out%add_scalar('Tg_grd',lg%Tg_grd)
         ! Output to ensight
         call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         integer :: nsc
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
         call mfile%add_column(x_itf,'x_itf')
         call mfile%add_column(u_itf,'u_itf')
         call mfile%add_column(x_itf_ext,'x_itf_ext')
         call mfile%add_column(u_itf_ext,'u_itf_ext')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%add_column(prhs_int,'prhs_int')
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
         call scfile%add_column(timeSC%n,'Timestep number')
         call scfile%add_column(timeSC%t,'Time')
         do nsc=1,sc%nscalar
           call scfile%add_column(sc%SCmin(nsc),trim(sc%SCname(nsc))//'_min')
           call scfile%add_column(sc%SCmax(nsc),trim(sc%SCname(nsc))//'_max')
           call scfile%add_column(sc%SCint(nsc),trim(sc%SCname(nsc))//'_int')
         end do
         call scfile%write()
         ! Create lgpc monitor
         lgpcfile=monitor(lg%cfg%amRoot,'lgpc')
         call lgpcfile%add_column(time%n,'Timestep number')
         call lgpcfile%add_column(time%t,'Time')
         call lgpcfile%add_column(lg%pseudo_time%dt,'Pseudo time step')
         call lgpcfile%add_column(lg%pseudo_time%cfl,'Maximum pseudo CFL')
         call lgpcfile%add_column(lg%pseudo_time%n,'No. pseudo steps')
         call lgpcfile%add_column(lg%mdot3p_int,'mdot3p int')
         call lgpcfile%add_column(lg%mdot3pL_int,'shifted mdot3pL int')
         call lgpcfile%add_column(lg%mdot3pG_int,'shifted mdot3pG int')
         call lgpcfile%add_column(lg%mdot3pL_int_err,'mdot3pL int err')
         call lgpcfile%add_column(lg%mdot3pG_int_err,'mdot3pG int err')
         call lgpcfile%add_column(lg%mdot3pL_err,'max mdot3pL err')
         call lgpcfile%add_column(lg%mdot3pG_err,'max mdot3pG err')
         call lgpcfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use tpns_class, only: static_contact,harmonic_visc
      use mathtools, only: Pi
      use tpscalar_class,  only: bcond
      implicit none
      type(bcond), pointer :: mybc
      integer  :: i,j,k,nsc,n,p
      real(WP) :: dt_sc

      ! Perform flow time integration
      do while (.not.time%done())

         ! Increment flow time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old VOF
         vf%VFold=vf%VF

         ! Remember old pahse change velocity divergence
         lg%div_vel_old=lg%div_vel
         
         ! Remember old SC
         sc%SCold =sc%SC
         sc%PVFold=sc%PVF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         call vf%apply_bcond(time%t,time%dt)

         ! Update the phasic VOF and face apertures
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         call sc%get_face_apt()

         ! ================== SCALAR ================== !
   
         ! Perform scalar time integration
         do while (timeSC%t.lt.time%t)
            
            ! Increment scalar time step
            if (timeSC%t+timeSC%dt.gt.time%t) then
               dt_sc=timeSC%dt
               timeSC%dt=time%t-timeSC%t
               call timeSC%increment()
               timeSC%dt=dt_sc
            else
               call timeSC%increment()
            end if
            
            ! Remember old scalar
            sc%SCold =sc%SC

            ! Perform sub-iterations
            do while (timeSC%it.le.timeSC%itmax)
               if (cfg%amRoot) print*,'Scalar sub-iteration ',timeSC%it

               ! Build mid-time scalar
               do nsc=1,sc%nscalar
                  p=sc%phase(nsc)
                  do k=cfg%kmino_,cfg%kmaxo_
                     do j=cfg%jmino_,cfg%jmaxo_
                        do i=cfg%imino_,cfg%imaxo_
                           if (sc%PVF(i,j,k,p).eq.1.0_WP) then
                              sc%SC(i,j,k,nsc)=0.5_WP*(sc%SC(i,j,k,nsc)+sc%SCold(i,j,k,nsc))
                           else if (sc%PVF(i,j,k,p).gt.0.0_WP) then
                              sc%SC(i,j,k,nsc)=T_sat
                           else
                              sc%SC(i,j,k,nsc)=0.0_WP
                           end if
                        end do
                     end do
                  end do
               end do

               ! Explicit calculation of dSC/dt
               call sc%get_dSCdt(dSCdt=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old)

               ! Assemble explicit residual
               do nsc=1,sc%nscalar
                  p=sc%phase(nsc)
                  do k=cfg%kmin_,cfg%kmax_
                     do j=cfg%jmin_,cfg%jmax_
                        do i=cfg%imin_,cfg%imax_
                           if (sc%PVF(i,j,k,p).eq.1.0_WP) then
                              resSC(i,j,k,nsc)=2.0_WP*(sc%SCold(i,j,k,nsc)-sc%SC(i,j,k,nsc))+timeSC%dt*resSC(i,j,k,nsc)
                           else
                              resSC(i,j,k,nsc)=0.0_WP
                           end if
                        end do
                     end do
                  end do
               end do

               ! Form implicit residual
               call sc%solve_implicit(dt=timeSC%dt,resSC=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old)

               ! Apply the residuals
               sc%SC=2.0_WP*sc%SC-sc%SCold+resSC

               ! Zero-out the scalar in the opposite phase
               do nsc=1,sc%nscalar
                  where (sc%PVF(:,:,:,sc%phase(nsc)).lt.VFlo) sc%SC(:,:,:,nsc)=0.0_WP
               end do

               ! One-field temperature
               T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)

               ! Apply boundary conditions
               where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
                  sc%SC(:,:,:,iTl)=T_sat
                  sc%SC(:,:,:,iTg)=T_sat
               end where
               call sc%apply_bcond(timeSC%t,timeSC%dt)
               call sc%get_bcond('xm',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  sc%SC(i,j,k,iTg)=2.0_WP*T_w-sc%SC(i+1,j,k,iTg)
               end do
         
               ! Increment scalar sub-iteration counter
               timeSC%it=timeSC%it+1

            end do

         end do


         ! ================== PHASE CHANGE ================== !

         ! Get temperature gradient
         call lg%get_temperature_grad()

         ! Phase change mass flux
         where ((vf%VF.gt.VFlo).and.(vf%VF.lt.VFhi))
            lg%mdot2p=(k_g*lg%Tg_grd-k_l*lg%Tl_grd)/h_lg
         else where
            lg%mdot2p=0.0_WP
         end where

         ! Get the volumetric phase change mass flux
         call lg%get_mdot3p()

         ! Shift the phase change mass flux
         call lg%shift_mdot3p()
         
         ! Get the phase change induced divergence
         call lg%get_div()


         ! ================== VELOCITY ================== !

         ! Perform velocity sub-iterations
         do while (time%it.le.time%itmax)
            if (cfg%amRoot) print*,'Flow sub iteration ',time%it
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Prepare new staggered viscosity (at n+1)
            call fs%get_viscosity(vf=vf,strat=harmonic_visc)

            ! Prepare old staggered density (at n)
            call fs%get_olddensity(vf=vf)

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
            call fs%correct_mfr(src=lg%div_vel)
            call fs%get_div(src=lg%div_vel)
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
            
            ! Increment velocity sub-iteration counter
            time%it=time%it+1
            
         end do

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(src=lg%div_vel)


         ! ================== OUTPUT ================== !

         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Update the interface location and velocity
         call get_interface()

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call sc%get_max()
         call mfile%write()
         call cflfile%write()
         call scfile%write()
         call lgpcfile%write()
         
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
      deallocate(resU,resV,resW,Ui,Vi,Wi,resSC,T)

   end subroutine simulation_final
   
   
end module simulation