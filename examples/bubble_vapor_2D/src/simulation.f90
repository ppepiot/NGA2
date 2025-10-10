!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,Lz
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
   use mathtools,         only: Pi
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ss,ps
   type(ddadi),       public :: vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(tpscalar),    public :: sc
   type(lgpc),        public :: lg
   type(timetracker), public :: time,timeSC
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt,T_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,scfile,lgfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resSC
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: T
   
   !> Problem definition
   real(WP) :: R0,R,R_ext,V_b
   real(WP), dimension(3) :: center
   integer  :: iTl,iTg
   real(WP) :: rho_l,rho_g,k_l,k_g,Cp_l,Cp_g,alpha_l,alpha_g,h_lg,T_sat
   real(WP) :: T_inf,beta,f_b_cnst,t0
   ! Debug
   real(WP) :: prhs_int
   real(WP) :: mfr_err
   real(WP) :: Tlgrd_min,Tlgrd_max,Tlgrd_avg,Tlgrd_ext
   type(monitor) :: tlgrd_file
   
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


   subroutine apply_dirichlet()
      use tpns_class, only: bcond,dirichlet
      use mathtools,  only: Pi
      type(bcond), pointer :: my_bc
      real(WP) :: Ub,Ux,Uy,vfr,myR
      integer  :: i,j,k,n,stag
      call cfg%integrate(lg%div_vel,vfr)
      my_bc=>fs%first_bc
      do while (associated(my_bc))
         if (my_bc%type.ne.dirichlet) cycle
         if (my_bc%itr%amIn) then
            select case (my_bc%face)
            case ('x')
               stag=min(my_bc%dir,0)
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  myR=sqrt(cfg%xm(i)**2+cfg%ym(j)**2+cfg%zm(k)**2)
                  Ub=vfr/(2.0_WP*myR*Pi*Lz)
                  Ux=cfg%xm(i)/myR*Ub
                  Uy=cfg%ym(j)/myR*Ub
                  fs%U(i     ,j    ,k    )=Ux
                  fs%V(i+stag,j:j+1,k    )=Uy
                  fs%W(i+stag,j    ,k:k+1)=0.0_WP
               end do
            case ('y')
               stag=min(my_bc%dir,0)
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  myR=sqrt(cfg%xm(i)**2+cfg%ym(j)**2+cfg%zm(k)**2)
                  Ub=vfr/(2.0_WP*myR*Pi*Lz)
                  Ux=cfg%xm(i)/myR*Ub
                  Uy=cfg%ym(j)/myR*Ub
                  fs%U(i:i+1,j+stag,k    )=Ux
                  fs%V(i    ,j     ,k    )=Uy
                  fs%W(i    ,j+stag,k:k+1)=0.0_WP
               end do
            ! case ('z')
            !    stag=min(my_bc%dir,0)
            !    do n=1,my_bc%itr%n_
            !       i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            !       myR=sqrt(cfg%xm(i)**2+cfg%ym(j)**2+cfg%zm(k)**2)
            !    end do
            end select
         end if
         my_bc=>my_bc%next
      end do
   end subroutine apply_dirichlet
   

   ! Integrand function for beta equation
   function beta_int(z,b)
      real(WP), intent(in) :: z,b
      real(WP) :: beta_int
      real(WP), parameter :: eps = 1.0e-12_WP
      real(WP) :: denom
      denom=max(1.0_WP-z,eps)
      beta_int=exp(-b**2*(denom**(-2)-2.0_WP*(1.0_WP-rho_g/rho_l)*z-1.0_WP))
   end function beta_int


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


   ! Simpson's rule integration
   function Simpson_adpt(integrand,b,zmin,zmax,tol) result(S2)
      interface
      function integrand(x1,x2)
            use precision, only: WP
            real(WP), intent(in) :: x1,x2
            real(WP) :: integrand
         end function integrand
      end interface
      real(WP), intent(in) :: b
      real(WP), intent(in) :: zmin,zmax
      real(WP), intent(in) :: tol
      integer  :: i,n
      real(WP) :: S1,S2,err
      n=2
      S2=Simpson(integrand,b,zmin,zmax,n)
      err=10.0_WP*tol
      do while (err.gt.tol)
         n=2*n
         S1=S2
         S2=Simpson(integrand,b,zmin,zmax,n)
         err=abs((S2-S1)/S1)
      end do
   end function Simpson_adpt


   !> Function that governs beta
   function f_beta(b)
      real(WP), intent(in) :: b
      real(WP) :: f_beta
      f_beta=2.0_WP*b**2*Simpson(beta_int,b,0.0_WP,1.0_WP,20)-f_b_cnst
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
      get_R=sqrt(V_b/(Pi*Lz))
   end function get_R


   ! Debug
   subroutine get_tlgrd()
      use mpi_f08,   only: MPI_MAX,MPI_MIN,MPI_SUM,MPI_ALLREDUCE
      use parallel,  only: MPI_REAL_WP
      integer :: index,i,j,k,ierr
      real(WP) :: my_area,area,my_Tlgrd_min,my_Tlgrd_max,my_Tlgrd_avg,tlgrd
      Tlgrd_ext=2.0_WP*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)/get_Rext(time%t)
      my_Tlgrd_min=1e6
      my_Tlgrd_max=0.0_WP
      my_Tlgrd_avg=0.0_WP
      my_area=0.0_WP
      do index=1,vf%band_count(0)
         i=vf%band_map(1,index)
         j=vf%band_map(2,index)
         k=vf%band_map(3,index)
         tlgrd=lg%Tl_grd(i,j,k)
         if (tlgrd.lt.my_Tlgrd_min) my_Tlgrd_min=tlgrd
         if (tlgrd.gt.my_Tlgrd_max) my_Tlgrd_max=tlgrd
         my_area=my_area+cfg%vol(i,j,k)*vf%SD(i,j,k)
         my_Tlgrd_avg=my_Tlgrd_avg+cfg%vol(i,j,k)*vf%SD(i,j,k)*abs(tlgrd)
      end do
      call MPI_ALLREDUCE(my_Tlgrd_min,Tlgrd_min,1,MPI_REAL_WP,MPI_MIN,cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Tlgrd_max,Tlgrd_max,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
      call MPI_ALLREDUCE(my_area,area,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Tlgrd_avg,Tlgrd_avg,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      Tlgrd_avg=Tlgrd_avg/area
   end subroutine get_tlgrd


   !> Calculate the asimuthally averaged liquid temperature and print it along with the analytical values
   subroutine get_Tl()
      use messager, only: die
      use mathtools, only: Pi
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_INT
      use parallel,  only: MPI_REAL_WP
      type(monitor) :: Tlfile
      real(WP) :: theta_s,theta_e,dtheta,theta
      real(WP) :: nx,ny
      real(WP) :: l,dR,R_e,x,y,z,vof,rad_num,rad_ext,Tl_num,Tl_ext
      real(WP), dimension(:), allocatable :: r_num,r_extt,T_l,T_l_avg
      integer  :: ntheta,itheta,nT_num,nT_ext,ind(3),i,j,k,T_ind,ierr
      character(len=str_medium) :: tstr
      integer :: found_local,found_global
      ! Create the monitors
      write(tstr,'(F6.3)') time%t
      Tlfile=monitor(fs%cfg%amRoot,'Tl_'//trim(adjustl(tstr)))
      call Tlfile%add_column(rad_num,'r_num')
      call Tlfile%add_column(Tl_num,'Tl_num')
      call Tlfile%add_column(rad_ext,'r_ext')
      call Tlfile%add_column(Tl_ext,'Tl_ext')
      ! Initialize
      theta_s=0.0_WP
      theta_e=2.0_WP*Pi
      ntheta=20
      dtheta=(theta_e-theta_s)/real(ntheta-1,WP)
      R_e=0.45_WP
      dR=cfg%dx(1)
      ! Create the radial grid
      nT_num=int((R_e-R)/dR)+1
      nT_ext=int((R_e-R_ext)/dR)+1
      allocate(r_num(nT_num))
      allocate(r_extt(nT_ext))
      allocate(T_l(nT_num)); T_l=0.0_WP
      allocate(T_l_avg(nT_num)); T_l_avg=0.0_WP
      r_num(1)=R
      do T_ind=2,nT_num
         r_num(T_ind)=r_num(T_ind-1)+dR
      end do
      r_extt(1)=R_ext
      do T_ind=2,nT_ext
         r_extt(T_ind)=r_extt(T_ind-1)+dR
      end do
      ! Loop over lines
      do itheta=1,ntheta
         ! Reset the flags
         found_local =0
         found_global=0
         ! Advance theta
         theta=theta_s+real(itheta-1,WP)*dtheta
         ! Get the normal vector
         nx=cos(theta)
         ny=sin(theta)
         ! Initilize the line marcher
         l=0.0_WP
         vof=0.0_WP
         i=(cfg%imax_+cfg%imin_)/2
         j=(cfg%jmax_+cfg%jmin_)/2
         k=(cfg%kmax_+cfg%kmin_)/2
         ! March the line until hit interfacial cell
         do while (found_global.eq.0)
            ! Advance the line segment
            l=l+dR
            x=l*nx
            y=l*ny
            z=cfg%zm(1)
            if (cfg%is_in_subdomain([x,y,z])) then
               ! Get the corresponding indices of the cell containing the line segment point
               ind=cfg%get_ijk_local(pos=[x,y,z],ind_guess=[i,j,k])
               i=ind(1)
               j=ind(2)
               k=ind(3)
               ! Get the VOF
               vof=vf%VF(i,j,k)
            else
               vof=0.0_WP
            end if
            if ((vof.gt.0.0_WP).and.(vof.lt.1.0_WP)) found_local=1
            ! Check if any process has found it
            call MPI_ALLREDUCE(found_local,found_global,1,MPI_INT,MPI_SUM,cfg%comm,ierr)
         end do
         if (found_global.eq.0) then
            if (cfg%amRoot) print*,'Warning: Line marcher could not find any interfacial cell'
            cycle
         end if
         ! Get the first temperature exactly at the interface
         if (found_local.eq.1) T_l(1)=T_l(1)+sc%SC(i,j,k,iTl)
         ! March the line from the interface to the end of the line
         do T_ind=2,nT_num
            ! Get the point coordinates
            x=r_num(T_ind)*nx
            y=r_num(T_ind)*ny
            z=cfg%zm(1)
            if (cfg%is_in_subdomain([x,y,z])) then
               ! Get the corresponding indices of the cell containing the line segment point
               ind=cfg%get_ijk_local(pos=[x,y,z],ind_guess=[i,j,k])
               i=ind(1)
               j=ind(2)
               k=ind(3)
               ! Get temperature
               T_l(T_ind)=T_l(T_ind)+cfg%get_scalar([x,y,z],i,j,k,sc%SC(:,:,:,iTl),'N')
            end if
         end do
      end do
      ! Calculate the mean temperature
      call MPI_ALLREDUCE(T_l,T_l_avg,nT_num,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      T_l_avg=T_l_avg/real(ntheta,WP)
      do T_ind=1,min(nT_num,nT_ext)
         rad_num=r_num(T_ind)
         Tl_num=T_l_avg(T_ind)
         rad_ext=r_extt(T_ind)
         Tl_ext=T_inf-2.0_WP*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)*Simpson(beta_int,beta,1.0_WP-R_ext/rad_ext,1.0_WP,20)
         call Tlfile%write()
      end do
      ! Deallocate the arrays
      deallocate(r_num,T_l,T_l_avg)
   end subroutine get_Tl


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
         time=timetracker(amRoot=cfg%amRoot,name='Main time')
         call param_read('Initial time',t0)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Sub-iterations',time%itmax)
         time%t=t0
         time%dt=time%dtmax
         timeSC=timetracker(amRoot=cfg%amRoot,name='Scalar time')
         call param_read('Scalar time step',timeSC%dtmax)
         call param_read('Scalar sub-iterations',timeSC%itmax)
         timeSC%t=time%t
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
         call param_read('Far-field temperature',T_inf);
         alpha_l=k_l/(rho_l*Cp_l)
         alpha_g=k_g/(rho_g*Cp_g)
         f_b_cnst=(rho_l*Cp_l*(T_inf-T_sat))/(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))
      end block read_inputs


      ! Analytical solution
      analytical_solution: block
         use messager, only: die
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
            betaR=1.0_WP
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
            it=min(it,itmax)
            print*,'Convergence = ',convergence
            print*,'beta =',beta
            print*,'Absolute error =',err
            print*,'Bi-section iterations =',it
            if (.not.convergence) call die('[simulation_init] Bi-section failed')
            ! Debug
            print*,'Analytical T grad = ',2.0_WP*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)/get_Rext(t0)
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
         use tpns_class,      only: bcond,dirichlet
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
         call fs%add_bcond(name='xm',type=dirichlet,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
         call fs%add_bcond(name='xp',type=dirichlet,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         call fs%add_bcond(name='ym',type=dirichlet,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call fs%add_bcond(name='yp',type=dirichlet,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
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
         use param,           only: param_read
         use tpscalar_class,  only: bcond,dirichlet
         use hypre_str_class, only: pcg_smg
         type(bcond), pointer :: my_bc
         real(WP) :: radius
         integer  :: i,j,k,n
         ! Create scalar solver
         call sc%initialize(cfg=cfg,nscalar=2,name='tpscalar')
         ! Boundary conditinos
         call sc%add_bcond(name='xm',type=dirichlet,locator=xm_locator_sc,dir='-x')
         call sc%add_bcond(name='xp',type=dirichlet,locator=xp_locator   ,dir='+x')
         call sc%add_bcond(name='ym',type=dirichlet,locator=ym_locator_sc,dir='-y')
         call sc%add_bcond(name='yp',type=dirichlet,locator=yp_locator   ,dir='+y')
         sc%SCname=[  'Tl',  'Tg']; iTl=1; iTg=2
         sc%phase =[Lphase,Gphase]
         sc%diff(:,:,:,iTl)=alpha_l
         sc%diff(:,:,:,iTg)=alpha_g
         ! Initialize the linear solver
         ! ss=ddadi(cfg=cfg,name='Scalar',nst=7)
         ss=hypre_str(cfg=cfg,name='Scalar',method=pcg_smg,nst=7)
         call param_read('Scalar iteration',ss%maxit)
         call param_read('Scalar tolerance',ss%rcvg)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
         ! Initialize scalar fields
         do i=sc%cfg%imino_,sc%cfg%imaxo_
            do j=sc%cfg%jmino_,sc%cfg%jmaxo_
               do k=sc%cfg%kmino_,sc%cfg%kmaxo_
                  radius=norm2([sc%cfg%xm(i),sc%cfg%ym(j),sc%cfg%zm(k)])
                  if (vf%VF(i,j,k).gt.VFlo) sc%SC(i,j,k,iTl)=T_inf-2.0_WP*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)*Simpson(beta_int,beta,1.0_WP-R0/radius,1.0_WP,20)
                  if (vf%VF(i,j,k).lt.VFhi) sc%SC(i,j,k,iTg)=T_sat
               end do
            end do
         end do
         radius=R0+0.001_WP
         ! print*,'T = ',T_inf-2.0_WP*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)*Simpson(beta_int,beta,1.0_WP-R0/radius,1.0_WP,20)
         ! print*,'T grad = ',((T_inf-2.0_WP*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)*Simpson(beta_int,beta,1.0_WP-R0/radius,1.0_WP,20))-T_sat)/0.001_WP
         ! print*,'dx = ',cfg%dx(1)
         ! Apply boundary conditions
         where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
            sc%SC(:,:,:,iTl)=T_sat
            sc%SC(:,:,:,iTg)=T_sat
         end where
         call sc%apply_bcond(time%t,time%dt)
         call sc%get_bcond('xm',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i+1,j,k,iTl)
         end do
         call sc%get_bcond('xp',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i-1,j,k,iTl)
         end do
         call sc%get_bcond('ym',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j+1,k,iTl)
         end do
         call sc%get_bcond('yp',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j-1,k,iTl)
         end do
         ! Initialize the phasic density and VOF
         sc%Prho(Lphase)=fs%rho_l
         sc%Prho(Gphase)=fs%rho_g
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         ! Get the phasic face apertures
         call sc%get_face_apt()
         ! Post process
         T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)
      end block create_scalar
      

      ! Create and initialize an lg object
      create_lgpc: block
         ! use lgpc_class, only: symmetry
         integer :: i,j,k,index!,ind(3)
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
         ! 
         ! ind=cfg%get_ijk_global([0.113_WP,0.001_WP,0.0_WP],[20,32,1])
         ! print*,'ind = ',ind
         ! print*,'xm = ',cfg%xm(ind(1))
         ! print*,'ym = ',cfg%ym(ind(2))
         ! print*,'zm = ',cfg%zm(ind(3))
         ! Get temperature gradient
         call lg%get_temperature_grad()
         ! Phase change mass flux
         lg%mdot2p=0.0_WP
         do k=lg%cfg%kmin_,lg%cfg%kmax_
            do j=lg%cfg%jmin_,lg%cfg%jmax_
               do i=lg%cfg%imin_,lg%cfg%imax_
                  if ((vf%VF(i,j,k).gt.VFlo).and.(vf%VF(i,j,k).lt.VFhi)) then
                     lg%mdot2p(i,j,k)=(-k_g*lg%Tg_grd(i,j,k)+k_l*lg%Tl_grd(i,j,k))/h_lg
                  end if
               end do
            end do
         end do
         call lg%cfg%sync(lg%mdot2p)
         ! Get the volumetric evaporation mass flux
         call lg%get_mdot3p()
         ! Initialize the liquid and gas side mass fluxes
         call lg%init_mdot3pLG()
         ! Get the interface normal
         call lg%get_normal()
      end block create_lgpc


      ! Get the bubble radius
      R=get_R()
      R_ext=get_Rext(time%t)


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         integer :: nsc
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='bubble_vapor_2D')
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
         call ens_out%add_scalar('mdot3p',lg%mdot3p)
         call ens_out%add_scalar('lgpc_div',lg%div_vel)
         call ens_out%add_scalar('mdot2p',lg%mdot2p)
         call ens_out%add_scalar('mdot3pL',lg%mdot3pLG(:,:,:,Lphase))
         call ens_out%add_scalar('mdot3pG',lg%mdot3pLG(:,:,:,Gphase))
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('Temperature',T)
         call ens_out%add_vector('normal',lg%normal(:,:,:,1),lg%normal(:,:,:,2),lg%normal(:,:,:,3))
         call ens_out%add_scalar('Tl_grd',lg%Tl_grd)
         call ens_out%add_scalar('Tg_grd',lg%Tg_grd)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         ! Create an event for temperature output
         T_evt=event(time=time,name='Temperature output')
         call param_read('Temperature file output period',T_evt%tper)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         integer :: nsc,i
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
         call lg%cfg%integrate(lg%div_vel,mfr_err)
         mfr_err=abs(mfr_err-sum(fs%mfr))
         call mfile%add_column(mfr_err,'mfr_err')
         call mfile%write()
         call get_tlgrd()
         tlgrd_file=monitor(cfg%amRoot,'Tl_grad')
         call tlgrd_file%add_column(time%n,'Timestep number')
         call tlgrd_file%add_column(time%t,'Time')
         call tlgrd_file%add_column(Tlgrd_min,'Tl grad min')
         call tlgrd_file%add_column(Tlgrd_max,'Tl grad max')
         call tlgrd_file%add_column(Tlgrd_avg,'Tl grad avg')
         call tlgrd_file%add_column(Tlgrd_ext,'Tl grad ext')
         call tlgrd_file%write()
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
         ! Create evaporation monitor
         lgfile=monitor(lg%cfg%amRoot,'lgpc')
         call lgfile%add_column(time%n,'Timestep number')
         call lgfile%add_column(time%t,'Time')
         call lgfile%add_column(lg%pseudo_time%dt,'Pseudo time step')
         call lgfile%add_column(lg%pseudo_time%cfl,'Maximum pseudo CFL')
         call lgfile%add_column(lg%pseudo_time%n,'No. pseudo steps')
         call lgfile%add_column(lg%mdot3p_int,'mdot3p int')
         call lgfile%add_column(lg%mdot3pL_int,'shifted mdot3pL int')
         call lgfile%add_column(lg%mdot3pG_int,'shifted mdot3pG int')
         call lgfile%add_column(lg%mdot3pL_int_err,'mdot3pL int err')
         call lgfile%add_column(lg%mdot3pG_int_err,'mdot3pG int err')
         call lgfile%add_column(lg%mdot3pL_err,'max mdot3pL err')
         call lgfile%add_column(lg%mdot3pG_err,'max mdot3pG err')
         call lgfile%write()
      end block create_monitor

      
      ! Output temperature profile
      call get_Tl()
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation-this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none

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
         
         ! Remember old scalar
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

         ! CN
         ! advance_scalar: block
         !    use tpscalar_class,  only: bcond
         !    type(bcond), pointer :: my_bc
         !    integer  :: i,j,k,nsc,n,p
         !    real(WP) :: dt_sc

         !    ! Perform scalar time integration
         !    do while (timeSC%t.lt.time%t)
               
         !       ! Increment scalar time step
         !       if (timeSC%t+timeSC%dt.gt.time%t) then
         !          dt_sc=timeSC%dt
         !          timeSC%dt=time%t-timeSC%t
         !          call timeSC%increment()
         !          timeSC%dt=dt_sc
         !       else
         !          call timeSC%increment()
         !       end if
               
         !       ! Remember old scalar
         !       sc%SCold =sc%SC

         !       ! Perform sub-iterations
         !       do while (timeSC%it.le.timeSC%itmax)
         !          if (cfg%amRoot) print*,'Scalar sub-iteration ',timeSC%it

         !          ! Build mid-time scalar
         !          do nsc=1,sc%nscalar
         !             p=sc%phase(nsc)
         !             do k=cfg%kmino_,cfg%kmaxo_
         !                do j=cfg%jmino_,cfg%jmaxo_
         !                   do i=cfg%imino_,cfg%imaxo_
         !                      if (sc%PVF(i,j,k,p).eq.1.0_WP) then
         !                         sc%SC(i,j,k,nsc)=0.5_WP*(sc%SC(i,j,k,nsc)+sc%SCold(i,j,k,nsc))
         !                      else if (sc%PVF(i,j,k,p).gt.0.0_WP) then
         !                         sc%SC(i,j,k,nsc)=T_sat
         !                      else
         !                         sc%SC(i,j,k,nsc)=0.0_WP
         !                      end if
         !                   end do
         !                end do
         !             end do
         !          end do

         !          ! Explicit calculation of dSC/dt
         !          call sc%get_dSCdt(dSCdt=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old)

         !          ! Assemble explicit residual
         !          do nsc=1,sc%nscalar
         !             p=sc%phase(nsc)
         !             do k=cfg%kmin_,cfg%kmax_
         !                do j=cfg%jmin_,cfg%jmax_
         !                   do i=cfg%imin_,cfg%imax_
         !                      if (sc%PVF(i,j,k,p).eq.1.0_WP) then
         !                         resSC(i,j,k,nsc)=2.0_WP*(sc%SCold(i,j,k,nsc)-sc%SC(i,j,k,nsc))+timeSC%dt*resSC(i,j,k,nsc)
         !                      else
         !                         resSC(i,j,k,nsc)=0.0_WP
         !                      end if
         !                   end do
         !                end do
         !             end do
         !          end do

         !          ! Form implicit residual
         !          call sc%solve_implicit(dt=timeSC%dt,resSC=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old)

         !          ! Apply the residuals
         !          sc%SC=2.0_WP*sc%SC-sc%SCold+resSC

         !          ! Zero-out the scalar in the opposite phase
         !          do nsc=1,sc%nscalar
         !             where (sc%PVF(:,:,:,sc%phase(nsc)).lt.VFlo) sc%SC(:,:,:,nsc)=0.0_WP
         !          end do

         !          ! Apply boundary conditions
         !          where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
         !             sc%SC(:,:,:,iTl)=T_sat
         !             sc%SC(:,:,:,iTg)=T_sat
         !          end where
         !          call sc%apply_bcond(timeSC%t,timeSC%dt)
         !          call sc%get_bcond('xm',my_bc)
         !          do n=1,my_bc%itr%no_
         !             i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !             sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i+1,j,k,iTl)
         !          end do
         !          call sc%get_bcond('xp',my_bc)
         !          do n=1,my_bc%itr%no_
         !             i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !             sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i-1,j,k,iTl)
         !          end do
         !          call sc%get_bcond('ym',my_bc)
         !          do n=1,my_bc%itr%no_
         !             i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !             sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j+1,k,iTl)
         !          end do
         !          call sc%get_bcond('yp',my_bc)
         !          do n=1,my_bc%itr%no_
         !             i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !             sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j-1,k,iTl)
         !          end do
            
         !          ! Increment scalar sub-iteration counter
         !          timeSC%it=timeSC%it+1

         !       end do

         !    end do

         !    ! One-field temperature
         !    T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)

         ! end block advance_scalar


         ! Backward Euler for diffusion and CN for advection
         advance_scalar: block
            use tpscalar_class,  only: bcond
            type(bcond), pointer :: my_bc
            integer  :: i,j,k,nsc,n,p
            real(WP) :: dt_sc

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

               ! Build mid-time scalar
               do nsc=1,sc%nscalar
                  p=sc%phase(nsc)
                  do k=cfg%kmino_,cfg%kmaxo_
                     do j=cfg%jmino_,cfg%jmaxo_
                        do i=cfg%imino_,cfg%imaxo_
                           if (sc%PVF(i,j,k,p).eq.0.0_WP) then
                              sc%SC(i,j,k,nsc)=0.0_WP
                           else if (sc%PVF(i,j,k,p).lt.1.0_WP) then
                              sc%SC(i,j,k,nsc)=T_sat
                           else
                           end if
                        end do
                     end do
                  end do
               end do

               ! Explicit calculation of dSC/dt
               call sc%get_dSCdt(dSCdt=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old)

               ! Assemble explicit residual
               resSC=resSC*timeSC%dt
               do nsc=1,sc%nscalar
                  p=sc%phase(nsc)
                  do k=cfg%kmin_,cfg%kmax_
                     do j=cfg%jmin_,cfg%jmax_
                        do i=cfg%imin_,cfg%imax_
                           if (sc%PVF(i,j,k,p).lt.1.0_WP) then
                              resSC(i,j,k,nsc)=0.0_WP
                           end if
                        end do
                     end do
                  end do
               end do

               ! Form implicit residual
               call sc%solve_implicit(dt=timeSC%dt,resSC=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old,w_adv=0.5_WP,w_dff=1.0_WP)

               ! Apply the residuals
               sc%SC=sc%SCold+resSC

               ! Apply boundary conditions
               where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
                  sc%SC(:,:,:,iTl)=T_sat
                  sc%SC(:,:,:,iTg)=T_sat
               end where
               call sc%apply_bcond(timeSC%t,timeSC%dt)
               call sc%get_bcond('xm',my_bc)
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i+1,j,k,iTl)
               end do
               call sc%get_bcond('xp',my_bc)
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i-1,j,k,iTl)
               end do
               call sc%get_bcond('ym',my_bc)
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j+1,k,iTl)
               end do
               call sc%get_bcond('yp',my_bc)
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j-1,k,iTl)
               end do
         
            end do

            ! One-field temperature
            T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)

         end block advance_scalar


         ! Analytical
         ! advance_scalar: block
         !    use tpscalar_class,  only: bcond
         !    type(bcond), pointer :: my_bc
         !    integer  :: i,j,k,nsc,n,p
         !    real(WP) :: dt_sc

         !    ! Perform scalar time integration
         !    do while (timeSC%t.lt.time%t)
               
         !       ! Increment scalar time step
         !       if (timeSC%t+timeSC%dt.gt.time%t) then
         !          dt_sc=timeSC%dt
         !          timeSC%dt=time%t-timeSC%t
         !          call timeSC%increment()
         !          timeSC%dt=dt_sc
         !       else
         !          call timeSC%increment()
         !       end if
               
         !       ! Remember old scalar
         !       sc%SCold =sc%SC

         !       ! Analytical temperature solutions
         !       do k=cfg%kmin_,cfg%kmax_
         !          do j=cfg%jmin_,cfg%jmax_
         !             do i=cfg%imin_,cfg%imax_
         !                if (vf%VF(i,j,k).eq.0.0_WP) then
         !                   sc%SC(i,j,k,iTl)=0.0_WP
         !                   sc%SC(i,j,k,iTg)=T_sat
         !                else if (vf%VF(i,j,k).eq.1.0_WP) then
         !                   R_ext=get_Rext(time%t)
         !                   r=sqrt(cfg%xm(i)**2+cfg%ym(j)**2)
         !                   ! sc%SC(i,j,k,iTl)=T_inf-2.0_WP*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)*Simpson_adpt(beta_int,beta,1.0_WP-R_ext/r,1.0_WP,real(1e-7,WP))
         !                   sc%SC(i,j,k,iTl)=T_inf-2.0_WP*beta**2*(rho_g*(h_lg+(Cp_l-Cp_g)*(T_inf-T_sat)))/(rho_l*Cp_l)*Simpson(beta_int,beta,1.0_WP-R_ext/r,1.0_WP,20)
         !                   sc%SC(i,j,k,iTg)=0.0_WP
         !                else
         !                   sc%SC(i,j,k,iTl)=T_sat
         !                   sc%SC(i,j,k,iTg)=T_sat
         !                end if
         !             end do
         !          end do
         !       end do

         !       ! Apply boundary conditions
         !       where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
         !          sc%SC(:,:,:,iTl)=T_sat
         !          sc%SC(:,:,:,iTg)=T_sat
         !       end where
         !       call sc%apply_bcond(timeSC%t,timeSC%dt)
         !       call sc%get_bcond('xm',my_bc)
         !       do n=1,my_bc%itr%no_
         !          i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !          sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i+1,j,k,iTl)
         !       end do
         !       call sc%get_bcond('xp',my_bc)
         !       do n=1,my_bc%itr%no_
         !          i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !          sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i-1,j,k,iTl)
         !       end do
         !       call sc%get_bcond('ym',my_bc)
         !       do n=1,my_bc%itr%no_
         !          i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !          sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j+1,k,iTl)
         !       end do
         !       call sc%get_bcond('yp',my_bc)
         !       do n=1,my_bc%itr%no_
         !          i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !          sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j-1,k,iTl)
         !       end do
         
         !    end do

         !    ! One-field temperature
         !    T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)

         ! end block advance_scalar


         ! ! RK2
         ! advance_scalar: block
         !    use tpscalar_class,  only: bcond
         !    type(bcond), pointer :: my_bc
         !    integer  :: i,j,k,nsc,n,p
         !    real(WP) :: dt_sc

         !    ! Perform scalar time integration
         !    do while (timeSC%t.lt.time%t)
               
         !       ! Increment scalar time step
         !       if (timeSC%t+timeSC%dt.gt.time%t) then
         !          dt_sc=timeSC%dt
         !          timeSC%dt=time%t-timeSC%t
         !          call timeSC%increment()
         !          timeSC%dt=dt_sc
         !       else
         !          call timeSC%increment()
         !       end if
               
         !       ! Remember old scalar
         !       sc%SCold =sc%SC

         !       ! Build mid-time scalar
         !       do nsc=1,sc%nscalar
         !          p=sc%phase(nsc)
         !          do k=cfg%kmino_,cfg%kmaxo_
         !             do j=cfg%jmino_,cfg%jmaxo_
         !                do i=cfg%imino_,cfg%imaxo_
         !                   if (sc%PVF(i,j,k,p).eq.0.0_WP) then
         !                      sc%SC(i,j,k,nsc)=0.0_WP
         !                   else if (sc%PVF(i,j,k,p).lt.1.0_WP) then
         !                      sc%SC(i,j,k,nsc)=T_sat
         !                   else
         !                   end if
         !                end do
         !             end do
         !          end do
         !       end do

         !       ! Explicit calculation of dSC/dt
         !       call sc%get_dSCdt(dSCdt=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old)

         !       ! Assemble explicit residual
         !       resSC=resSC*timeSC%dt
         !       do nsc=1,sc%nscalar
         !          p=sc%phase(nsc)
         !          do k=cfg%kmin_,cfg%kmax_
         !             do j=cfg%jmin_,cfg%jmax_
         !                do i=cfg%imin_,cfg%imax_
         !                   if (sc%PVF(i,j,k,p).lt.1.0_WP) then
         !                      resSC(i,j,k,nsc)=0.0_WP
         !                   end if
         !                end do
         !             end do
         !          end do
         !       end do

         !       sc%SC=sc%SCold+0.5_WP*resSC

         !       call sc%get_dSCdt(dSCdt=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old)

         !       ! Assemble explicit residual
         !       resSC=resSC*timeSC%dt
         !       do nsc=1,sc%nscalar
         !          p=sc%phase(nsc)
         !          do k=cfg%kmin_,cfg%kmax_
         !             do j=cfg%jmin_,cfg%jmax_
         !                do i=cfg%imin_,cfg%imax_
         !                   if (sc%PVF(i,j,k,p).lt.1.0_WP) then
         !                      resSC(i,j,k,nsc)=0.0_WP
         !                   end if
         !                end do
         !             end do
         !          end do
         !       end do

         !       ! Apply the residuals
         !       sc%SC=sc%SCold+resSC

         !       ! Apply boundary conditions
         !       where (vf%VF.gt.VFlo.and.vf%VF.lt.VFhi)
         !          sc%SC(:,:,:,iTl)=T_sat
         !          sc%SC(:,:,:,iTg)=T_sat
         !       end where
         !       call sc%apply_bcond(timeSC%t,timeSC%dt)
         !       call sc%get_bcond('xm',my_bc)
         !       do n=1,my_bc%itr%no_
         !          i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !          sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i+1,j,k,iTl)
         !       end do
         !       call sc%get_bcond('xp',my_bc)
         !       do n=1,my_bc%itr%no_
         !          i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !          sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i-1,j,k,iTl)
         !       end do
         !       call sc%get_bcond('ym',my_bc)
         !       do n=1,my_bc%itr%no_
         !          i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !          sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j+1,k,iTl)
         !       end do
         !       call sc%get_bcond('yp',my_bc)
         !       do n=1,my_bc%itr%no_
         !          i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
         !          sc%SC(i,j,k,iTl)=2.0_WP*T_inf-sc%SC(i,j-1,k,iTl)
         !       end do
         
         !    end do

         !    ! One-field temperature
         !    T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)

         ! end block advance_scalar


         ! ================== PHASE CHANGE ================== !

         advance_lgpc: block
            integer :: i,j,k
            ! Get temperature gradient
            call lg%get_temperature_grad()
            ! Phase change mass flux
            lg%mdot2p=0.0_WP
            do k=lg%cfg%kmin_,lg%cfg%kmax_
               do j=lg%cfg%jmin_,lg%cfg%jmax_
                  do i=lg%cfg%imin_,lg%cfg%imax_
                     if ((vf%VF(i,j,k).gt.VFlo).and.(vf%VF(i,j,k).lt.VFhi)) then
                        lg%mdot2p(i,j,k)=(-k_g*lg%Tg_grd(i,j,k)+k_l*lg%Tl_grd(i,j,k))/h_lg
                     end if
                  end do
               end do
            end do
            call lg%cfg%sync(lg%mdot2p)
            ! Get the volumetric phase change mass flux
            call lg%get_mdot3p()
            ! Shift the phase change mass flux
            call lg%shift_mdot3p()
            ! Get the phase change induced divergence
            call lg%get_div()
         end block advance_lgpc


         ! ================== VELOCITY ================== !

         advance_flow: block
            use tpns_class, only: static_contact,harmonic_visc

            ! Prepare old staggered density (at n)
            call fs%get_olddensity(vf=vf)

            ! Prepare new staggered viscosity (at n+1)
            call fs%get_viscosity(vf=vf,strat=harmonic_visc)
            
            ! Perform velocity sub-iterations
            do while (time%it.le.time%itmax)
               if (cfg%amRoot) print*,'Flow sub iteration ',time%it
               
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
               
               ! Apply boundary conditions
               call fs%apply_bcond(time%t,time%dt)
               call apply_dirichlet()
               
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
               
               ! Increment sub-iteration counter
               time%it=time%it+1
               
            end do
            
            ! Recompute interpolated velocity and divergence
            call fs%interp_vel(Ui,Vi,Wi)
            call fs%get_div(src=lg%div_vel)

            ! Debug
            call fs%get_mfr()
            call lg%cfg%integrate(lg%div_vel,mfr_err)
            mfr_err=abs(mfr_err-sum(fs%mfr))

         end block advance_flow
         

         ! ================== OUTPUT ================== !

         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Update bubble radius
         R=get_R()
         R_ext=get_Rext(time%t)

         ! Temperature profile
         if (T_evt%occurs()) call get_Tl()

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call sc%get_max()
         call mfile%write()
         call cflfile%write()
         call scfile%write()
         call lgfile%write()
         ! Debug
         call get_tlgrd()
         call tlgrd_file%write()
         
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