!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use ddadi_class,       only: ddadi
   use hypre_str_class,   only: hypre_str
   use lowmach_class,     only: lowmach
   use vdscalar_class,    only: vdscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single low Mach flow solver and scalar solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs
   type(ddadi),       public :: ss
   type(lowmach),     public :: fs
   type(vdscalar),    public :: sc
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resSC
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   !> Equation of state
   real(WP) :: Zst,rho0,rho1,rhost
   
   !> Initial mass
   real(WP) :: mass=0.0_WP
   
   !> Initial perturbations
   integer :: nwaveX,nwaveZ
   real(WP) :: wamp
   real(WP), dimension(:), allocatable :: wnumbX,wshiftX,wnumbZ,wshiftZ

   !> Quadrature rule
   integer :: nq
   real(WP), dimension(:), allocatable :: xq,wq
   
contains
   
   
   !> Obtain density from equation of state based on Burke-Schumann
   subroutine get_rho()
      implicit none
      integer  :: i,j,k
      real(WP) :: coeff,mySC
      !real(WP), dimension(3) :: grad
      integer :: iq,jq,kq
      ! Evaluate rho
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               ! Direct evaluation
               !fs%RHO(i,j,k)=rho(sc%SC(i,j,k))
               ! Quadrature rule
               fs%RHO(i,j,k)=0.0_WP
               !grad(1)=minmod((sc%SC(i+1,j,k)-sc%SC(i,j,k))*cfg%dxmi(i+1),(sc%SC(i,j,k)-sc%SC(i-1,j,k))*cfg%dxmi(i))
               !grad(2)=minmod((sc%SC(i,j+1,k)-sc%SC(i,j,k))*cfg%dymi(j+1),(sc%SC(i,j,k)-sc%SC(i,j-1,k))*cfg%dymi(j))
               !grad(3)=minmod((sc%SC(i,j,k+1)-sc%SC(i,j,k))*cfg%dzmi(k+1),(sc%SC(i,j,k)-sc%SC(i,j,k-1))*cfg%dzmi(k))
               do kq=1,nq
                  do jq=1,nq
                     do iq=1,nq
                        ! Use MUSCL to get scalar in the cell
                        !mySC=sc%SC(i,j,k)+grad(1)*xq(iq)*cfg%dx(i)+grad(2)*xq(jq)*cfg%dy(j)+grad(3)*xq(kq)*cfg%dz(k)
                        ! Use trilinear interpolation to get scalar in the cell
                        mySC=cfg%get_scalar([cfg%x(i)+xq(iq)*cfg%dx(i),cfg%y(j)+xq(jq)*cfg%dy(j),cfg%z(k)+xq(kq)*cfg%dz(k)],i,j,k,sc%SC,'n')
                        ! Accumulate integral
                        fs%RHO(i,j,k)=fs%RHO(i,j,k)+wq(iq)*wq(jq)*wq(kq)*rho(mySC)
                     end do
                  end do
               end do
            end do
         end do
      end do
      ! Synchronize and apply Neumann top and bottom
      call cfg%sync(fs%RHO)
      if (cfg%jproc.eq.1) then
         do j=cfg%jmino,cfg%jmin-1
            fs%RHO(:,j,:)=fs%RHO(:,cfg%jmin,:)
         end do
      end if
      if (cfg%jproc.eq.cfg%npy) then
         do j=cfg%jmax+1,cfg%jmaxo
            fs%RHO(:,j,:)=fs%RHO(:,cfg%jmax,:)
         end do
      end if
      ! Rescale rho to preserve its integral
      if (mass.eq.0.0_WP) then
         ! First time, just calculate mass
         call cfg%integrate(fs%RHO,integral=mass)
      else
         ! Mass is non-zero, rescale to conserve it
         call cfg%integrate(fs%RHO,integral=coeff); coeff=mass/coeff
         fs%RHO=fs%RHO*coeff
      end if
   contains
      real(WP) function rho(Z)
         implicit none
         real(WP), intent(in) :: Z
         real(WP) :: Zclip
         ! Clip mixture fraction between - and 1
         Zclip=min(max(Z,0.0_WP),1.0_WP)
         ! Pure mixing flamelet
         !rho=rho0*rho1/(rho1*(1.0_WP-Zclip)+rho0*Zclip)
         ! Burke-Schumann flamelet
         if (Zclip.le.Zst) then
            rho=Zst*rho0*rhost/(rhost*(Zst-Zclip)+rho0*Zclip)
         else
            rho=(1.0_WP-Zst)*rhost*rho1/(rho1*(1.0_WP-Zclip)+rhost*(Zclip-Zst))
         end if
      end function rho
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
   end subroutine get_rho
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Read in EOS parameters
      initialize_eos: block
         use mathtools, only: quadrature_rule
         ! EOS parameters
         call param_read('rho0',rho0)
         call param_read('rho1',rho1)
         call param_read('rhost',rhost)
         call param_read('Zst',Zst)
         if (Zst.le.0.0_WP) Zst=0.0_WP+epsilon(Zst)
         if (Zst.ge.1.0_WP) Zst=1.0_WP-epsilon(Zst)
         ! Prepare quadrature for integrating density
         nq=5; allocate(xq(nq),wq(nq))
         call quadrature_rule(nq,xq,wq)
      end block initialize_eos
      
      ! Initialize disturbances
      initialize_disturbances: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         use parallel,  only: MPI_REAL_WP
         use mpi_f08,   only: MPI_BCAST
         use string,    only: str_long
         use messager,  only: log
         character(str_long) :: message
         integer :: ierr,n
         ! Prepare interface disturbance in X
         call param_read('Wave amplitude',wamp)
         call param_read('NwaveX',NwaveX,default=-1)
         if (NwaveX.gt.0) then
            allocate(wnumbX(nwaveX),wshiftX(nwaveX))
            call param_read('wnumbX',wnumbX)
            call param_read('wshiftX',wshiftX)
         else
            nwaveX=6
            allocate(wnumbX(nwaveX),wshiftX(nwaveX))
            wnumbX=[3.0_WP,4.0_WP,5.0_WP,6.0_WP,7.0_WP,8.0_WP]*twoPi/cfg%xL
            if (cfg%amRoot) then
               do n=1,nwaveX
                  wshiftX(n)=random_uniform(lo=-0.5_WP*cfg%xL,hi=+0.5_WP*cfg%xL)
               end do
            end if
            call MPI_BCAST(wshiftX,nwaveX,MPI_REAL_WP,0,cfg%comm,ierr)
         end if
         ! Prepare interface disturbance in Z
         call param_read('NwaveZ',NwaveZ,default=-1)
         if (NwaveZ.gt.0) then
            allocate(wnumbZ(nwaveZ),wshiftZ(nwaveZ))
            call param_read('wnumbZ',wnumbZ)
            call param_read('wshiftZ',wshiftZ)
         else
            nwaveZ=6
            allocate(wnumbZ(nwaveZ),wshiftZ(nwaveZ))
            wnumbZ=[3.0_WP,4.0_WP,5.0_WP,6.0_WP,7.0_WP,8.0_WP]*twoPi/cfg%zL
            ! Handle 2D case
            if (cfg%nz.eq.1) wnumbZ=0.0_WP
            if (cfg%amRoot) then
               do n=1,nwaveZ
                  wshiftZ(n)=random_uniform(lo=-0.5_WP*cfg%zL,hi=+0.5_WP*cfg%zL)
               end do
            end if
            call MPI_BCAST(wshiftZ,nwaveZ,MPI_REAL_WP,0,cfg%comm,ierr)
         end if
         ! Print out initial disturbance
         if (cfg%amRoot) then
            write(message,'("[Initial conditions] =>  NwaveX =",i6)')              nwaveX; call log(message)
            write(message,'("[Initial conditions] =>  wnumbX =",1000(es12.5,x))')  wnumbX; call log(message)
            write(message,'("[Initial conditions] => wshiftX =",1000(es12.5,x))') wshiftX; call log(message)
            write(message,'("[Initial conditions] =>  NwaveZ =",i6)')              nwaveZ; call log(message)
            write(message,'("[Initial conditions] =>  wnumbZ =",1000(es12.5,x))')  wnumbZ; call log(message)
            write(message,'("[Initial conditions] => wshiftZ =",1000(es12.5,x))') wshiftZ; call log(message)
         end if
      end block initialize_disturbances
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
      end block initialize_timetracker
      
      ! Create a scalar solver
      create_scalar: block
         use vdscalar_class, only: dirichlet,neumann
         real(WP) :: diffusivity
         ! Create scalar solver
         call sc%initialize(cfg=cfg,name='MixFrac')
         ! Add slight backward bias to CN scheme
         sc%theta=sc%theta!+1.0e-2_WP
         ! Define boundary conditions
         call sc%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         call sc%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         ! Assign constant diffusivity
         call param_read('Dynamic diffusivity',diffusivity)
         sc%diff=diffusivity
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Scalar',nst=13)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
      end block create_scalar
      
      ! Create a low-Mach flow solver with bconds
      create_velocity_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use lowmach_class,   only: slip
         real(WP) :: visc
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Variable density low Mach NS')
         ! Add slight backward bias to CN scheme
         fs%theta=fs%theta+1.0e-2_WP
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define boundary conditions
         call fs%add_bcond(name='yp',type=slip,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
         call fs%add_bcond(name='ym',type=slip,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=18
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_velocity_solver
      
      ! Initialize our mixture fraction field
      initialize_scalar: block
         integer :: i,j,k,nX,nZ
         real(WP) :: y
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Distance to the nominal interface position
                  y=cfg%ym(j)
                  ! Perturb interface position
                  do nX=1,nwaveX
                     do nZ=1,nwaveZ
                        y=y+wamp*cos(wnumbX(nX)*(cfg%xm(i)-wshiftX(nX)))*cos(wnumbZ(nZ)*(cfg%zm(k)-wshiftZ(nZ)))
                     end do
                  end do
                  ! Hyperbolic tangent
                  sc%SC(i,j,k)=0.5_WP*(1.0_WP+tanh(y))
               end do
            end do
         end do
      end block initialize_scalar
      
      ! Initialize our velocity field
      initialize_velocity: block
         integer :: i,j,k,nX,nZ
         real(WP) :: y
         ! Initialize density from scalar
         call get_rho(); call fs%update_sRHO()
         fs%RHOold=fs%RHO; fs%sRHOXold=fs%sRHOX; fs%sRHOYold=fs%sRHOY; fs%sRHOZold=fs%sRHOZ
         ! Initialize velocity field
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Distance to the nominal interface position
                  y=cfg%ym(j)
                  ! Perturb interface position
                  do nX=1,nwaveX
                     do nZ=1,nwaveZ
                        y=y+wamp*cos(wnumbX(nX)*(cfg%x(i)-wshiftX(nX)))*cos(wnumbZ(nZ)*(cfg%zm(k)-wshiftZ(nZ)))
                     end do
                  end do
                  ! Hyperbolic tangent
                  fs%U(i,j,k)=tanh(y)
               end do
            end do
         end do
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Get Umid/Vmid/Wmid
         call fs%get_Umid()
         ! Calculate cell-centered velocities and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(dt=time%dt)
      end block initialize_velocity
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='mixing')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('density',fs%rho)
         call ens_out%add_scalar('mixfrac',sc%SC)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call sc%get_max(fs%RHO)
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
         call mfile%add_column(sc%SCmax,'SCmax')
         call mfile%add_column(sc%SCmin,'SCmin')
         call mfile%add_column(sc%rhoSCint,'rhoSCint')
         call mfile%add_column(fs%RHOmax,'RHOmax')
         call mfile%add_column(fs%RHOmin,'RHOmin')
         call mfile%add_column(fs%RHOint,'RHOint')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
   contains
      
      !> Function that localizes y- boundary
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin) isIn=.true.
      end function ym_locator
      
      !> Function that localizes y- boundary
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
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old scalar, velocity, and density
         sc%SCold=sc%SC
         fs%RHOold=fs%RHO
         fs%sRHOXold=fs%sRHOX
         fs%sRHOYold=fs%sRHOY
         fs%sRHOZold=fs%sRHOZ
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! ============= SCALAR SOLVER =======================
            ! Explicit calculation of drhoSC/dt from scalar equation
            call sc%get_drhoSCdt(resSC,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Assemble explicit residual
            resSC=time%dt*resSC-(fs%RHO*sc%SC-fs%RHOold*sc%SCold)
            
            ! Form implicit residual
            call sc%solve_implicit(time%dt,resSC,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Apply this residual
            sc%SC=sc%SC+resSC
            
            ! Apply other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! ============ UPDATE DENSITY========================
            call get_rho(); call fs%update_sRHO()
            ! ===================================================
            
            ! ============ VELOCITY SOLVER ======================
            ! Compute Umid
            call fs%get_Umid()
            
            ! Get rhoU/rhoV/rhoW
            call fs%rho_multiply()
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=time%dt*resU-(fs%sRHOX**2*fs%U-fs%sRHOXold**2*fs%Uold)
            resV=time%dt*resV-(fs%sRHOY**2*fs%V-fs%sRHOYold**2*fs%Vold)
            resW=time%dt*resW-(fs%sRHOZ**2*fs%W-fs%sRHOZold**2*fs%Wold)
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=fs%U+resU
            fs%V=fs%V+resV
            fs%W=fs%W+resW
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! ============ PRESSURE SOLVER ======================
            ! Compute Umid
            call fs%get_Umid()
            ! Get rhoU/rhoV/rhoW
            call fs%rho_multiply()
            ! Correct MFR
            call fs%correct_mfr(dt=time%dt)
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%get_div(dt=time%dt)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dt*resU*((1.0_WP-fs%theta)*fs%sRHOXold**2+fs%theta*fs%sRHOX**2)/((fs%sRHOX+fs%sRHOXold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOX)
            fs%rhoV=fs%rhoV-time%dt*resV*((1.0_WP-fs%theta)*fs%sRHOYold**2+fs%theta*fs%sRHOY**2)/((fs%sRHOY+fs%sRHOYold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOY)
            fs%rhoW=fs%rhoW-time%dt*resW*((1.0_WP-fs%theta)*fs%sRHOZold**2+fs%theta*fs%sRHOZ**2)/((fs%sRHOZ+fs%sRHOZold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOZ)
            ! Recover Umid
            call fs%rho_divide()
            ! Recover U
            call fs%get_U()
            ! ===================================================
            
            ! Increment sub-iteration
            time%it=time%it+1
            
         end do
         
         ! Could consider doing final scalar step here
         !call sc%get_drhoSCdt(resSC,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
         !resSC=time%dt*resSC-(fs%RHO*sc%SC-fs%RHOold*sc%SCold)
         !call sc%solve_implicit(time%dt,resSC,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
         !sc%SC=sc%SC+resSC
         !call sc%apply_bcond(time%t,time%dt)
         
         ! Recompute interpolated velocity and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(dt=time%dt)
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call sc%get_max(fs%RHO)
         call mfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(resSC,resU,resV,resW,Ui,Vi,Wi)
   end subroutine simulation_final
   
   
end module simulation
