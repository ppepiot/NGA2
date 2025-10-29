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
   use string,            only: str_short,str_medium
   use YAMLRead,          only: YAMLElement
   use chem_sys_class,    only: chem_sys
   use chem_state_class,  only: chem_state,fixed_PH
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
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

   !> The array of the species. Eeach stored as a YAMLElement object
   type(YAMLElement), dimension(:), allocatable :: species

   !> Species names
   character(len=str_medium), dimension(:), allocatable :: sp_names

   !> Chemical system and state
   type(chem_sys)   :: sys
   type(chem_state) :: state
   
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
   real(WP), dimension(:),       allocatable :: MM,sc_itf
   ! Debug
   real(WP), dimension(:,:,:),   allocatable :: dbg_flg
   integer :: i_itf
   real(WP) :: mfr_err
   
   !> Problem definition
   real(WP) :: H0
   integer  :: iWv,iWl,iO2,iN2,iTl,iTg
   real(WP) :: rho_l,rho_g,k_l,k_g,Cp_l,Cp_g,alpha_l,alpha_g,h_lg,T_sat,T_g
   real(WP) :: T_w,beta,t0
   real(WP) :: x_itf,u_itf,x_itf_ext,u_itf_ext
   real(WP) :: prhs_int
   real(WP) :: pressure,air2wv_rat,N2O_rat,YO2,YN2
   integer  :: ns,np=2
   
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
      use irl_fortran_interface, only: calculateCentroid
      integer :: i,j,k
      real(WP), dimension(3) :: posI
      ! Get the interface location
      i=vf%band_map(1,1)
      j=vf%band_map(2,1)
      k=vf%band_map(3,1)
      posI=calculateCentroid(vf%interface_polygon(1,i,j,k))
      x_itf=posI(1)
      i_itf=i
      ! Get the interface velocity
      u_itf=fs%U(i+1,j,k)
      ! Get the analytical solution
      x_itf_ext=2.0_WP*beta*sqrt(alpha_g*(time%t))
      u_itf_ext=beta*sqrt(alpha_g/(time%t))
   end subroutine get_interface


   !> Function that returns the index of an input species name
   function get_sp_ind(name)
      use messager, only: die
      implicit none
      character(len=*), intent(in) :: name
      integer :: isc,get_sp_ind
      do isc=1,ns
         if (trim(name).eq.trim(sp_names(isc))) then
            get_sp_ind=isc
            return
         end if
      end do
      call die('[water_drop_2D_ceq get_sp_ind] Unknown species')
   end function get_sp_ind


   subroutine interface_jump()
      use messager, only: die
      implicit none
      real(WP), dimension(:),     allocatable :: vol_new,vol_old,mp,N,phasicHoR,Y
      logical,  dimension(:,:,:), allocatable :: clustered
      logical,  dimension(:),     allocatable :: active
      integer,  dimension(:,:),   allocatable :: cell_indices
      real(WP), dimension(:),     allocatable :: Vscaled,vof_old,vof_new,w
      real(WP) :: Vnew,Vold,Nsum,vof,itf_area,Tl,Tg
      integer  :: i,j,k,index,isc,p,n_clustered,m
      integer  :: in,jn,kn
      integer  :: stx,sty,stz
      real(WP) :: mdot2p
      real(WP), parameter :: trsh_lo=0.05_WP,trsh_hi=1.0_WP-trsh_lo,wmin=1.0e-16_WP,dVlmin=1.0e-16_WP
      integer,  parameter :: nc_max=27
      real(WP) :: dVl,dVl_i,dVl_rem,Vref,vof_tmp,interfaceness,wsum

      ! Debug
      dbg_flg=0.0_WP

      Vref=minval(cfg%vol)

      ! Zero out the phase change mass flux
      lg%mdot2p=0.0_WP

      ! Allocate arrays
      allocate(vol_new(Lphase:Gphase))
      allocate(vol_old(Lphase:Gphase))
      allocate(mp(Lphase:Gphase))
      allocate(N(ns))
      allocate(phasicHoR(Lphase:Gphase))
      allocate(Y(ns))
      allocate(clustered(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); clustered=.false.
      allocate(cell_indices(3,nc_max)); cell_indices=0
      allocate(Vscaled(nc_max)); Vscaled=0.0_WP
      allocate(vof_old(nc_max)); vof_old=0.0_WP
      allocate(vof_new(nc_max)); vof_new=0.0_WP
      allocate(w(nc_max)); w=0.0_WP
      allocate(active(nc_max)); active=.false.

      ! Clustering stencil
      if (cfg%nx.gt.1) then
         stx=1
      else
         stx=0
      end if
      if (cfg%ny.gt.1) then
         sty=1
      else
         sty=0
      end if
      if (cfg%nz.gt.1) then
         stz=1
      else
         stz=0
      end if

      ! Loop over the interfacial cells
      do index=1,vf%band_count(0)

         ! Get the interfacial cell indices
         i=vf%band_map(1,index)
         j=vf%band_map(2,index)
         k=vf%band_map(3,index)

         ! Skip if already clustered
         if (clustered(i,j,k)) cycle
         
         ! Add current cell to the potential cluster
         n_clustered=1
         cell_indices(:,1)=[i,j,k]

         ! Initialize the interfacial area and old volumes
         itf_area=cfg%vol(i,j,k)*vf%SD(i,j,k)
         vol_old=sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) then
            print*,'before jump:'
            print*,'Index of liquid and gas = ',Lphase,Gphase
            print*,'Density of liquid and gas = ',sc%Prho
            print*,'VOF(i_itf) = ',vf%VF(i,j,k)
            print*,'Volume of liquid and gase = ',vol_old
            print*,'sc_itf = ',sc_itf
         end if

         ! Pre-evaluate the equilibrium
         mp=sc%Prho*sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
         Y =sc%SC(i,j,k,1:ns)
         Tl=sc%SC(i,j,k,iTl)
         Tg=sc%SC(i,j,k,iTg)
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) then
            print*,'Initial mass = ',mp
            print*,'Y = ',Y
            print*,'Tl = ',Tl
            print*,'Tg = ',Tg
         end if
         call get_equilibrium()
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) then
            print*,'state%success = ',state%success
            ! call die('')
         end if
         ! Cluster cells if the equilibrium failed
         if (.not.state%success) then

            ! Mark it as clustered
            clustered(i,j,k)=.true.

            ! Initialize the species mass
            do isc=1,ns
               p=sc%phase(isc)
               Y(isc)=sc%Prho(p)*sc%PVF(i,j,k,p)*cfg%vol(i,j,k)*sc%SC(i,j,k,isc)
            end do

            ! Initialize the mass-averaged temperatures
            Tl=sc%Prho(Lphase)*sc%PVF(i,j,k,Lphase)*cfg%vol(i,j,k)*sc%SC(i,j,k,iTl)
            Tg=sc%Prho(Gphase)*sc%PVF(i,j,k,Gphase)*cfg%vol(i,j,k)*sc%SC(i,j,k,iTg)

            ! Loop over the cluster stencil skipping the ghost cells
            z_loop: do kn=k-stz,k+stz
               if (kn.lt.cfg%kmin_.or.kn.gt.cfg%kmax_) cycle
               y_loop: do jn=j-sty,j+sty
                  if (jn.lt.cfg%jmin_.or.jn.gt.cfg%jmax_) cycle
                  x_loop: do in=i-stx,i+stx
                     if (in.lt.cfg%imin_.or.in.gt.cfg%imax_) cycle

                     ! Neighbor must be interfacial and not clustered yet
                     if (vf%VF(in,jn,kn).gt.VFlo.and.vf%VF(in,jn,kn).lt.VFhi.and..not.clustered(in,jn,kn)) then

                        ! Mark it as clustered
                        n_clustered=n_clustered+1
                        cell_indices(:,n_clustered)=[in,jn,kn]
                        clustered(in,jn,kn)=.true.

                        ! Accumulate old volumes
                        vol_old=vol_old+sc%PVF(in,jn,kn,:)*cfg%vol(in,jn,kn)

                        ! Accumulate mass*SC and mass*temperature
                        do isc=1,ns
                           p=sc%phase(isc)
                           Y(isc)=Y(isc)+sc%Prho(p)*sc%PVF(in,jn,kn,p)*cfg%vol(in,jn,kn)*sc%SC(in,jn,kn,isc)
                        end do
                        Tl=Tl+sc%Prho(Lphase)*sc%PVF(in,jn,kn,Lphase)*cfg%vol(in,jn,kn)*sc%SC(in,jn,kn,iTl)
                        Tg=Tg+sc%Prho(Gphase)*sc%PVF(in,jn,kn,Gphase)*cfg%vol(in,jn,kn)*sc%SC(in,jn,kn,iTg)

                        ! Accumulate interface area
                        itf_area=itf_area+cfg%vol(in,jn,kn)*vf%SD(in,jn,kn)
                     end if

                  end do x_loop
               end do y_loop
            end do z_loop

            ! Cluster-level phase masses
            mp=sc%Prho*vol_old

            ! Cluster-averaged mass fractions and temperatures
            do isc=1,ns
               Y(isc)=Y(isc)/mp(sc%phase(isc))
            end do
            Tl=Tl/mp(Lphase)
            Tg=Tg/mp(Gphase)

            ! Get the equilibrium state of the cluster
            call get_equilibrium()
            
         end if

         ! Calculate the cluster VOF
         vof=vol_old(Lphase)/sum(vol_old)
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) then
            print*,'cluster VOF = ',vof
            print*,'Initial volume = ',vol_old
         end if

         ! Store the old total volume
         Vold=sum(vol_old)

         ! Debug
         if (.not.state%success) then
            print*,'Cluster VOF = ',vof
            print*,'N initial scaled and fed into ceq = ',N
            print*,'N initial actual = ',N*Nsum
            print*,'HoR = ',sum(phasicHoR)
            print*,'T_g = ',T_g
            print*,'Clustered cells info:'
            print*,'n_clustered = ',n_clustered
            do m=1,n_clustered
               i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
               print*,'i,j,k = ',i,j,k
               print*,'VOF = ',vf%VF(i,j,k)
               print*,'SD = ',vf%SD(i,j,k)
            end do
            call die('line 302')
         end if

         ! Update the phase masses
         mp=0.0_WP
         do isc=1,ns
            p=sc%phase(isc)
            mp(p)=mp(p)+N(isc)*MM(isc)
         end do
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) print*,'Equilibrium phase masses = ',mp

         ! Get the phase volumes
         vol_new=mp/sc%Prho
         Vnew=sum(vol_new)

         ! Get the phase change mass flux
         mdot2p=(Vnew-Vold)/(time%dt*(1.0_WP/sc%Prho(Gphase)-1.0_WP/sc%Prho(Lphase))*itf_area)

         ! Gather geometry and current VOF per clustered cell
         do m=1,n_clustered
            i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
            active(m)=.true.
            Vscaled(m)=cfg%vol(i,j,k)/Vref ! Scale it for more accurate calculations
            vof_old(m)=vf%VF(i,j,k)
            vof_new(m)=vof_old(m)
         end do

         ! Total liquid volume change ( > 0 condensation, < 0 vaporization)
         dVl=(vol_new(Lphase)-vol_old(Lphase))/Vref ! Scale it for more accurate calculations
         dVl_rem=dVl

         if (abs(dVl).gt.dVlmin) then

            ! Build weights
            do m=1,n_clustered
               i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
               itf_area=cfg%vol(i,j,k)*vf%SD(i,j,k)
               interfaceness=minval(sc%PVF(i,j,k,:))
               w(m)=max(itf_area*interfaceness,wmin)
            end do

            ! Iteratively redistribute liquid
            do

               ! Update weights sum
               wsum=0.0_WP
               do m=1,n_clustered
                  if (active(m)) wsum=wsum+w(m)
               end do

               ! Terminate if succssesd
               if (abs(dVl_rem).le.dVlmin.or.wsum.le.0.0_WP) exit

               ! Distribute
               dVl=dVl_rem
               do m=1,n_clustered
                  if (.not.active(m)) cycle

                  ! Estimate VOF
                  dVl_i =(w(m)/wsum)*dVl
                  vof_tmp=vof_new(m)+dVl_i/Vscaled(m)

                  ! Clip it
                  if (vof_tmp.gt.1.0_WP) then
                     dVl_i=(1.0_WP-vof_new(m))*Vscaled(m)
                     vof_new(m)=1.0_WP
                     active(m)=.false.
                  else if (vof_tmp.lt.0.0_WP) then
                     dVl_i=(0.0_WP-vof_new(m))*Vscaled(m)
                     vof_new(m)=0.0_WP
                     active(m)=.false.
                  else
                     vof_new(m)=vof_tmp
                  end if

                  ! Correct the liquid volume change
                  dVl_rem=dVl_rem-dVl_i

               end do

            end do

         end if

         ! Assign per-cell fields (Need to treat cells with VOF=0 and 1, differently)
         do m=1,n_clustered

            ! Get the cell indices
            i=cell_indices(1,m)
            j=cell_indices(2,m)
            k=cell_indices(3,m)

            ! Assign VOF
            if (vof_new(m).lt.VFlo) then
               vf%VF(i,j,k)=0.0_WP
            else if (vof_new(m).gt.VFhi) then
               vf%VF(i,j,k)=1.0_WP
            else
               vf%VF(i,j,k)=vof_new(m)
            end if
            sc%PVF(i,j,k,Lphase)=vf%VF(i,j,k)
            sc%PVF(i,j,k,Gphase)=1.0_WP-vf%VF(i,j,k)

            ! Composition (Assuming the same mass fraction for all non-empty cells in the cluster)
            do isc=1,ns
               p=sc%phase(isc)
               if(sc%PVF(i,j,k,p).gt.0.0_WP) then
                  sc%SC(i,j,k,isc)=MM(isc)*N(isc)/mp(sc%phase(isc))
               else
                  sc%SC(i,j,k,isc)=0.0_WP
               end if
            end do

            ! Temperature and phase change mass flux
            if (vf%VF(i,j,k).eq.1.0_WP) then
               sc%SC(i,j,k,iTl)=state%T ! Not sure if this is good enough. Should I do linear interpolation instead?
               sc%SC(i,j,k,iTg)=0.0_WP
               lg%mdot2p(i,j,k)=0.0_WP
            else if (vf%VF(i,j,k).eq.0.0_WP) then
               sc%SC(i,j,k,iTl)=0.0_WP
               sc%SC(i,j,k,iTg)=state%T ! Not sure if this is good enough. Should I do linear interpolation instead?
               lg%mdot2p(i,j,k)=0.0_WP
            else
               sc%SC(i,j,k,iTl)=state%T
               sc%SC(i,j,k,iTg)=state%T
               lg%mdot2p(i,j,k)=mdot2p
            end if

         end do

      end do

      ! Update the interface (Do I need it? I don't think so)
      ! call vf%advect_interface(0.0_WP,fs%U,fs%V,fs%W)

      ! Remove flotsams and thin structures if needed
      call vf%remove_flotsams()
      call vf%remove_thinstruct()
      
      ! Synchronize and clean-up barycenter fields
      call vf%sync_and_clean_barycenters()
      
      ! Update the interface band
      call vf%update_band()
      
      ! Perform interface reconstruction from transported moments
      call vf%build_interface()
      
      ! Create discontinuous polygon mesh from IRL interface
      call vf%polygonalize_interface()
      
      ! Perform interface sensing (Do I need it?)
      ! if (vf%two_planes) call vf%sense_interface()
      
      ! Calculate distance from polygons (I don't think it's needed anywhere)
      ! call vf%distance_from_polygon()
      
      ! Calculate subcell phasic volumes (I don't think it's needed anywhere)
      ! call vf%subcell_vol()
      
      ! Calculate curvature
      call vf%get_curvature()
      
      ! Reset moments to guarantee compatibility with interface reconstruction
      call vf%reset_moments()

      ! Sync fields
      do isc=1,sc%nscalar
         call cfg%sync(sc%SC(:,:,:,isc))
      end do
      call cfg%sync(vf%VF)
      call cfg%sync(lg%mdot2p)

      ! Apply boundary conditions
      call sc%apply_bcond(time%t,time%dt)
      call vf%apply_bcond(time%t,time%dt)

      ! Debug
      where (clustered) dbg_flg=1.0_WP
      call cfg%sync(dbg_flg)

      ! Deallocate arrays
      deallocate(vol_new,vol_old,mp,N,phasicHoR,Y,clustered,cell_indices,Vscaled,vof_old,vof_new,w,active)

      contains

      subroutine get_equilibrium()
         implicit none
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) print*,'Inside get_equilibrium'
         ! Calculate and normalize the mole numbers

         do isc=1,ns
            N(isc)=Y(isc)*mp(sc%phase(isc))/MM(isc)
         end do
         Nsum=sum(N)
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) then
            print*,'Y = ',Y
            print*,'N before scaling up = ',N
         end if
         if (Nsum.gt.0.0_WP) N=N/Nsum
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) print*,'N after scaling up = ',N
         ! Get the phasic enthalpies
         call state%get_phasic_HoR(Lphase,N,Tl,phasicHoR(Lphase))
         call state%get_phasic_HoR(Gphase,N,Tg,phasicHoR(Gphase))
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) then
            print*,'Tl = ',Tl
            print*,'Tg = ',Tg
            print*,'phasicHoR = ',phasicHoR
            print*,'sum(phasicHoR) = ',sum(phasicHoR)
         end if
         ! Reinitialize the mole numbers
         call state%N_init(N=N,HoR=sum(phasicHoR),T_g=T_g)
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) then
            print*,'The computed enthalpy of the mixture cell = ',sum(phasicHoR)
            ! print*,'The enthalpy that corresponds to 373.15 K = ',-31764
         end if
         ! call state%N_init(N=N,HoR=-31764.252928818278_WP,T_g=T_g)
         ! call state%N_init(T=T_sat,N=N)
         if (.not.state%success) then
            print*,'Cluster VOF = ',vof
            print*,'N = ',N
            print*,'N*Nsum = ',N*Nsum
            print*,'HoR = ',sum(phasicHoR)
            print*,'T_g = ',T_g
            print*,'Clustered cells info:'
            print*,'n_clustered = ',n_clustered
            do m=1,n_clustered
               i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
               print*,'i,j,k = ',i,j,k
               print*,'VOF = ',vf%VF(i,j,k)
               print*,'SD = ',vf%SD(i,j,k)
            end do
            call die('line 517')
         end if
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) then
            print*,'Initialized N = ',state%N
            print*,'Initialized Ndu = ',state%Ndu
         end if
         ! Get the chemical equilibrium
         call state%equilibrate()
         if (i.eq.i_itf.and.j.eq.1.and.k.eq.1) then
            print*,'Equilibrium N = ',state%N
            print*,'Equilibrium T = ',state%T
         end if

         ! Re-scale the mole numbers
         ! N=state%N*Nsum
         ! Debug: If not successful, don't assign N so I know what initial moles caused this
         if (state%success) N=state%N*Nsum

      end subroutine get_equilibrium

   end subroutine interface_jump


   ! Linear interpolation of the scalar at the interfacial cell center using the interfacial value and the adjacent cells
   ! Right now only works for this stephan case (1D + Interface goes from left to right)
   ! Requires the most recent interfacial values
   subroutine apply_itf_bcond()
      integer :: isc,p,i,j,k,ishift,in
      real    :: w_n,w_i,w_sum
      ! Update the interface location
      call get_interface()
      ! Loop over the scalars
      do isc=1,sc%nscalar
         p=sc%phase(isc)
         ishift=1-2*p
         in=i_itf+ishift
         ! Prepare the interpolation weights
         w_n=x_itf-cfg%xm(i_itf)
         w_i=cfg%xm(i_itf)-cfg%xm(in)
         w_sum=x_itf-cfg%xm(in)
         w_n=w_n/w_sum
         w_i=w_i/w_sum
         ! Loop over the cells
         do k=sc%cfg%kmin_,sc%cfg%kmax_
            do j=sc%cfg%jmin_,sc%cfg%jmax_
               do i=sc%cfg%imin_,sc%cfg%imax_
                  if (sc%PVF(i,j,k,p).eq.0.0_WP) then
                     sc%SC(i,j,k,isc)=0.0_WP
                  else if (sc%PVF(i,j,k,p).lt.1.0_WP) then
                     sc%SC(i,j,k,isc)=w_n*sc%SC(in,j,k,isc)+w_i*sc_itf(isc)
                  else
                  end if
               end do
            end do
         end do
         ! Sync it
         call cfg%sync(sc%SC(:,:,:,isc))
      end do
   end subroutine apply_itf_bcond


   !> Initialization of problem solver
   subroutine simulation_init
      use param,    only: param_exists,param_read,param_getsize
      use messager, only: die
      implicit none
      integer :: ne,ncs
      character(len=str_short), dimension(:), allocatable :: e_names
      real(WP), dimension(:,:), allocatable :: elem_mat
      real(WP), dimension(:,:), allocatable :: phse_mat
      real(WP), allocatable :: nasa_coef(:,:)
      character(len=str_medium), dimension(:), allocatable :: const_sp
      integer,  dimension(:), allocatable :: CS
      
      
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
         call param_read('Pressure',Pressure)
         alpha_l=k_l/(rho_l*Cp_l)
         alpha_g=k_g/(rho_g*Cp_g)
      end block read_inputs


      ! Parse the mechanism file
      parse_mech: block
         use chem_sys_class, only: ncof
         use YAMLRead,       only: YAMLHandler,YAMLSequence,YAMLMap,yaml_open_file,yaml_start_from_sequence,yaml_close_file
         character(len=str_medium) :: mch_file
         character(len=str_short), dimension(:), allocatable :: sp_names_copy,const_sp_copy
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
         if (param_exists('Constrained species')) then
            ncs=param_getsize('Constrained species')
            allocate(const_sp(1:ncs))
            allocate(const_sp_copy(1:ncs))
            call param_read('Constrained species',const_sp)
            const_sp_copy=const_sp
         else
            ncs=0
         end if
         allocate(sp_names(1:ns))
         allocate(sp_names_copy(1:ns))
         allocate(CS(ncs))
         call param_read('Species',sp_names)
         sp_names_copy=sp_names
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
         deallocate(sp_names_copy)
         if (allocated(const_sp_copy)) deallocate(const_sp_copy)
         print*,'sp_names = ',sp_names
      end block parse_mech


      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resSC (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:ns+2))
         allocate(resU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(T     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(MM(ns)); MM=[32.0_WP,18.0_WP,18.0_WP,28.0_WP]; MM=0.001_WP*MM
         ! Debug
         allocate(dbg_flg(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); dbg_flg=0.0_WP
      end block allocate_work_arrays
      
      
      ! Initialize the chemical equilibrium framework
      ceq_init: block
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
         print*,'CS = ',CS
         call sys%initialize(np=np,ns=ns,ne=ne,ncs=ncs,ng=ng,P=phse_mat,Ein=elem_mat,CS=CS,Bg=Bg,thermo_in=nasa_coef,diag=5)
         ! Initialize the chemical state
         call state%initialize(sys=sys,cond=fixed_PH,p=pressure)
         call param_read('Newton tolerance',state%tol_N)
         call param_read('Newton max iterations',state%iter_N_max)
         if (state%cond.eq.fixed_PH) then
            call param_read('T tolerance',state%tol_T)
            call param_read('T max iterations',state%iter_T_max)
            ! call param_read('Temperature initial guess',T_g)
            T_g=0.5_WP*(T_w+T_sat)
            ! T_g=T_sat
         end if
         ! Deallocate arrays
         deallocate(Bg)
      end block ceq_init


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
         timeSC%t=t0
         timeSC%tmax=time%tmax
         timeSC%dt=timeSC%dtmax
      end block initialize_timetracker


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
         ! Correct mask outside physical domain
         if (cfg%iproc.eq.1) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imin_-1
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%iproc.eq.cfg%npx) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imax_+1,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%jproc.eq.1) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmin_-1
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%jproc.eq.cfg%npy) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmax_+1,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%kproc.eq.1) then
            do k=cfg%kmino_,cfg%kmin_-1
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%kproc.eq.cfg%npz) then
            do k=cfg%kmax_+1,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
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
         use hypre_str_class, only: pcg_smg,gmres_smg,pcg_pfmg2,gmres_pfmg2
         use mathtools,       only: Pi
         use tpns_class,      only: clipped_neumann
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
         ! where (vf%VF.gt.0.0_WP) fs%U=beta*sqrt(alpha_g/(time%t))
         ! Apply boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      

      ! Initialize the interface location and velocity
      call get_interface()

      
      ! Create a one-sided scalar solver
      create_scalar: block
         use param,           only: param_read
         use tpscalar_class,  only: bcond,neumann
         use hypre_str_class, only: pcg_smg
         type(bcond), pointer :: mybc
         real(WP) :: N_init(ns),spDiff
         real(WP), dimension(:), allocatable :: mp
         integer :: i,j,k,n,isc,p
         integer  :: pos_open,pos_close
         allocate(mp(Lphase:Gphase)); mp=0.0_WP
         ! Read-in inputs
         call param_read('Air to water vapor mole ratio',air2wv_rat)
         call param_read('Nitrogen to oxygen mole ratio',N2O_rat)
         call param_read('Species diffusivity',spDiff)
         ! Create scalar solver
         call sc%initialize(cfg=cfg,nscalar=ns+2,name='tpscalar')
         sc%skip(get_sp_ind('H2O(L)'))=.true.
         ! Boundary conditinos
         call sc%add_bcond(name='xm',type=neumann,locator=xm_locator_sc,dir='-x')
         call sc%add_bcond(name='xp',type=neumann,locator=xp_locator   ,dir='+x')
         ! Assign scalar names and phases
         sc%SCname=[sp_names,'Tl','Tg']; iTl=ns+1; iTg=ns+2
         do isc=1,ns
            pos_open =index(sc%SCname(isc),'(')
            pos_close=index(sc%SCname(isc),')')
            if (pos_open.gt.0.and.pos_close.gt.pos_open) then
               sc%SCname(isc)(pos_open:pos_open)  ='_'
               sc%SCname(isc)(pos_close:pos_close)=' '
               sc%SCname(isc)=adjustl(trim(sc%SCname(isc)))
            end if
            sc%phase(isc)=sys%get_pind(isc)
         end do
         sc%SCname(get_sp_ind('H2O'))='H2O_g'
         sc%phase(iTl)=Lphase
         sc%phase(iTg)=Gphase
         ! Initialize the phasic density and VOF
         sc%Prho(Lphase)=fs%rho_l
         sc%Prho(Gphase)=fs%rho_g
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         ! Assign diffusivities
         sc%diff(:,:,:,1:ns)=spDiff
         sc%diff(:,:,:,iTl)=alpha_l
         sc%diff(:,:,:,iTg)=alpha_g
         ! Initialize the linear solver
         ss=hypre_str(cfg=cfg,name='Scalar',method=pcg_smg,nst=7)
         call param_read('Scalar iteration',ss%maxit)
         call param_read('Scalar tolerance',ss%rcvg)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
         ! Initialize mole numbers
         iWv=get_sp_ind('H2O')
         iWl=get_sp_ind('H2O(L)')
         ! iO2=get_sp_ind('O2')
         ! iN2=get_sp_ind('N2')
         print*,'iWv = ',iWv
         print*,'iWl = ',iWl
         ! print*,'iO2 = ',iO2
         ! print*,'iN2 = ',iN2
         N_init(iWl)=1.0_WP
         N_init(iWv)=1.0_WP
         ! N_init(iO2)=air2wv_rat*N_init(iWv)
         ! N_init(iN2)=N2O_rat*N_init(iO2)
         ! Get the phase mass
         mp=0.0_WP
         do isc=1,ns
            p=sc%phase(isc)
            mp(p)=mp(p)+N_init(isc)*MM(isc)
         end do
         ! Initialize interfacial scalar values
         allocate(sc_itf(sc%nscalar))
         do isc=1,ns
            p=sc%phase(isc)
            sc_itf(isc)=MM(isc)*N_init(isc)/mp(p)
         end do
         sc_itf(iTl:iTg)=T_sat
         print*,'Initial sc_itf = ',sc_itf
         ! Initialize scalar fields
         do isc=1,sc%nscalar
            p=sc%phase(isc)
            do k=sc%cfg%kmin_,sc%cfg%kmax_
               do j=sc%cfg%jmin_,sc%cfg%jmax_
                  do i=sc%cfg%imin_,sc%cfg%imax_
                     if (sc%PVF(i,j,k,p).gt.0.0_WP) then
                        sc%SC(i,j,k,isc)=sc_itf(isc)
                     else
                        sc%SC(i,j,k,isc)=0.0_WP
                     end if
                  end do
               end do
            end do
         end do
         ! Correct for the gas temperature
         do k=sc%cfg%kmin_,sc%cfg%kmax_
            do j=sc%cfg%jmin_,sc%cfg%jmax_
               do i=sc%cfg%imin_,sc%cfg%imax_
                  if (vf%VF(i,j,k).lt.1.0_WP) sc%SC(i,j,k,iTg)=T_w-(sc%cfg%xm(i)-cfg%x(cfg%imin))*(T_w-T_sat)/x_itf
               end do
            end do
         end do
         ! Apply boundary conditions at the interface
         call apply_itf_bcond()
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
         ! Get the phasic face apertures
         call sc%get_face_apt()
         ! Post process
         T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)
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
         ! Apply the interface jump conditions
         call interface_jump()
         call apply_itf_bcond()
         ! Get the volumetric pahse change mass flux
         call lg%get_mdot3p()
         ! Initialize the liquid and gas side mass fluxes
         call lg%init_mdot3pLG()
         ! Get the interface normal
         call lg%get_normal()
      end block create_lgpc


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         integer :: isc
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
         do isc=1,sc%nscalar
           call ens_out%add_scalar(trim(sc%SCname(isc)),sc%SC(:,:,:,isc))
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
         call mfile%add_column(x_itf,'x_itf')
         call mfile%add_column(u_itf,'u_itf')
         call mfile%add_column(x_itf_ext,'x_itf_ext')
         call mfile%add_column(u_itf_ext,'u_itf_ext')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%add_column(prhs_int,'prhs_int')
         ! Debug
         call fs%get_mfr()
         call cfg%integrate(lg%div_vel,mfr_err)
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
         call scfile%add_column(timeSC%n,'Timestep number')
         call scfile%add_column(timeSC%t,'Time')
         do isc=1,sc%nscalar
           call scfile%add_column(sc%SCmin(isc),trim(sc%SCname(isc))//'_min')
           call scfile%add_column(sc%SCmax(isc),trim(sc%SCname(isc))//'_max')
           call scfile%add_column(sc%SCint(isc),trim(sc%SCname(isc))//'_int')
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
      integer :: isc

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
         
         ! Update interface location
         call get_interface()
         print*,'i_itf = ',i_itf

         ! Update the phasic VOF and face apertures
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         call sc%get_face_apt()

         ! ================== SCALAR ================== !

         ! Backward Euler for diffusion and CN for advection
         advance_scalar: block
            use tpscalar_class,  only: bcond
            type(bcond), pointer :: mybc
            integer  :: i,j,k,n,p
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
               print*,'Old Tg = ',sc%SC(0:i_itf,1,1,iTg)

               ! Re-apply interface boundary conditions and zero out emptied cells
               call apply_itf_bcond()
               print*,'Old Tg after itf bcond = ',sc%SC(0:i_itf,1,1,iTg)

               ! Explicit calculation of dSC/dt
               call sc%get_dSCdt(dSCdt=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old)

               ! Assemble explicit residual
               resSC=resSC*timeSC%dt
               print*,'RHS Tg = ',resSC(0:i_itf,1,1,iTg)

               ! Form implicit residual
               call sc%solve_implicit(dt=timeSC%dt,resSC=resSC,U=fs%Uold,V=fs%Vold,W=fs%Wold,divU=lg%div_vel_old,w_adv=0.5_WP,w_dff=1.0_WP)
               print*,'sol Tg = ',resSC(0:i_itf,1,1,iTg)

               ! Apply the residuals
               do isc=1,sc%nscalar
                  p=sc%phase(isc)
                  where (sc%PVF(:,:,:,p).eq.1.0_WP) sc%SC(:,:,:,isc)=sc%SCold(:,:,:,isc)+resSC(:,:,:,isc)
               end do
               print*,'new Tg before bcond = ',sc%SC(0:i_itf,1,1,iTg)

               ! Apply boundary conditions
               call sc%apply_bcond(timeSC%t,timeSC%dt)
               call sc%get_bcond('xm',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  sc%SC(i,j,k,iTg)=2.0_WP*T_w-sc%SC(i+1,j,k,iTg)
               end do
               print*,'new Tg after bcond = ',sc%SC(0:i_itf,1,1,iTg)
            end do

         end block advance_scalar


         ! ================== PHASE CHANGE ================== !

         ! Apply the interface jump conditions
         call interface_jump()

         ! Get the interface location
         call get_interface()
         print*,'after jump:'
         print*,'i_itf = ',i_itf
         print*,'VOF(i_itf) = ',vf%VF(i_itf,1,1)

         ! Update the interface scalars
         do isc=1,sc%nscalar
            sc_itf(isc)=sc%SC(i_itf,1,1,isc)
         end do
         print*,'sc_itf = ',sc_itf

         ! Get the volumetric phase change mass flux
         call lg%get_mdot3p()

         ! Shift the phase change mass flux
         call lg%shift_mdot3p()
         
         ! Get the phase change induced divergence
         call lg%get_div()


         ! ================== VELOCITY ================== !

         advance_flow: block
            use tpns_class, only: static_contact,harmonic_visc
            
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
            ! Debug
            call fs%get_mfr()
            call cfg%integrate(lg%div_vel,mfr_err)
            mfr_err=abs(mfr_err-sum(fs%mfr))

         end block advance_flow


         ! ================== OUTPUT ================== !

         ! Output to ensight
         if (ens_evt%occurs()) then
            ! One-field temperature
            T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)
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