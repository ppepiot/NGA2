!> Two-dimensional laminar slot burner, following L. Selle, T. Poinsot and B. Ferret,
!> "Experimental and numerical study of the accuracy of flame-speed measurements for methane/air
!> combustion in a slot burner", Combustion and Flame 158 (2011) 146-154.
!>
!> The burner is a rectangular slot of width h = 10 mm and length L = 100 mm. The paper's DNS resolves the
!> central plane y = L/2, where the flow is two-dimensional, and that is what this driver reproduces: a slot
!> of width h at the bottom of the domain issuing premixed reactants into a slow air co-flow, with the
!> water-cooled slot lips imposed at T_lip, symmetry on the lateral boundaries and an outlet at the top
!> (paper Fig. 6a). A steady flame anchors on the lips and closes at a tip, and its length L_f gives the
!> mean consumption speed through s_L = U h / L_f (paper Eq. 1).
!>
!> Fig. 6 is the methane/air case; Fig. 11 is the same burner run on hydrogen/air at phi = 0.6 and 0.8.
!> Both are reachable here because the Igni73 skeletal mechanism carries a complete H2/O2 sub-mechanism.
!>
!> Chemistry is integrated explicitly inside the iterated-midpoint scheme of amrcomp; the NASA7 mixture
!> material carries the formation enthalpies so no energy source is needed.
module simulation
   use precision,         only: WP
   use string,            only: str_medium
   use amrgrid_class,     only: amrgrid
   use amrcomp_class,     only: amrcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use nasa7mix_class,    only: nasa7mix
   use amrchem_class,     only: amrchem
   use fcmech,            only: nS,nA
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid, time, solver, chemistry
   type(amrgrid), target :: amr
   type(timetracker) :: time
   type(amrcomp), target :: fs
   type(amrdata) :: dQdt
   type(amrdata) :: Yclose   !< Closure-species mass fraction, derived for output only
   type(nasa7mix), target :: gas
   type(amrchem) :: chem

   !> Output and monitoring
   type(event) :: viz_evt,regrid_evt
   type(amrviz) :: viz
   type(monitor) :: mfile,cflfile,consfile,chemfile,gridfile,flamefile

   !> Species bookkeeping: perm(n) is the fcmech index of material species n (closure species last)
   integer, dimension(nS) :: perm
   character(len=str_medium), dimension(nS) :: names   !< fcmech order
   integer :: ifuel,iO2,iOH,iH2O                       !< material indices (0 if absent or closure)
   integer :: kfuel                                    !< fcmech index of the fuel

   !> Thermochemical states (fcmech order)
   real(WP) :: P0,T_u,T_b,rho_u,q_comb
   real(WP), dimension(nS) :: Y_u,Y_b,Y_co
   logical :: have_burnt=.false.

   !> Burner geometry and operating point (paper Fig. 1 and Section 2)
   real(WP) :: Lx,Ly,Lz                     !< Domain: x transverse, y vertical, z collapsed
   real(WP) :: h_slot=0.01_WP               !< Slot width h [m]
   real(WP) :: w_lip=0.001_WP               !< Cooled lip thickness on each side of the slot [m]
   real(WP) :: T_lip=350.0_WP               !< Lip temperature (water cooled, paper Section 2) [K]
   real(WP) :: u_bulk=1.9_WP                !< Bulk velocity in the slot U [m/s]
   real(WP) :: u_co=0.1_WP                  !< Co-flow velocity [m/s]
   real(WP) :: delta_in=2.0e-4_WP           !< Smoothing thickness of the inlet profiles [m]
   logical  :: poiseuille=.false.           !< Parabolic rather than plug velocity profile across the slot
   real(WP) :: Lf_init=0.02_WP              !< Height of the initial flame cone [m]
   real(WP) :: t_ramp=0.0_WP                !< Inflow ramp time (0: full mass flux from t = 0) [s]
   real(WP) :: L_spg=0.0_WP                 !< Damping-layer thickness at the top outflow [m]
   real(WP) :: L_spg_side=0.0_WP            !< Damping-layer thickness at each lateral boundary [m]
   real(WP) :: spg_coeff=0.5_WP             !< Damping-layer strength: nu = spg_coeff * c * dx
   real(WP) :: delta_flame=5.0e-4_WP        !< Thickness of the initial flame front [m]

   !> Tagging thresholds
   real(WP) :: tag_grad=huge(1.0_WP),tag_hrr=huge(1.0_WP)

   !> LES (laminar by default: the paper's flames are laminar)
   logical :: use_sgs=.false.
   real(WP) :: Cs_vreman=0.17_WP

   !> Flame diagnostics (paper Eq. 1 and Section 5)
   real(WP) :: H_flame=0.0_WP               !< Flame tip height [m]
   real(WP) :: L_flame=0.0_WP               !< Flame length, TRIANG estimate [m]
   real(WP) :: sL_bar=0.0_WP                !< Mean consumption speed U h / L_f [m/s]
   real(WP) :: hrr_ratio=0.0_WP             !< Integrated heat release / expected at steady state

   !> Conservation diagnostics
   real(WP), dimension(:), allocatable, target :: elem_mass
   integer,  dimension(:,:), allocatable :: comp
   real(WP), dimension(:), allocatable :: atom_W
   character(len=2), dimension(:), allocatable :: atom_names

contains

   !> fcmech index of a species by name (0 if absent)
   integer function species_index(name) result(k)
      character(len=*), intent(in) :: name
      integer :: n
      k=0
      do n=1,nS
         if (trim(adjustl(names(n))).eq.trim(adjustl(name))) then; k=n; return; end if
      end do
   end function species_index

   !> Material index of a species by name (0 if absent or closure species)
   integer function material_index(name) result(n)
      character(len=*), intent(in) :: name
      integer :: k,m
      n=0
      k=species_index(name)
      if (k.eq.0) return
      do m=1,nS-1
         if (perm(m).eq.k) then; n=m; return; end if
      end do
   end function material_index

   !> Fresh mixture of Fuel with standard air (O2:N2 = 0.21:0.79 by mole) at equivalence ratio phi
   subroutine mixture_from_phi(fuel,phi,Y)
      use fcmech,   only: W_sp
      use messager, only: die
      character(len=*), intent(in) :: fuel
      real(WP), intent(in) :: phi
      real(WP), dimension(nS), intent(out) :: Y
      real(WP) :: nC,nH,nO,nO2st,FAst,YO2air,YF
      integer :: kf,kO2,kN2,a
      kf=species_index(fuel); kO2=species_index('O2'); kN2=species_index('N2')
      if (kf.eq.0.or.kO2.eq.0.or.kN2.eq.0) call die('[simulation mixture_from_phi] fuel, O2 and N2 must be in the mechanism')
      nC=0.0_WP; nH=0.0_WP; nO=0.0_WP
      do a=1,nA
         select case (trim(adjustl(atom_names(a))))
         case ('C'); nC=real(comp(a,kf),WP)
         case ('H'); nH=real(comp(a,kf),WP)
         case ('O'); nO=real(comp(a,kf),WP)
         end select
      end do
      nO2st=(2.0_WP*nC+0.5_WP*nH-nO)/2.0_WP
      if (nO2st.le.0.0_WP) call die('[simulation mixture_from_phi] fuel needs oxygen to burn')
      YO2air=0.21_WP*W_sp(kO2)/(0.21_WP*W_sp(kO2)+0.79_WP*W_sp(kN2))
      FAst=W_sp(kf)*YO2air/(nO2st*W_sp(kO2))
      YF=phi*FAst/(1.0_WP+phi*FAst)
      Y=0.0_WP
      Y(kf)=YF
      Y(kO2)=(1.0_WP-YF)*YO2air
      Y(kN2)=1.0_WP-Y(kf)-Y(kO2)
   end subroutine mixture_from_phi

   !> Conserved state of one cell from (T, fcmech-ordered Y, velocity) at the thermodynamic pressure P0
   subroutine fill_Q(Qcell,T,Yf,u,v,w)
      real(WP), dimension(:), intent(inout) :: Qcell
      real(WP), intent(in) :: T,u,v,w
      real(WP), dimension(nS), intent(in) :: Yf
      real(WP), dimension(nS) :: Ym
      real(WP) :: rho,e
      integer :: n
      do n=1,nS; Ym(n)=Yf(perm(n)); end do
      rho=gas%get_rho_from_p_T(p=P0,T=T,y=Ym)
      e  =gas%get_e_from_p_T(p=P0,T=T,y=Ym)
      Qcell(1)=rho
      Qcell(2)=rho*u
      Qcell(3)=rho*v
      Qcell(4)=rho*w
      Qcell(5)=rho*e
      Qcell(fs%Y_lo:fs%Y_hi)=rho*Ym(1:nS-1)
   end subroutine fill_Q

   !> Smooth step from 0 (left) to 1 (right) over a thickness delta
   real(WP) function step(x,x0,delta)
      real(WP), intent(in) :: x,x0,delta
      step=0.5_WP*(1.0_WP+tanh(4.0_WP*(x-x0)/delta))
   end function step

   !> Fraction of the cooled lip at transverse position x: 1 over the rim, 0 in the slot and the co-flow
   real(WP) function lip_frac(x)
      real(WP), intent(in) :: x
      real(WP) :: ax
      ax=abs(x)
      lip_frac=step(ax-0.5_WP*h_slot,0.0_WP,delta_in)*step(0.5_WP*h_slot+w_lip-ax,0.0_WP,delta_in)
   end function lip_frac

   !> Inflow ramp factor: 0 at t = 0 rising to 1 over t_ramp, or 1 throughout when no ramp is requested
   real(WP) function ramp(t)
      real(WP), intent(in) :: t
      ramp=1.0_WP
      if (t_ramp.gt.0.0_WP) ramp=1.0_WP-exp(-3.0_WP*t/t_ramp)
   end function ramp

   !> Inlet state at transverse position x (paper Fig. 1b and Section 2): premixed reactants over the slot
   !> |x| < h/2, a no-slip band at the water-cooled lip temperature over the lip, and air co-flow beyond.
   subroutine inlet_state(x,T,Yf,v)
      real(WP), intent(in) :: x
      real(WP), intent(out) :: T,v
      real(WP), dimension(nS), intent(out) :: Yf
      real(WP) :: s,sl,ax,prof
      ax=abs(x)
      ! Slot indicator: 1 inside the slot, 0 outside
      s=step(0.5_WP*h_slot-ax,0.0_WP,delta_in)
      ! Lip indicator: 1 over the cooled rim, 0 in the slot and in the co-flow
      sl=lip_frac(x)
      ! Transverse velocity profile in the slot
      prof=1.0_WP
      if (poiseuille) prof=max(0.0_WP,1.5_WP*(1.0_WP-(2.0_WP*x/h_slot)**2))
      T=T_u+(T_lip-T_u)*sl
      Yf=Y_u*s+Y_co*(1.0_WP-s)
      v=(u_bulk*prof*s+u_co*(1.0_WP-s))*(1.0_WP-sl)
   end subroutine inlet_state

   !> Initial condition: a conical flame anchored on the slot lips and closing at height Lf_init. Fresh gas
   !> inside the cone, burnt gas above/around it, air in the co-flow. Starting from the converged shape
   !> rather than from a cold jet saves the (long) ignition transient.
   subroutine burner_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP), dimension(nS) :: Yf,Ym
      real(WP) :: x,y,sflame,sjet,slip,T,rho,v
      integer :: i,j,k,n
      call amrex_mfiter_build(mfi,ba,dm,tiling=.true.)
      do while (mfi%next())
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            x=solver%amr%xlo+(real(i,WP)+0.5_WP)*solver%amr%dx(lvl)
            y=solver%amr%ylo+(real(j,WP)+0.5_WP)*solver%amr%dy(lvl)
            ! Signed distance to the conical flame sheet, whose two branches run from the lips
            ! (|x| = h/2, y = 0) to the tip (x = 0, y = Lf_init): negative in the fresh core, positive in
            ! the burnt gas. A distance is needed rather than |x| - x_cone(y): above the tip x_cone
            ! collapses to zero, so an |x|-based indicator becomes independent of y and leaves a ribbon of
            ! partly unburnt mixture streaming from the tip to the outlet (half the fresh fuel on the axis
            ! itself), which is an artefact of the initial condition and not a flame.
            sflame=step((abs(x)*Lf_init+y*0.5_WP*h_slot-0.5_WP*h_slot*Lf_init) &
            &           /sqrt(Lf_init**2+(0.5_WP*h_slot)**2),0.0_WP,delta_flame)
            ! Column above the slot (burnt plume) versus the co-flow
            sjet=step(0.5_WP*h_slot+w_lip-abs(x),0.0_WP,delta_in)
            T =T_u+(T_b-T_u)*sflame*sjet
            Yf=Y_u+(Y_b-Y_u)*sflame
            Yf=Yf*sjet+Y_co*(1.0_WP-sjet)
            ! Near the inlet, impose the same cooled no-slip lip the boundary condition will apply, decaying
            ! with height. Without this the initial field carries flow and fresh-gas temperature across the
            ! lip while the boundary sets u = 0 and T = T_lip there, and that mismatch fires an acoustic
            ! pulse at t = 0 exactly where the oscillations are seen.
            slip=lip_frac(x)*exp(-(y/max(delta_in,2.0_WP*solver%amr%dx(lvl)))**2)
            T=T+(T_lip-T)*slip
            do n=1,nS; Ym(n)=Yf(perm(n)); end do
            rho=gas%get_rho_from_p_T(p=P0,T=T,y=Ym)
            ! Keep the mass flux of the slot stream: the burnt gas accelerates by the expansion ratio
            v=(u_bulk*rho_u/rho)*sjet+u_co*(1.0_WP-sjet)
            v=v*(1.0_WP-slip)*ramp(0.0_WP)
            call fill_Q(pQ(i,j,k,:),T,Yf,0.0_WP,v,0.0_WP)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine burner_init

   !> Boundary conditions: inflow at the bottom face (face 3). The lateral faces are symmetry and the top
   !> is an outflow, both handled by the reflect/foextrap BC types set in simulation_init.
   subroutine burner_bc(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP), dimension(nS) :: Yf
      real(WP) :: x,T,v
      integer :: i,j,k
      if (face.ne.3) return
      do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
         ! Cell-centred x for Q and for the tangential components; the V face sits at the same x
         x=solver%amr%xlo+(real(i,WP)+0.5_WP)*solver%amr%dx(lvl)
         call inlet_state(x,T,Yf,v)
         v=v*ramp(time)
         select case (comp)
         case ('U','W'); p(i,j,k,1)=0.0_WP
         case ('V');     p(i,j,k,1)=v
         case ('Q');     call fill_Q(p(i,j,k,:),T,Yf,0.0_WP,v,0.0_WP)
         end select
      end do; end do; end do
   end subroutine burner_bc

   !> Tagging: relative density jump and heat release rate
   subroutine burner_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pH
      real(WP) :: grad
      integer :: i,j,k
      tags=tags_ptr
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         pH=>chem%hrr%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Relative density jump per cell (density is available when the initial levels are built, T is not)
            grad=max(abs(pQ(i+1,j,k,1)-pQ(i-1,j,k,1)),abs(pQ(i,j+1,k,1)-pQ(i,j-1,k,1)))*0.5_WP &
            &   /max(pQ(i,j,k,1),solver%rho_floor)
            if (grad.gt.tag_grad) tagarr(i,j,k,1)=SETtag
            if (pH(i,j,k,1).gt.tag_hrr) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine burner_tagger

   !> Flame length and mean consumption speed, following the paper's IMAGE/TRIANG post-processing:
   !> the reaction layer is the region where the heat release rate exceeds half its maximum (which
   !> excludes the quenched flame base), its highest point gives the flame tip H_f, and the flame length
   !> follows from the triangle fitted between the lips and the tip. Eq. (1) then gives s_L = U h / L_f.
   subroutine get_flame()
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
      use parallel,         only: MPI_REAL_WP
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pH
      real(WP) :: y,thr
      integer :: lvl,i,j,k,ierr
      H_flame=0.0_WP
      thr=0.5_WP*chem%hrr_max
      if (thr.gt.0.0_WP) then
         do lvl=0,amr%clvl()
            call amrex_mfiter_build(mfi,chem%hrr%mf(lvl),tiling=.true.)
            do while (mfi%next())
               bx=mfi%tilebox()
               pH=>chem%hrr%mf(lvl)%dataptr(mfi)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  if (pH(i,j,k,1).gt.thr) then
                     y=amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl)
                     H_flame=max(H_flame,y)
                  end if
               end do; end do; end do
            end do
            call amrex_mfiter_destroy(mfi)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,H_flame,1,MPI_REAL_WP,MPI_MAX,amr%comm,ierr)
      end if
      ! TRIANG estimate: both sides of the triangle fitted between the lips and the tip
      L_flame=2.0_WP*sqrt((0.5_WP*h_slot)**2+H_flame**2)
      sL_bar=0.0_WP
      if (L_flame.gt.0.0_WP) sL_bar=u_bulk*h_slot/L_flame
      ! At steady state the integrated heat release must balance the fuel fed through the slot
      hrr_ratio=0.0_WP
      if (q_comb.gt.0.0_WP) hrr_ratio=chem%hrr_int/(rho_u*u_bulk*h_slot*Lz*q_comb)
   end subroutine get_flame

   !> Closure-species mass fraction for the output: the solver transports only nS-1 species, the last one
   !> is recovered as 1 - sum of the others. Recomputed before every plotfile write, so it never goes stale.
   subroutine update_closure()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pY,pYc
      integer :: lvl,i,j,k
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pY =>fs%Y%mf(lvl)%dataptr(mfi)
            pYc=>Yclose%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pYc(i,j,k,1)=max(0.0_WP,1.0_WP-sum(pY(i,j,k,1:nS-1)))
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine update_closure

   !> Element masses and total mass (conservation check)
   subroutine get_conservation()
      real(WP), dimension(nS) :: Mk
      real(WP) :: Mclose
      integer :: k,a
      Mclose=fs%Qint(1)
      do k=1,nS
         if (chem%iQ(k).gt.0) then
            Mk(k)=fs%Qint(chem%iQ(k)); Mclose=Mclose-Mk(k)
         end if
      end do
      Mk(chem%kclose)=Mclose
      do a=1,nA
         elem_mass(a)=0.0_WP
         do k=1,nS
            elem_mass(a)=elem_mass(a)+real(comp(a,k),WP)*atom_W(a)/gas%W(material_of(k))*Mk(k)
         end do
      end do
   contains
      integer function material_of(kf)
         integer, intent(in) :: kf
         integer :: n
         material_of=0
         do n=1,nS
            if (perm(n).eq.kf) then; material_of=n; return; end if
         end do
      end function material_of
   end subroutine get_conservation

   !> Damping layer at the artificial far boundaries (top outflow and, optionally, the lateral symmetry
   !> planes). Every non-periodic face of this set-up reflects: ext_dir at the inlet pins the pressure,
   !> reflect_* on the sides are rigid walls and foextrap at the top is zero-gradient, not characteristic.
   !> Nothing here makes them non-reflecting -- that needs NSCBC -- but raising the viscosities over a layer
   !> in front of them absorbs part of what would otherwise bounce back. The bulk viscosity is the term that
   !> acts on the dilatational (acoustic) field, so it is raised alongside the shear viscosity.
   !>
   !> nu is set from the local sound speed and cell size, NOT from 1/dt: a dt-based coefficient pins the
   !> viscous CFL at a fixed value and drives the time step geometrically to zero whenever Max CFL is set
   !> below it.
   subroutine apply_sponge()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pC,pT,pVisc,pBeta,pDiff,pDiffY
      real(WP) :: x,y,f,nu,rho
      integer :: lvl,i,j,k
      if (L_spg.le.0.0_WP.and.L_spg_side.le.0.0_WP) return
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%growntilebox(fs%nover)
            pQ    =>fs%Q%mf(lvl)%dataptr(mfi)
            pC    =>fs%C%mf(lvl)%dataptr(mfi)
            pT    =>fs%T%mf(lvl)%dataptr(mfi)
            pVisc =>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta =>fs%beta%mf(lvl)%dataptr(mfi)
            pDiff =>fs%diff%mf(lvl)%dataptr(mfi)
            pDiffY=>fs%diffY%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=amr%xlo+(real(i,WP)+0.5_WP)*amr%dx(lvl)
               y=amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl)
               ! Quadratic ramp, zero at the start of the layer so the layer itself does not reflect
               f=0.0_WP
               if (L_spg.gt.0.0_WP)      f=max(f,max(0.0_WP,(y-(amr%yhi-L_spg))/L_spg)**2)
               if (L_spg_side.gt.0.0_WP) f=max(f,max(0.0_WP,(abs(x)-(0.5_WP*(amr%xhi-amr%xlo)-L_spg_side))/L_spg_side)**2)
               if (f.le.0.0_WP) cycle
               rho=max(pQ(i,j,k,1),fs%rho_floor)
               nu=spg_coeff*f*pC(i,j,k,1)*amr%dx(lvl)
               pVisc (i,j,k,1)=pVisc (i,j,k,1)+rho*nu
               pBeta (i,j,k,1)=pBeta (i,j,k,1)+rho*nu
               pDiff (i,j,k,1)=pDiff (i,j,k,1)+rho*nu*chem%tr%get_cp(pT(i,j,k,1),Y_u)
               pDiffY(i,j,k,1)=pDiffY(i,j,k,1)+rho*nu
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine apply_sponge

   !> Transport coefficients (and SGS if requested)
   subroutine update_transport()
      if (use_sgs) then
         call chem%update_properties(fs=fs,dt=time%dt,Cs=Cs_vreman)
      else
         call chem%update_properties(fs=fs,dt=time%dt)
      end if
      call apply_sponge()
   end subroutine update_transport

   !> Initialization
   subroutine simulation_init
      use param,    only: param_read,param_exists
      use messager, only: die
      use string,   only: lowercase
      use fcmech,   only: W_sp,T_mid,thermo_coeffs,Rcst,P_ref,fcmech_get_speciesnames,fcmech_get_composition,fcmech_get_atommasses,fcmech_get_atomnames
      implicit none

      ! Species bookkeeping and the closure species (permuted last for amrcomp)
      species_setup: block
         character(len=str_medium) :: closure
         integer :: k,n,kc
         call fcmech_get_speciesnames(names)
         allocate(comp(nA,nS),atom_W(nA),atom_names(nA),elem_mass(nA))
         call fcmech_get_atomnames(atom_names)
         call fcmech_get_atommasses(atom_W)
         call fcmech_get_composition(comp)
         call param_read('Closure species',closure,default='N2')
         kc=species_index(closure)
         if (kc.eq.0) call die('[simulation init] Closure species not found in the mechanism')
         n=0
         do k=1,nS
            if (k.ne.kc) then; n=n+1; perm(n)=k; end if
         end do
         perm(nS)=kc
      end block species_setup

      ! Operating point (paper Section 2: atmospheric pressure, fresh gas at 300 K)
      thermo_setup: block
         character(len=str_medium) :: fuel
         real(WP) :: phi,val
         integer :: k
         call param_read('Fuel',fuel)
         call param_read('Equivalence ratio',phi)
         call param_read('Temperature',T_u,default=300.0_WP)
         call param_read('Pressure',P0,default=101325.0_WP)
         kfuel=species_index(fuel)
         call mixture_from_phi(trim(fuel),phi,Y_u)
         ! Co-flow: air at the fresh-gas temperature
         call mixture_from_phi(trim(fuel),0.0_WP,Y_co)
         ! Burnt state: read from the input (Cantera reference) or fall back to the fresh state
         call param_read('Burnt temperature',T_b,default=T_u)
         Y_b=Y_u
         have_burnt=param_exists('Burnt temperature')
         if (have_burnt) then
            Y_b=0.0_WP
            do k=1,nS
               if (param_exists('Burnt Y '//trim(adjustl(names(k))))) then
                  call param_read('Burnt Y '//trim(adjustl(names(k))),val); Y_b(k)=val
               end if
            end do
            if (sum(Y_b).le.0.0_WP) call die('[simulation init] Burnt temperature given without any Burnt Y')
            Y_b=Y_b/sum(Y_b)
         end if
      end block thermo_setup

      ! NASA7 material built from the mechanism's thermo table, closure species last
      create_material: block
         real(WP), dimension(nS,14) :: coeffs
         character(len=str_medium), dimension(nS) :: mnames
         integer :: n
         do n=1,nS
            coeffs(n,:)=thermo_coeffs(perm(n),:)
            mnames(n)=names(perm(n))
         end do
         call gas%initialize(W=W_sp(perm),T_mid=T_mid(perm),coeffs=coeffs,species_names=mnames, &
         &                   name='Igni73 mixture',Ru=Rcst,p_ref=P_ref)
      end block create_material

      ! AMR grid: x transverse (slot centred on x = 0), y vertical, z collapsed to a 2D plane
      create_amrgrid: block
         real(WP) :: dxf
         amr%name='amrcomp_slotburner'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz,default=2)
         call param_read('Max level',amr%maxlvl)
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         call param_read('Lz',Lz,default=Lx)
         dxf=Lx/real(amr%nx*2**amr%maxlvl,WP)
         amr%xlo=-0.5_WP*Lx; amr%xhi=+0.5_WP*Lx
         amr%ylo=0.0_WP;     amr%yhi=Ly
         amr%zlo=-0.5_WP*Lz; amr%zhi=+0.5_WP*Lz
         ! Collapsed spanwise direction: 2 cells of the finest size with unit refinement ratio
         if (amr%nz.le.2) then
            amr%rrefz=[1]; amr%zlo=-0.5_WP*real(amr%nz,WP)*dxf; amr%zhi=-amr%zlo; Lz=2.0_WP*amr%zhi
         end if
         if (amr%nz.eq.1.and.amr%maxlvl.gt.0) call die('[simulation init] the collapsed direction needs 2 cells when Max level > 0')
         amr%xper=.false.; amr%yper=.false.; amr%zper=.true.
         call amr%initialize()
      end block create_amrgrid

      ! Time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot,name='slotburner')
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
      end block initialize_timetracker

      ! Compressible solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap,amrex_bc_reflect_even,amrex_bc_reflect_odd
         use amrdata_class,    only: interp_face_lin,interp_blin,interp_ublin,interp_const
         character(len=str_medium) :: qinterp
         call param_read('Use projection',fs%use_projection,default=.false.)
         fs%mat=>gas
         call param_read('Q interpolation',qinterp,default='blin')
         select case (lowercase(trim(qinterp)))
         case ('blin');  fs%interp_Q=interp_blin
         case ('ublin'); fs%interp_Q=interp_ublin
         case ('const'); fs%interp_Q=interp_const
         case default; call die('[simulation init] Unknown Q interpolation (blin, ublin, const)')
         end select
         call fs%initialize(amr=amr,name='amrcomp slotburner')
         if (amr%nz.le.2) fs%interp_vel=interp_face_lin
         fs%psolver%max_iter=20
         fs%psolver%tol_rel=1.0e-6_WP
         fs%user_init=>burner_init
         fs%user_tagging=>burner_tagger
         fs%user_bc=>burner_bc
         ! Lateral boundaries: symmetry (paper Fig. 6a). x-momentum and the x face velocity are odd.
         fs%Q%lo_bc(1,:)=amrex_bc_reflect_even; fs%Q%hi_bc(1,:)=amrex_bc_reflect_even
         fs%Q%lo_bc(1,2)=amrex_bc_reflect_odd;  fs%Q%hi_bc(1,2)=amrex_bc_reflect_odd
         fs%U%lo_bc(1,1)=amrex_bc_reflect_odd;  fs%U%hi_bc(1,1)=amrex_bc_reflect_odd
         fs%V%lo_bc(1,1)=amrex_bc_reflect_even; fs%V%hi_bc(1,1)=amrex_bc_reflect_even
         fs%W%lo_bc(1,1)=amrex_bc_reflect_even; fs%W%hi_bc(1,1)=amrex_bc_reflect_even
         ! Bottom: inflow through user_bc. Top: outflow.
         fs%Q%lo_bc(2,:)=amrex_bc_ext_dir;  fs%Q%hi_bc(2,:)=amrex_bc_foextrap
         fs%U%lo_bc(2,1)=amrex_bc_ext_dir;  fs%U%hi_bc(2,1)=amrex_bc_foextrap
         fs%V%lo_bc(2,1)=amrex_bc_ext_dir;  fs%V%hi_bc(2,1)=amrex_bc_foextrap
         fs%W%lo_bc(2,1)=amrex_bc_ext_dir;  fs%W%hi_bc(2,1)=amrex_bc_foextrap
         ! Reference fresh state and heat of combustion per unit mass of mixture
         reference_state: block
            real(WP), dimension(nS) :: Ym,Ymb
            integer :: n
            do n=1,nS; Ym(n)=Y_u(perm(n)); Ymb(n)=Y_b(perm(n)); end do
            rho_u=gas%get_rho_from_p_T(p=P0,T=T_u,y=Ym)
            q_comb=gas%get_h_from_p_T(p=P0,T=T_u,y=Ym)-gas%get_h_from_p_T(p=P0,T=T_u,y=Ymb)
         end block reference_state
      end block create_solver

      ! Chemistry and transport
      create_chemistry: block
         character(len=str_medium) :: model
         real(WP) :: val
         call chem%initialize(amr=amr,fs=fs,name='finite-rate chemistry')
         call param_read('Chemical CFL',chem%CFLchem_max,default=1.0_WP)
         call param_read('Chemistry T min',chem%T_min,default=0.0_WP)
         call param_read('Transport model',model,default='mixavg')
         select case (lowercase(trim(model)))
         case ('mixavg');     call chem%tr%set_mixavg()
         case ('powerlaw');   call chem%tr%set_powerlaw()
         case ('sutherland'); call chem%tr%set_sutherland()
         case default; call die('[simulation init] Unknown Transport model (mixavg, powerlaw, sutherland)')
         end select
         call param_read('Prandtl number',val,default=0.68_WP); chem%tr%Pr=val
         call param_read('Lewis number',val,default=1.0_WP);    chem%tr%Le=val
         if (param_exists('Schmidt number')) then
            call param_read('Schmidt number',val); chem%tr%Sc=val
         end if
         call param_read('Use SGS',use_sgs,default=.false.)
         call param_read('Vreman constant',Cs_vreman,default=0.17_WP)
      end block create_chemistry

      ! Burner geometry and operating point
      burner_setup: block
         call param_read('Slot width',h_slot)
         call param_read('Lip thickness',w_lip,default=0.1_WP*h_slot)
         call param_read('Lip temperature',T_lip,default=350.0_WP)
         call param_read('Bulk velocity',u_bulk)
         call param_read('Coflow velocity',u_co,default=0.05_WP*u_bulk)
         ! The inlet profile has to be resolved on the BASE grid, not just on the finest level: a sub-cell
         ! step in velocity and temperature at the lip is a strong acoustic source right at the inlet.
         call param_read('Inlet smoothing',delta_in,default=max(0.02_WP*h_slot,3.0_WP*amr%dx(0)))
         call param_read('Inflow ramp time',t_ramp,default=0.0_WP)
         call param_read('Sponge length',L_spg,default=0.0_WP)
         call param_read('Side sponge length',L_spg_side,default=0.0_WP)
         call param_read('Sponge strength',spg_coeff,default=0.5_WP)
         call param_read('Poiseuille inlet',poiseuille,default=.false.)
         call param_read('Initial flame height',Lf_init,default=2.0_WP*h_slot)
         call param_read('Initial flame thickness',delta_flame,default=0.05_WP*h_slot)
         call param_read('Tag grad',tag_grad,default=huge(1.0_WP))
         call param_read('Tag HRR',tag_hrr,default=huge(1.0_WP))
      end block burner_setup

      ! Workspace, initial grid and fields
      init_fields: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call Yclose%initialize(amr,name='Yclose',ncomp=1,ng=0,interp=interp_none); call Yclose%register()
         ifuel=material_index(trim(adjustl(names(kfuel))))
         iO2=material_index('O2'); iOH=material_index('OH'); iH2O=material_index('H2O')
         call amr%init_from_scratch(time=time%t)
         call fs%get_primitive(Q=fs%Q); call fs%get_face_velocity()
         call update_transport()
         ! Prime the chemical time-scale limiter so the first step is already chemistry-limited
         call chem%add_source(fs=fs,dQdt=dQdt,dt=time%dt)
         call chem%get_stats(fs=fs)
      end block init_fields

      ! Events
      create_events: block
         real(WP) :: out_period
         integer :: nregrid,out_nsteps
         ! Plotfiles are written on elapsed simulated time, on elapsed steps, or on whichever comes first.
         ! With explicit chemistry the step is a few nanoseconds, so a purely time-based period can mean
         ! tens of thousands of steps before the first file appears: 'Output nsteps' is the way to see
         ! something early in a long run.
         call param_read('Output period',out_period,default=0.0_WP)
         call param_read('Output nsteps',out_nsteps,default=0)
         if (out_period.le.0.0_WP.and.out_nsteps.le.0) &
         &   call die('[simulation init] give Output period (in seconds) or Output nsteps')
         viz_evt=event(time=time,name='Output')
         viz_evt%tper=out_period
         viz_evt%nper=out_nsteps
         call param_read('Regrid nsteps',nregrid,default=20)
         regrid_evt=event(time=time,name='Regrid')
         regrid_evt%nper=nregrid
      end block create_events

      ! Visualization
      create_viz: block
         integer :: n
         call viz%initialize(amr=amr,name='slotburner',use_hdf5=.false.)
         ! Thermodynamic state
         call viz%add_scalar(fs%Q,1,'RHO')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%T,1,'T')
         call viz%add_scalar(fs%C,1,'SoS')
         ! Velocity (cell-centred)
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         ! Every species: the solver transports nS-1 of them, the closure species is the remainder
         do n=1,nS-1
            call viz%add_scalar(fs%Y,n,'Y_'//trim(adjustl(names(perm(n)))))
         end do
         call viz%add_scalar(Yclose,1,'Y_'//trim(adjustl(names(perm(nS)))))
         ! Chemistry and transport
         call viz%add_scalar(chem%hrr,1,'HRR')
         call viz%add_scalar(chem%tauchem,1,'tau_chem')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%diff,1,'lambda')
         call viz%add_scalar(fs%diffY,1,'rhoD')
         call update_closure()
         call viz%write(time=time%t)
      end block create_viz

      ! Monitors
      create_monitors: block
         integer :: a
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call chem%get_stats(fs=fs)
         call get_flame()
         call get_conservation()
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Pmin,'Pmin')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%Tmin,'Tmin')
         call mfile%add_column(fs%Tmax,'Tmax')
         if (ifuel.gt.0) call mfile%add_column(fs%Ymax(ifuel),'Yfuel max')
         if (iOH.gt.0)   call mfile%add_column(fs%Ymax(iOH),'YOH max')
         call mfile%write()
         flamefile=monitor(amRoot=amr%amRoot,name='flame')
         call flamefile%add_column(time%n,'Timestep number')
         call flamefile%add_column(time%t,'Time')
         call flamefile%add_column(H_flame,'Flame height')
         call flamefile%add_column(L_flame,'Flame length')
         call flamefile%add_column(sL_bar,'sL bar')
         call flamefile%add_column(chem%hrr_int,'HRR integral')
         call flamefile%add_column(hrr_ratio,'HRR/expected')
         call flamefile%add_column(fs%Tmax,'Tmax')
         call flamefile%write()
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(chem%CFLchem,'CFLchem')
         call cflfile%add_column(chem%tau_min,'tau_chem min')
         call cflfile%write()
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Mass')
         call consfile%add_column(fs%Qint(5),'Internal energy')
         do a=1,nA
            call consfile%add_column(elem_mass(a),'Mass '//trim(atom_names(a)))
         end do
         call consfile%write()
         chemfile=monitor(amRoot=amr%amRoot,name='chemistry')
         call chemfile%add_column(time%n,'Timestep number')
         call chemfile%add_column(time%t,'Time')
         call chemfile%add_column(chem%ncell_react,'Reacting cells')
         call chemfile%add_column(chem%hrr_int,'HRR integral')
         call chemfile%add_column(chem%hrr_max,'HRR max')
         call chemfile%add_column(chem%CFLchem,'CFLchem')
         call chemfile%add_column(chem%tau_min,'tau_chem min')
         call chemfile%add_column(chem%wtime,'Chemistry time')
         call chemfile%write()
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%write()
      end block create_monitors

   end subroutine simulation_init

   !> Time advancement
   subroutine simulation_run
      implicit none
      do while (.not.time%done())
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         if (chem%CFLchem_max.gt.0.0_WP) time%cfl=max(time%cfl,time%cflmax*chem%CFLchem/chem%CFLchem_max)
         call time%adjust_dt()
         call time%increment()
         call fs%Qold%copy(src=fs%Q)
         call fs%Uold%copy(src=fs%U)
         call fs%Vold%copy(src=fs%V)
         call fs%Wold%copy(src=fs%W)
         do while (time%it.le.time%itmax)
            call fs%Q%lincomb(a=0.5_WP,src1=fs%Qold,b=0.5_WP,src2=fs%Q)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)
            call fs%get_primitive(Q=fs%Q)
            call fs%get_dQdt(dQdt=dQdt)
            call chem%add_source(fs=fs,dQdt=dQdt,dt=time%dt)
            call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
            call fs%Q%average_down(); call fs%Q%fill(time%t)
            call fs%get_face_velocity(); call fs%average_down_velocity()
            call fs%get_primitive(Q=fs%Q)
            call fs%add_pressure(scale=time%dt,phi=fs%P)
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
            if (fs%use_projection) then
               call fs%get_div(); call fs%div%mult(val=1.0_WP/time%dt)
               call fs%prepare_psolver(dt=time%dt)
               call fs%psolver%solve(rhs=fs%div)
               call fs%add_pressure(scale=time%dt)
               call fs%Q%average_down(); call fs%Q%fill(time=time%t)
               call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
            end if
            time%it=time%it+1
         end do
         call fs%get_primitive(Q=fs%Q)
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if
         call update_transport()
         if (viz_evt%occurs()) then
            call update_closure()
            call viz%write(time%t)
         end if
         call fs%get_info()
         call chem%get_stats(fs=fs)
         call get_flame()
         call get_conservation()
         call mfile%write()
         call flamefile%write()
         call cflfile%write()
         call consfile%write()
         call chemfile%write()
      end do
   end subroutine simulation_run

   !> Finalization
   subroutine simulation_final
      implicit none
      call time%finalize()
      call amr%finalize()
      call regrid_evt%finalize()
      call fs%finalize()
      call dQdt%finalize()
      call Yclose%finalize()
      call chem%finalize()
      call gas%finalize()
      call viz%finalize()
      call viz_evt%finalize()
      call mfile%finalize()
      call flamefile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call chemfile%finalize()
      call gridfile%finalize()
   end subroutine simulation_final

end module simulation
