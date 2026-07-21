!==================================================================================================
!  PROGRAM thermo0D
!
!  0D constant-volume reactor. Supports:
!    - pure liquid water (NASG), Liquid volume fraction = 1
!    - pure nitrogen (ideal gas), Liquid volume fraction = 0
!    - mixture (e.g. half liquid water / half pure N2), 0 < alpha < 1
!
!  Specific internal energy of each present phase rises at Heating rate [J/(kg.s)].
!  Total volume is fixed. For mixtures, mechanical (p) or mechanical+thermal (pT)
!  relaxation equalizes the phases after each heat step.
!
!  Build:  make
!  Run:    mpirun -np 1 ./thermo0D.dp.gnu.opt.mpi.exe -i input
!  Plot:   python3 plot_thermo_interactive.py data/thermo.csv
!==================================================================================================
program thermo0D
   use precision,           only: WP
   use string,              only: str_medium,str_long,lowercase
   use param
   use parallel
   use messager,            only: messager_init,messager_final,die,warn,log
   use nasg_class,          only: nasg
   use ideal_gas_class,     only: ideal_gas
   use relax_ig_nasg_class, only: relax_ig_nasg,Prelax,PTrelax
   use thermorelax_class,   only: RELAX_OK
   implicit none

   ! Water IAPWS critical point (reference thresholds for flags; NASG itself has no critical point)
   real(WP), parameter :: Tcrit_H2O = 647.096_WP      !< K
   real(WP), parameter :: Pcrit_H2O = 22.064e6_WP     !< Pa
   real(WP), parameter :: VF_PURE = 1.0e-12_WP

   type(nasg),      target :: water
   type(ideal_gas), target :: nitrogen
   type(relax_ig_nasg)     :: relaxer
   real(WP), dimension(1)  :: y
   real(WP), dimension(7)  :: Q   !< [mL,mG,EL,EG,mx,my,mz] partial densities / energy densities

   ! Liquid NASG parameters
   real(WP) :: gammaL,pinfL,bL,cvL,qL,qpL
   ! Gas (N2) ideal-gas parameters
   real(WP) :: gammaG,cvG,qG,qpG

   ! Reactor / ICs
   real(WP) :: Vtot,alpha
   real(WP) :: P0,T0,PL0,TL0,PG0,TG0
   real(WP) :: rhoL,rhoG,eL,eG,massL,massG,mass
   real(WP) :: PL,PG,P,TL,TG,Temp
   real(WP) :: hL,hG,sL,sG,gL,gG,cL,cG,cvL_loc,cvG_loc
   real(WP) :: e_mix,rho_mix,vL,vG,brho
   real(WP) :: heating_rate,Q_added
   real(WP) :: time,t_end,dt

   ! Critical / packing flags
   real(WP) :: Tcrit,Pcrit,brhomax_warn
   integer  :: flag_Tcrit,flag_Pcrit,flag_critical,flag_packing
   logical  :: hit_Tcrit,hit_Pcrit,hit_critical,hit_packing
   logical  :: stop_at_critical,announced_T,announced_P,announced_pack
   logical  :: has_liq,has_gas,do_relax

   ! Relaxation
   character(len=str_medium) :: relax_type
   integer :: relax_ierr

   ! I/O
   character(len=str_medium) :: outfile
   character(len=str_long) :: msg
   integer :: iunit,ierr,nstep,nout
   logical :: am_root

   !----------------------------------------------------------------------------------------------
   call parallel_init
   call messager_init
   call param_init
   am_root = amRoot
   y = 1.0_WP
   Q = 0.0_WP

   ! ---- Liquid NASG (defaults: Le Métayer & Saurel 2016, Table V) --------------------------------
   gammaL = 1.19_WP
   pinfL  = 7.028e8_WP
   bL     = 6.61e-4_WP
   cvL    = 3610.0_WP
   qL     =-1.177788e6_WP
   qpL    = 0.0_WP
   call param_read('Gamma', gammaL, default=gammaL)
   call param_read('Pinf',  pinfL,  default=pinfL)
   call param_read('b',     bL,     default=bL)
   call param_read('Cv',    cvL,    default=cvL)
   call param_read('q',     qL,     default=qL)
   call param_read('qp',    qpL,    default=qpL)
   call water%initialize(gamma=gammaL,pinf=pinfL,b=bL,cv=cvL,q=qL,qp=qpL,name='water')

   ! ---- Gas N2 ideal gas (R=296.8 J/(kg.K) -> cv=R/(gamma-1)=742 at gamma=1.4) -------------------
   gammaG = 1.4_WP
   cvG    = 742.0_WP
   qG     = 0.0_WP
   qpG    = 0.0_WP
   call param_read('Gas Gamma', gammaG, default=gammaG)
   call param_read('Gas Cv',    cvG,    default=cvG)
   call param_read('Gas q',     qG,     default=qG)
   call param_read('Gas qp',    qpG,    default=qpG)
   call nitrogen%initialize(gamma=gammaG,cv=cvG,q=qG,qp=qpG,name='N2')

   ! ---- Geometry and phase fractions ------------------------------------------------------------
   Vtot  = 1.0_WP
   alpha = 0.5_WP
   call param_read('Volume', Vtot, default=Vtot)
   call param_read('Liquid volume fraction', alpha, default=alpha)
   if (Vtot.le.0.0_WP) call die('[thermo0D] Volume must be > 0')
   if (alpha.lt.0.0_WP.or.alpha.gt.1.0_WP) call die('[thermo0D] Liquid volume fraction must be in [0,1]')
   has_liq = (alpha.gt.VF_PURE)
   has_gas = (alpha.lt.1.0_WP-VF_PURE)
   if (alpha.le.VF_PURE) alpha = 0.0_WP
   if (alpha.ge.1.0_WP-VF_PURE) alpha = 1.0_WP

   ! ---- Initial pressures / temperatures (shared defaults, optional per-phase overrides) --------
   P0 = 101325.0_WP
   T0 = 298.15_WP
   call param_read('Pressure',    P0, default=P0)
   call param_read('Temperature', T0, default=T0)
   ! Per-phase overrides (defaults must not alias the intent(out) argument)
   call param_read('Liquid pressure',    PL0, default=P0)
   call param_read('Liquid temperature', TL0, default=T0)
   call param_read('Gas pressure',       PG0, default=P0)
   call param_read('Gas temperature',    TG0, default=T0)
   if (has_liq.and.(PL0.le.0.0_WP.or.TL0.le.0.0_WP)) call die('[thermo0D] Liquid P/T must be > 0')
   if (has_gas.and.(PG0.le.0.0_WP.or.TG0.le.0.0_WP)) call die('[thermo0D] Gas P/T must be > 0')

   ! Phasic densities and energies from (p,T)
   rhoL = 0.0_WP; eL = 0.0_WP; massL = 0.0_WP
   rhoG = 0.0_WP; eG = 0.0_WP; massG = 0.0_WP
   if (has_liq) then
      rhoL  = water%get_rho_from_p_T(p=PL0,T=TL0,y=y)
      eL    = water%get_e_from_p_rho(p=PL0,rho=rhoL,y=y)
      massL = rhoL*alpha*Vtot
   end if
   if (has_gas) then
      rhoG  = nitrogen%get_rho_from_p_T(p=PG0,T=TG0,y=y)
      eG    = nitrogen%get_e_from_p_rho(p=PG0,rho=rhoG,y=y)
      massG = rhoG*(1.0_WP-alpha)*Vtot
   end if
   mass = massL+massG
   if (mass.le.0.0_WP) call die('[thermo0D] Total mass must be > 0')

   ! Conserved multiphase state (unit-volume partial densities)
   Q(1) = rhoL*alpha
   Q(2) = rhoG*(1.0_WP-alpha)
   Q(3) = Q(1)*eL
   Q(4) = Q(2)*eG

   brho = 0.0_WP
   if (has_liq) then
      brho = water%b*rhoL
      if (brho.ge.water%brhomax) then
         write(msg,'(a,es12.5,a,es12.5)') &
              '[thermo0D] Initial b*rho=',brho,' exceeds packing clamp brhomax=',water%brhomax
         call warn(msg)
      end if
   end if

   ! ---- Relaxation (mixtures only) --------------------------------------------------------------
   relax_type = 'pT'
   call param_read('Relaxation type', relax_type, default=relax_type)
   relax_type = lowercase(relax_type)
   do_relax = has_liq.and.has_gas
   if (do_relax) then
      call relaxer%initialize(gas=nitrogen,liq=water)
      select case (trim(relax_type))
      case ('none','off','.false.','false')
         do_relax = .false.
      case ('p','mechanical')
         relaxer%model = Prelax
      case ('pt','pT','thermal')
         relaxer%model = PTrelax
      case default
         call die('[thermo0D] Relaxation type must be none, p, or pT')
      end select
   end if

   ! ---- Time integration / heating --------------------------------------------------------------
   heating_rate = 1.0e4_WP
   t_end = 200.0_WP
   dt    = 0.1_WP
   call param_read('Heating rate', heating_rate, default=heating_rate)  ! J/(kg.s)
   call param_read('End time',     t_end,        default=t_end)
   call param_read('Time step',    dt,           default=dt)
   if (dt.le.0.0_WP)    call die('[thermo0D] Time step must be > 0')
   if (t_end.lt.0.0_WP) call die('[thermo0D] End time must be >= 0')

   ! ---- Critical thresholds ---------------------------------------------------------------------
   Tcrit = Tcrit_H2O
   Pcrit = Pcrit_H2O
   brhomax_warn = water%brhomax
   call param_read('Critical temperature', Tcrit, default=Tcrit)
   call param_read('Critical pressure',    Pcrit, default=Pcrit)
   call param_read('Packing warn fraction', brhomax_warn, default=brhomax_warn)
   stop_at_critical = .false.
   call param_read('Stop at critical', stop_at_critical, default=stop_at_critical)

   outfile = 'data/thermo.csv'
   call param_read('Output file', outfile, default=outfile)
   nout = 1
   call param_read('Output every', nout, default=nout)
   if (nout.lt.1) nout = 1

   ! ---- Banner ----------------------------------------------------------------------------------
   if (am_root) then
      write(*,'(a)') 'thermo0D: constant-volume heat addition'
      call log('thermo0D: constant-volume heat addition')
      write(msg,'(a,es12.5,a)') '  Volume                 = ',Vtot,' m^3'; call log(msg); write(*,'(a)') trim(msg)
      write(msg,'(a,es12.5)')   '  Liquid volume fraction = ',alpha;          call log(msg); write(*,'(a)') trim(msg)
      if (has_liq) then
         write(msg,'(a,es12.5,a,es12.5,a)') '  Liquid P0,T0           = ',PL0,' Pa, ',TL0,' K'
         call log(msg); write(*,'(a)') trim(msg)
         write(msg,'(a,es12.5,a)') '  rhoL                   = ',rhoL,' kg/m^3'; call log(msg); write(*,'(a)') trim(msg)
         write(msg,'(a,es12.5,a)') '  massL                  = ',massL,' kg';    call log(msg); write(*,'(a)') trim(msg)
      end if
      if (has_gas) then
         write(msg,'(a,es12.5,a,es12.5,a)') '  Gas P0,T0              = ',PG0,' Pa, ',TG0,' K'
         call log(msg); write(*,'(a)') trim(msg)
         write(msg,'(a,es12.5,a)') '  rhoG                   = ',rhoG,' kg/m^3'; call log(msg); write(*,'(a)') trim(msg)
         write(msg,'(a,es12.5,a)') '  massG                  = ',massG,' kg';    call log(msg); write(*,'(a)') trim(msg)
      end if
      write(msg,'(a,es12.5,a)') '  mass (total)           = ',mass,' kg';           call log(msg); write(*,'(a)') trim(msg)
      write(msg,'(a,es12.5,a)') '  Heating rate           = ',heating_rate,' J/(kg.s)'; call log(msg); write(*,'(a)') trim(msg)
      write(msg,'(a,es12.5,a)') '  End time               = ',t_end,' s';            call log(msg); write(*,'(a)') trim(msg)
      write(msg,'(a,es12.5,a)') '  Time step              = ',dt,' s';               call log(msg); write(*,'(a)') trim(msg)
      if (has_liq.and.has_gas) then
         write(msg,'(2x,a,a)') 'Relaxation type         = ',trim(relax_type); call log(msg); write(*,'(a)') trim(msg)
      end if
      write(msg,'(2x,a,a)') 'Output file             = ',trim(outfile); call log(msg); write(*,'(a)') trim(msg)
      call water%print()
      call nitrogen%print()
   end if

   ! ---- Open CSV --------------------------------------------------------------------------------
   if (am_root) then
      open(newunit=iunit,file=trim(outfile),form='formatted',status='replace',action='write',iostat=ierr)
      if (ierr.ne.0) call die('[thermo0D] Could not open output file: '//trim(outfile))
      write(iunit,'(a)') '# thermo0D constant-volume heat addition'
      write(iunit,'(a,es16.8)') '# Volume = ',Vtot
      write(iunit,'(a,es16.8)') '# Liquid volume fraction = ',alpha
      write(iunit,'(a,es16.8)') '# Heating rate [J/(kg.s)] = ',heating_rate
      write(iunit,'(a,es16.8)') '# Tcrit = ',Tcrit
      write(iunit,'(a,es16.8)') '# Pcrit = ',Pcrit
      write(iunit,'(a)') 't,alpha,eL,eG,e,PL,PG,P,TL,TG,T,rhoL,rhoG,rho,'// &
           'vL,vG,hL,hG,sL,sG,gL,gG,cL,cG,cvL,cvG,massL,massG,mass,V,Q_added,b_rho,'// &
           'flag_Tcrit,flag_Pcrit,flag_critical,flag_packing'
   end if

   ! ---- Time loop -------------------------------------------------------------------------------
   time = 0.0_WP
   Q_added = 0.0_WP
   nstep = 0
   hit_Tcrit = .false.; hit_Pcrit = .false.; hit_critical = .false.; hit_packing = .false.
   announced_T = .false.; announced_P = .false.; announced_pack = .false.

   do
      call update_thermo()

      ! Critical flags from liquid when present, else gas
      if (has_liq) then
         if (TL.ge.Tcrit) hit_Tcrit = .true.
         if (PL.ge.Pcrit) hit_Pcrit = .true.
         if (brho.ge.brhomax_warn) hit_packing = .true.
      else
         if (TG.ge.Tcrit) hit_Tcrit = .true.
         if (PG.ge.Pcrit) hit_Pcrit = .true.
      end if
      if (hit_Tcrit.and.hit_Pcrit) hit_critical = .true.

      flag_Tcrit     = merge(1,0,hit_Tcrit)
      flag_Pcrit     = merge(1,0,hit_Pcrit)
      flag_critical  = merge(1,0,hit_critical)
      flag_packing   = merge(1,0,hit_packing)

      if (am_root) then
         if (hit_Tcrit.and.(.not.announced_T)) then
            write(msg,'(a,es12.5,a,es12.5,a)') &
                 '  >>> Reached critical temperature at t=',time,' s  (T=',Temp,' K)'
            call log(msg); write(*,'(a)') trim(msg)
            announced_T = .true.
         end if
         if (hit_Pcrit.and.(.not.announced_P)) then
            write(msg,'(a,es12.5,a,es12.5,a)') &
                 '  >>> Reached critical pressure at t=',time,' s  (P=',P,' Pa)'
            call log(msg); write(*,'(a)') trim(msg)
            announced_P = .true.
         end if
         if (hit_packing.and.(.not.announced_pack)) then
            write(msg,'(a,es12.5,a,es12.5)') &
                 '  >>> Packing warning: b*rho=',brho,' >= ',brhomax_warn
            call log(msg); write(*,'(a)') trim(msg)
            announced_pack = .true.
         end if
      end if

      if (am_root .and. mod(nstep,nout).eq.0) call write_row()

      if (time.ge.t_end) exit
      if (stop_at_critical.and.hit_critical) then
         if (am_root) then
            call log('  Stopping: both critical T and P thresholds reached.')
            write(*,'(a)') '  Stopping: both critical T and P thresholds reached.'
         end if
         exit
      end if
      if (Temp.le.0.0_WP .or. P.ne.P) then
         if (am_root) then
            call warn('[thermo0D] Non-physical thermo state; aborting time loop.')
            write(*,'(a)') '[thermo0D] Non-physical thermo state; aborting time loop.'
         end if
         exit
      end if

      ! Heat: same specific rate on each present phase  (de_k/dt = heating_rate)
      if (has_liq) then
         eL = eL + heating_rate*dt
         Q(3) = Q(1)*eL
      end if
      if (has_gas) then
         eG = eG + heating_rate*dt
         Q(4) = Q(2)*eG
      end if
      Q_added = Q_added + heating_rate*dt

      ! Mechanical / thermal relaxation for mixtures (conserves mass and total U)
      if (do_relax) then
         call relaxer%apply(dt=dt,VF=alpha,Q=Q,Pjump=0.0_WP,ierr=relax_ierr)
         if (relax_ierr.ne.RELAX_OK .and. am_root .and. nstep.eq.0) then
            write(msg,'(a,i0)') '[thermo0D] relaxation returned status ',relax_ierr
            call warn(msg)
         end if
         ! Refresh phase presence after VF update
         has_liq = (alpha.gt.VF_PURE)
         has_gas = (alpha.lt.1.0_WP-VF_PURE)
         if (has_liq) then
            eL = Q(3)/Q(1)
         end if
         if (has_gas) then
            eG = Q(4)/Q(2)
         end if
      end if

      time = time + dt
      nstep = nstep + 1
   end do

   if (am_root) then
      if (mod(nstep,nout).ne.0) then
         call update_thermo()
         call write_row()
      end if
      close(iunit)
      write(msg,'(a,i0,a)') 'Done. Wrote ',nstep/nout+1,' rows to '//trim(outfile)
      call log(msg); write(*,'(a)') trim(msg)
      write(msg,'(a,es12.5,a,es12.5,a,es12.5)') 'Final state: T=',Temp,' K,  P=',P,' Pa,  alpha=',alpha
      call log(msg); write(*,'(a)') trim(msg)
      write(msg,'(a,l1,a,l1,a,l1)') &
           'Flags: Tcrit=',hit_Tcrit,'  Pcrit=',hit_Pcrit,'  packing=',hit_packing
      call log(msg); write(*,'(a)') trim(msg)
   end if

   call water%finalize()
   call nitrogen%finalize()
   call param_final
   call messager_final
   call parallel_final

contains

   subroutine update_thermo()
      implicit none
      ! Recover phasic state from (alpha,Q) or pure-phase (rho,e)
      PL=0.0_WP; PG=0.0_WP; TL=0.0_WP; TG=0.0_WP
      hL=0.0_WP; hG=0.0_WP; sL=0.0_WP; sG=0.0_WP; gL=0.0_WP; gG=0.0_WP
      cL=0.0_WP; cG=0.0_WP; cvL_loc=0.0_WP; cvG_loc=0.0_WP
      vL=0.0_WP; vG=0.0_WP; brho=0.0_WP
      rhoL=0.0_WP; rhoG=0.0_WP

      if (has_liq) then
         if (alpha.gt.VF_PURE) then
            rhoL = Q(1)/alpha
            eL   = Q(3)/Q(1)
         end if
         PL      = water%get_p_from_rho_e(rho=rhoL,e=eL,y=y)
         TL      = water%get_T_from_rho_e(rho=rhoL,e=eL,y=y)
         cL      = water%get_c_from_rho_e(rho=rhoL,e=eL,y=y)
         cvL_loc = water%get_cv_from_rho_e(rho=rhoL,e=eL,y=y)
         hL      = water%get_h_from_p_T(p=PL,T=TL,y=y)
         sL      = water%get_s_from_p_T(p=PL,T=TL,y=y)
         gL      = water%get_g_from_p_T(p=PL,T=TL,y=y)
         vL      = 1.0_WP/rhoL
         brho    = water%b*rhoL
         massL   = Q(1)*Vtot
      else
         massL = 0.0_WP
         eL = 0.0_WP
      end if

      if (has_gas) then
         if (alpha.lt.1.0_WP-VF_PURE) then
            rhoG = Q(2)/(1.0_WP-alpha)
            eG   = Q(4)/Q(2)
         end if
         PG      = nitrogen%get_p_from_rho_e(rho=rhoG,e=eG,y=y)
         TG      = nitrogen%get_T_from_rho_e(rho=rhoG,e=eG,y=y)
         cG      = nitrogen%get_c_from_rho_e(rho=rhoG,e=eG,y=y)
         cvG_loc = nitrogen%get_cv_from_rho_e(rho=rhoG,e=eG,y=y)
         hG      = nitrogen%get_h_from_p_T(p=PG,T=TG,y=y)
         sG      = nitrogen%get_s_from_p_T(p=PG,T=TG,y=y)
         gG      = nitrogen%get_g_from_p_T(p=PG,T=TG,y=y)
         vG      = 1.0_WP/rhoG
         massG   = Q(2)*Vtot
      else
         massG = 0.0_WP
         eG = 0.0_WP
      end if

      mass = massL+massG
      rho_mix = mass/Vtot
      e_mix = 0.0_WP
      if (mass.gt.0.0_WP) e_mix = (massL*eL+massG*eG)/mass

      ! Mixture P,T for flags / summary
      if (has_liq.and.has_gas) then
         P = 0.5_WP*(PL+PG)
         Temp = 0.5_WP*(TL+TG)
      else if (has_liq) then
         P = PL; Temp = TL
      else
         P = PG; Temp = TG
      end if
   end subroutine update_thermo

   subroutine write_row()
      implicit none
      write(iunit,'(31(es16.8,","),es16.8,",",3(i0,","),i0)') &
           time,alpha,eL,eG,e_mix,PL,PG,P,TL,TG,Temp,rhoL,rhoG,rho_mix, &
           vL,vG,hL,hG,sL,sG,gL,gG,cL,cG,cvL_loc,cvG_loc,massL,massG,mass,Vtot,Q_added,brho, &
           flag_Tcrit,flag_Pcrit,flag_critical,flag_packing
   end subroutine write_row

end program thermo0D
