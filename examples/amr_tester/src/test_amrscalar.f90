!> Test amrscalar solver with time integration and Ensight output
module mod_test_amrscalar
   use precision,         only: WP
   use amrgrid_class,     only: amrgrid
   use amrscalar_class,   only: amrscalar
   use amrensight_class,  only: amrensight
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use messager,          only: log,warn
   use amrex_amr_module,  only: amrex_boxarray,amrex_distromap,amrex_multifab
   implicit none
   private
   public :: test_amrscalar

   ! Module-level pointers for callbacks
   type(amrgrid), pointer :: amr_ptr=>null()
   type(amrscalar), pointer :: sc_ptr=>null()
   real(WP), parameter :: PI=3.14159265358979323846_WP
   real(WP), parameter :: SC_REFINE_THRESH=0.01_WP  !< Refine where SC > this value

contains

   !> Geometric tagger - refine center of domain
   !> NOTE: SC-based tagging requires that the tags BoxArray and SC MultiFab
   !> BoxArray be identical. During regrid, they may differ, causing crashes.
   !> Future improvement: implement coordinate-based lookup into SC data.
   subroutine box_tagger(lvl,tags_ptr,time)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_tagboxarray,amrex_mfiter,amrex_box
      implicit none
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags_ptr
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      character(kind=c_char), parameter :: SET=char(1)
      real(WP) :: x,y,z,dx,dy,dz
      integer :: i,j,k
      tags=tags_ptr
      dx=amr_ptr%geom(lvl)%dx(1)
      dy=amr_ptr%geom(lvl)%dx(2)
      dz=amr_ptr%geom(lvl)%dx(3)
      call amr_ptr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tags%dataPtr(mfi)
         do k=bx%lo(3),bx%hi(3)
            z=amr_ptr%zlo+(real(k,WP)+0.5_WP)*dz
            do j=bx%lo(2),bx%hi(2)
               y=amr_ptr%ylo+(real(j,WP)+0.5_WP)*dy
               do i=bx%lo(1),bx%hi(1)
                  x=amr_ptr%xlo+(real(i,WP)+0.5_WP)*dx
                  ! Refine center of domain
                  if (x>0.25_WP.and.x<0.75_WP.and.&
                     y>0.25_WP.and.y<0.75_WP.and.&
                     z>0.25_WP.and.z<0.75_WP) then
                     tagarr(i,j,k,1)=SET
                  end if
               end do
            end do
         end do
      end do
      call amr_ptr%mfiter_destroy(mfi)
   end subroutine box_tagger


   subroutine test_amrscalar()
      use iso_c_binding,    only: c_associated
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      type(amrgrid), target :: amr
      type(amrscalar), target :: sc
      type(amrensight) :: ens
      type(timetracker) :: time
      type(event) :: regrid_evt,ensight_evt
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: U,V,W,dSCdt,SCfill
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSC,pU,pV,pW
      integer :: lvl,i,j,k
      real(WP) :: x,y,z,dx,dy,dz,int0,intF

      call log("---------------------------------------------------")
      call log("Running Test: amrscalar Time Integration")
      call log("---------------------------------------------------")

      ! Setup grid (periodic box)
      amr%nx  =32
      amr%ny  =32
      amr%nz  =32
      amr%xlo =0.0_WP
      amr%xhi =1.0_WP
      amr%ylo =0.0_WP
      amr%yhi =1.0_WP
      amr%zlo =0.0_WP
      amr%zhi =1.0_WP
      amr%xper=.true.
      amr%yper=.true.
      amr%zper=.true.
      amr%nlvl=1   ! 2 levels total (0,1)
      amr%nmax=16

      call amr%initialize("sc_amr")
      amr_ptr=>amr

      ! Initialize scalar solver first (registers callbacks)
      call sc%initialize(amr,nscalar=1,name="test_scalar")
      sc_ptr=>sc

      ! Set tagging - now sc_ptr is available for gradient-based refinement
      call amr%add_tagging(box_tagger)

      ! Build grid
      call amr%initialize_grid(0.0_WP)
      call amr%get_info()
      call log("Grid built: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")

      ! Initialize ensight output
      call ens%initialize(amr=amr,name="sc_advect")
      call ens%add_scalar(data=sc%SC,comp=1,name="SC")

      ! Initialize scalar with Gaussian blob centered at (0.5, 0.5, 0.5)
      do lvl=0,amr%clvl()
         dx=amr%geom(lvl)%dx(1)
         dy=amr%geom(lvl)%dx(2)
         dz=amr%geom(lvl)%dx(3)
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pSC=>sc%SC%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3)
               z=amr%zlo+(real(k,WP)+0.5_WP)*dz
               do j=bx%lo(2),bx%hi(2)
                  y=amr%ylo+(real(j,WP)+0.5_WP)*dy
                  do i=bx%lo(1),bx%hi(1)
                     x=amr%xlo+(real(i,WP)+0.5_WP)*dx
                     ! Gaussian blob
                     pSC(i,j,k,1)=exp(-50.0_WP*((x-0.5_WP)**2+(y-0.5_WP)**2+(z-0.5_WP)**2))
                  end do
               end do
            end do
         end do
         call amr%mfiter_destroy(mfi)
      end do

      ! Regrid now that we have data - creates refined levels where SC > threshold
      call amr%regrid(baselvl=0,time=0.0_WP)
      call amr%get_info()
      call log("After initial regrid: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")

      ! Average down to make levels consistent (only if we have fine levels)
      if (amr%clvl().ge.1) call amr%average_downto(sc%SC%mf,0)

      ! Get initial integral
      call sc%get_info()
      int0=sc%SCint(1)
      call log("Initial: SCint="//trim(rtoa(int0)))

      ! Write initial output
      call ens%write_data(time=0.0_WP)

      ! Setup timetracker
      time=timetracker(amRoot=amr%amRoot)
      time%dt=0.005_WP
      time%dtmax=time%dt
      time%tmax=1.0_WP
      time%t=0.0_WP
      time%n=0

      ! Setup regrid event (every 10 steps)
      regrid_evt=event(time=time,name='Regrid')
      regrid_evt%nper=10

      ! Setup ensight event (every 20 steps)
      ensight_evt=event(time=time,name='Ensight')
      ensight_evt%nper=20

      call log("Advancing to t="//trim(rtoa(time%tmax))//", dt="//trim(rtoa(time%dt)))

      ! Time loop
      do while (.not.time%done())
         call time%increment()

         ! Copy to old
         call sc%copy2old()

         ! Loop over levels (uniform dt)
         do lvl=0,amr%clvl()
            dx=amr%dx(lvl)
            dy=amr%dy(lvl)
            dz=amr%dz(lvl)

            ! Build velocity (rotating vortex in xy plane)
            call amr%mfab_build(lvl=lvl,mfab=U,ncomp=1,nover=0,atface=[.true.,.false.,.false.])
            call amr%mfab_build(lvl=lvl,mfab=V,ncomp=1,nover=0,atface=[.false.,.true.,.false.])
            call amr%mfab_build(lvl=lvl,mfab=W,ncomp=1,nover=0,atface=[.false.,.false.,.true.])
            call amr%mfab_build(lvl=lvl,mfab=dSCdt,ncomp=1,nover=0,atface=[.false.,.false.,.false.])
            call amr%mfab_build(lvl=lvl,mfab=SCfill,ncomp=1,nover=2,atface=[.false.,.false.,.false.])

            ! Set velocity: rotating vortex in z-plane (U=-y+0.5, V=x-0.5, W=0)
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               bx=mfi%tilebox()
               pU=>U%dataptr(mfi)
               pV=>V%dataptr(mfi)
               pW=>W%dataptr(mfi)
               ! U on x-faces
               do k=lbound(pU,3),ubound(pU,3)
                  do j=lbound(pU,2),ubound(pU,2)
                     y=amr%ylo+(real(j,WP)+0.5_WP)*dy
                     do i=lbound(pU,1),ubound(pU,1)
                        pU(i,j,k,1)=-(y-0.5_WP)*cos(PI*time%t/time%tmax)
                     end do
                  end do
               end do
               ! V on y-faces
               do k=lbound(pV,3),ubound(pV,3)
                  do j=lbound(pV,2),ubound(pV,2)
                     do i=lbound(pV,1),ubound(pV,1)
                        x=amr%xlo+(real(i,WP)+0.5_WP)*dx
                        pV(i,j,k,1)=(x-0.5_WP)*cos(PI*time%t/time%tmax)
                     end do
                  end do
               end do
               ! W on z-faces
               pW=0.0_WP
            end do
            call amr%mfiter_destroy(mfi)

            ! Copy interior from SCold and fill ghost cells
            call SCfill%copy(sc%SCold%mf(lvl),1,1,1,0)
            call SCfill%fill_boundary(amr%geom(lvl))

            ! Calculate dSC/dt
            call sc%get_dSCdt(lvl,dSCdt,SCfill,U,V,W)

            ! Forward Euler step: SC = SCold + dt*dSC/dt
            call sc%SC%mf(lvl)%lincomb(a=1.0_WP,srcmf1=sc%SCold%mf(lvl),srccomp1=1,&
            &                          b=time%dt,srcmf2=dSCdt            ,srccomp2=1,&
            &                          dstcomp=1,nc=1,ng=0)

            ! Cleanup level work arrays
            call amr%mfab_destroy(U)
            call amr%mfab_destroy(V)
            call amr%mfab_destroy(W)
            call amr%mfab_destroy(dSCdt)
            call amr%mfab_destroy(SCfill)
         end do

         ! Reflux and average down
         call sc%reflux_avg(time%dt)

         ! Regrid if needed
         if (regrid_evt%occurs()) then
            call log("Regridding at step "//trim(itoa(time%n)))
            call amr%regrid(baselvl=0,time=time%t)
            call amr%get_info()
            if (amr%clvl().ge.1) call amr%average_downto(sc%SC%mf,0)
            call log("  Grid: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")
         end if

         ! Write ensight output
         if (ensight_evt%occurs().or.time%done()) then
            call ens%write_data(time=time%t)
            call sc%get_info()
            call log("Step "//trim(itoa(time%n))//": t="//trim(rtoa(time%t))//", SCint="//trim(rtoa(sc%SCint(1))))
         end if
      end do

      ! Final conservation check
      call sc%get_info()
      intF=sc%SCint(1)
      call log("Final: SCint="//trim(rtoa(intF)))
      call log("Conservation error: "//trim(rtoa(abs(intF-int0)/int0*100.0_WP))//"%")

      if (abs(intF-int0)/int0<1.0e-10_WP) then
         call log("PASS: Scalar conserved to machine precision")
      else
         call warn("FAIL: Conservation error too large!")
      end if

      ! Cleanup
      amr_ptr=>null()
      sc_ptr=>null()
      call sc%finalize()
      call amr%finalize()
      call log("PASS: amrscalar test complete!")
      call log("View output in ensight/sc_advect/")

   contains
      function itoa(ii) result(str)
         integer, intent(in) :: ii
         character(len=32) :: str
         write(str,'(i0)') ii
      end function itoa
      function rtoa(r) result(str)
         real(WP), intent(in) :: r
         character(len=32) :: str
         write(str,'(es10.3)') r
      end function rtoa
   end subroutine test_amrscalar

end module mod_test_amrscalar
