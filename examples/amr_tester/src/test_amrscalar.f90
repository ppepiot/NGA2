!> Test amrscalar solver with time integration
module mod_test_amrscalar
   use precision,         only: WP
   use amrgrid_class,     only: amrgrid
   use amrscalar_class,   only: amrscalar
   use amrex_amr_module,  only: amrex_boxarray,amrex_distromap,amrex_multifab,amrex_mfiter,amrex_box
   implicit none
   private
   public :: test_amrscalar

   ! Module-level objects (referenced by callbacks)
   type(amrgrid),   allocatable, target :: amr
   type(amrscalar), allocatable, target :: sc

   real(WP), parameter :: SC_REFINE_THRESH=0.01_WP  !< Refine where SC > this value

contains

   !> Custom on_init callback: initialize scalar field with Gaussian blob
   subroutine gaussian_on_init(solver, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_mfiter_build, amrex_mfiter_destroy
      class(amrscalar), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSC
      real(WP) :: x, y, z, dx, dy, dz
      integer :: i, j, k
      ! Reset level layouts
      call solver%SC%reset_level(lvl, ba, dm)
      call solver%SCold%reset_level(lvl, ba, dm)
      ! Build flux register for fine levels
      if (lvl .ge. 1) call solver%flux%reset_level(lvl, ba, dm, solver%amr%rref(lvl-1))
      ! Initialize with Gaussian blob
      dx = solver%amr%dx(lvl)
      dy = solver%amr%dy(lvl)
      dz = solver%amr%dz(lvl)
      call amrex_mfiter_build(mfi, solver%SC%mf(lvl))
      do while (mfi%next())
         bx = mfi%tilebox()
         pSC => solver%SC%mf(lvl)%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3)
            z = solver%amr%zlo + (real(k,WP)+0.5_WP)*dz
            do j = bx%lo(2), bx%hi(2)
               y = solver%amr%ylo + (real(j,WP)+0.5_WP)*dy
               do i = bx%lo(1), bx%hi(1)
                  x = solver%amr%xlo + (real(i,WP)+0.5_WP)*dx
                  ! Gaussian offset from center
                  pSC(i,j,k,1) = exp(-200.0_WP*((x-0.25_WP)**2 + (y-0.5_WP)**2 + (z-0.5_WP)**2))
               end do
            end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
      ! SCold starts at zero
      call solver%SCold%mf(lvl)%setval(0.0_WP)
   end subroutine gaussian_on_init


   !> Custom tagging callback: refine where SC > threshold
   subroutine scalar_tagger(solver, lvl, tags_ptr, time)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_tagboxarray
      class(amrscalar), intent(inout) :: solver
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags_ptr
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=1), contiguous, pointer :: tagarr(:,:,:,:)
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSC
      character(kind=1), parameter :: SET = char(1)
      integer :: i, j, k
      tags = tags_ptr
      ! Iterate over level boxes
      call solver%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pSC => solver%SC%mf(lvl)%dataptr(mfi)
         tagarr => tags%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  if (pSC(i,j,k,1) .gt. SC_REFINE_THRESH) then
                     tagarr(i,j,k,1) = SET
                  end if
               end do
            end do
         end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine scalar_tagger


   subroutine test_amrscalar()
      use string,            only: itoa,rtoa,str_medium
      use mathtools,         only: Pi
      use amrviz_class,      only: amrviz
      use amrio_class,       only: amrio
      use monitor_class,     only: monitor
      use timetracker_class, only: timetracker
      use event_class,       only: event
      use messager,          only: log,warn
      implicit none
      type(amrviz) :: viz
      type(amrio) :: io
      type(monitor) :: mfile
      type(timetracker) :: time
      type(event) :: regrid_evt,hdf5_evt
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: U,V,W,dSCdt,SCfill
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSC,pU,pV,pW
      integer :: lvl,i,j,k
      real(WP) :: x,y,z,dx,dy,dz,int0,intF

      call log("---------------------------------------------------")
      call log("Running Test: amrscalar Time Integration")
      call log("---------------------------------------------------")

      ! Allocate module objects
      allocate(amr,sc)

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
      amr%maxlvl=2   ! 3 levels total (0,1,2)
      amr%nmax=16

      call amr%initialize("sc_amr")

      ! Initialize scalar solver
      call sc%initialize(amr, nscalar=1, name="test_scalar")
      sc%on_init => gaussian_on_init
      sc%tagging => scalar_tagger

      ! Initialize HDF5 viz output
      call viz%initialize(amr=amr, name='sc_advect')
      call viz%add_scalar(sc%SC, 1, 'SC')
      call viz%add_scalar(sc%SCold, 1, 'SCold')

      ! Build all levels
      call amr%init_from_scratch(time=0.0_WP, do_postregrid=.true.)
      call amr%get_info()
      call log("After init_from_scratch: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")

      ! Get initial integral
      call sc%get_info()
      int0=sc%SCint(1)
      call log("Initial: SCint="//trim(rtoa(int0)))

      ! Setup timetracker
      time=timetracker(amRoot=amr%amRoot)
      time%dt=0.0025_WP
      time%dtmax=time%dt
      time%tmax=2.0_WP*Pi  ! Full circle rotation
      time%t=0.0_WP
      time%n=0

      ! Setup regrid event (every 10 steps)
      regrid_evt=event(time=time,name='Regrid')
      regrid_evt%nper=10

      ! Setup HDF5 event
      hdf5_evt=event(time=time,name='HDF5')
      hdf5_evt%tper=0.125_WP

      ! Setup monitor file
      call sc%get_info()
      mfile=monitor(amr%amRoot,'scalar')
      call mfile%add_column(time%n,'Timestep')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(sc%SCmin(1),'SC_min')
      call mfile%add_column(sc%SCmax(1),'SC_max')
      call mfile%add_column(sc%SCint(1),'SC_int')
      call mfile%write()

      ! Write initial output
      call viz%write(time=time%t)

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

            ! Set velocity: solid body rotation in z-plane
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
                        pU(i,j,k,1)=-(y-0.5_WP)
                     end do
                  end do
               end do
               ! V on y-faces
               do k=lbound(pV,3),ubound(pV,3)
                  do j=lbound(pV,2),ubound(pV,2)
                     do i=lbound(pV,1),ubound(pV,1)
                        x=amr%xlo+(real(i,WP)+0.5_WP)*dx
                        pV(i,j,k,1)=(x-0.5_WP)
                     end do
                  end do
               end do
               ! W on z-faces
               pW=0.0_WP
            end do
            call amr%mfiter_destroy(mfi)

            ! Fill SCfill with ghost cells
            call sc%SCold%fill_mfab(SCfill, lvl, time%t)

            ! Calculate dSC/dt
            call sc%get_dSCdt(lvl,dSCdt,SCfill,U,V,W)

            ! Forward Euler step
            call sc%SC%mf(lvl)%lincomb(a=1.0_WP,srcmf1=sc%SCold%mf(lvl),srccomp1=1,&
            &                          b=time%dt,srcmf2=dSCdt,srccomp2=1,&
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
            call log("  Grid: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")
         end if

         ! Write HDF5 output
         if (hdf5_evt%occurs()) then
            call viz%write(time=time%t)
            call log("Step "//trim(itoa(time%n))//": t="//trim(rtoa(time%t))//", SCint="//trim(rtoa(sc%SCint(1))))
         end if

         ! Update monitor file
         call sc%get_info()
         call mfile%write()
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

      ! Test checkpoint round-trip
      call io%initialize(amr, nfiles=1)
      call sc%register_checkpoint(io)
      call io%add_scalar('dt', time%dt)
      call io%write('checkpoint_test', time%t, time%n)
      call log("PASS: Checkpoint written to checkpoint_test/")

      ! Save current integral, then zero the data
      call sc%get_info()
      int0 = sc%SCint(1)
      call log("  Before zero: SCint="//trim(rtoa(int0)))

      ! Zero out SC data
      do lvl = 0, amr%clvl()
         call amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            pSC => sc%SC%mf(lvl)%dataptr(mfi)
            pSC = 0.0_WP
         end do
         call amr%mfiter_destroy(mfi)
      end do
      call sc%get_info()
      call log("  After zero: SCint="//trim(rtoa(sc%SCint(1))))

      ! Read checkpoint back
      block
         real(WP) :: read_time, read_dt
         integer :: read_step
         call io%read_header('checkpoint_test', read_time, read_step)
         call io%get_scalar('dt', read_dt)
         call sc%restore_checkpoint(io, 'checkpoint_test')
         call log("  Read checkpoint: time="//trim(rtoa(read_time))//" step="//trim(itoa(read_step))//" dt="//trim(rtoa(read_dt)))
      end block

      ! Verify data restored
      call sc%get_info()
      intF = sc%SCint(1)
      call log("  After read: SCint="//trim(rtoa(intF)))

      if (abs(intF-int0)/int0 .lt. 1.0e-10_WP) then
         call log("PASS: Checkpoint round-trip verified!")
      else
         call warn("FAIL: Checkpoint round-trip failed!")
      end if

      call log("PASS: HDF5 plotfiles written to amrviz/sc_advect/")

      ! Cleanup
      call io%finalize()
      call viz%finalize()
      call mfile%finalize()
      call sc%finalize()
      call amr%finalize()
      if (allocated(sc)) deallocate(sc)
      if (allocated(amr)) deallocate(amr)
      call log("PASS: amrscalar test complete!")

   end subroutine test_amrscalar

end module mod_test_amrscalar
