!> AMR Visualization class: HDF5 plotfile output with field registration
!> Outputs Chombo-compatible HDF5 files with all registered fields in single file
!> readable by VisIt (Chombo reader) and ParaView (VisItChomboReader)
module amrviz_class
   use precision,      only: WP
   use string,         only: str_medium
   use iso_c_binding,  only: c_ptr,c_char,c_int,c_double,c_null_char,c_loc
   use amrgrid_class,  only: amrgrid
   use amrdata_class,  only: amrdata
   use amrex_interface, only: amrplotfile_write_hdf5,amrplotfile_read_time
   use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box, &
   &                           amrex_multifab_build,amrex_multifab_destroy, &
   &                           amrex_mfiter_build,amrex_mfiter_destroy
   implicit none
   private

   public :: amrviz

   ! Scalar field registration node
   type :: scl
      type(scl), pointer :: next => null()
      character(len=str_medium) :: name
      type(amrdata), pointer :: ptr => null()
      integer :: comp
      logical :: recenter = .true.  !< If true, interpolate to cell centers
   end type scl

   ! Vector field registration node
   type :: vct
      type(vct), pointer :: next => null()
      character(len=str_medium) :: name
      type(amrdata), pointer :: ptrx => null()
      type(amrdata), pointer :: ptry => null()
      type(amrdata), pointer :: ptrz => null()
      integer :: compx, compy, compz
      logical :: recenter = .true.  !< If true, interpolate to cell centers
   end type vct

   !> AMR visualization handler
   type :: amrviz
      class(amrgrid), pointer :: amr => null()   !< Pointer to AMR grid
      character(len=str_medium) :: name          !< Subdirectory name for output
      type(scl), pointer :: first_scl => null()  !< Registered scalars
      type(vct), pointer :: first_vct => null()  !< Registered vectors
      integer :: ntime = 0                       !< File counter for output
      real(WP), allocatable :: time(:)           !< Time values for each file
   contains
      procedure :: initialize                    !< Initialize with grid
      procedure :: add_scalar                    !< Register a scalar field
      procedure :: add_vector                    !< Register a vector field
      procedure :: write                         !< Write all registered fields to single HDF5
      procedure :: finalize                      !< Clean up registered field lists
   end type amrviz

contains

   !> Initialize visualization handler with AMR grid
   !> If output directory already exists with files, reads their times to continue series
   subroutine initialize(this, amr, name)
      use filesys, only: makedir, isdir, isfile
      use mpi_f08, only: MPI_BCAST, MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrviz), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in) :: name

      character(len=str_medium) :: filename, dirname
      integer :: ierr, n
      real(WP) :: file_time

      this%amr => amr
      this%name = trim(adjustl(name))
      this%ntime = 0
      this%first_scl => null()
      this%first_vct => null()

      ! Create output directory: amrviz/<name>/
      dirname = 'amrviz/'//trim(this%name)
      if (this%amr%amRoot) then
         if (.not.isdir('amrviz')) call makedir('amrviz')
         if (.not.isdir(trim(dirname))) call makedir(trim(dirname))
      end if

      ! Check for existing files and read their times (root only, then broadcast)
      ! Scan for any centering pattern - try each until we find files
      if (this%amr%amRoot) then
         centering_names: block
            character(len=8), parameter :: centerings(8) = &
               ['cell   ', 'xface  ', 'yface  ', 'zface  ', 'node   ', 'xyedge ', 'xzedge ', 'yzedge ']
            character(len=8) :: found_centering
            integer :: ic

            ! Find first centering type that has files
            found_centering = ''
            find_centering: do ic = 1, 8
               write(filename,'(a,"/nga2.",a,".",i6.6,".h5")') trim(dirname), trim(centerings(ic)), 1
               if (isfile(trim(filename))) then
                  found_centering = centerings(ic)
                  exit find_centering
               end if
            end do find_centering

            if (len_trim(found_centering) > 0) then
               ! Count existing files using found centering
               n = 0
               do
                  n = n + 1
                  write(filename,'(a,"/nga2.",a,".",i6.6,".h5")') trim(dirname), trim(found_centering), n
                  if (.not.isfile(trim(filename))) exit
                  file_time = amrplotfile_read_time(trim(filename)//c_null_char)
                  if (file_time .lt. 0.0_WP) exit  ! File exists but unreadable
               end do
               this%ntime = n - 1

               ! Allocate and fill time array
               if (this%ntime .gt. 0) then
                  allocate(this%time(this%ntime))
                  do n = 1, this%ntime
                     write(filename,'(a,"/nga2.",a,".",i6.6,".h5")') trim(dirname), trim(found_centering), n
                     this%time(n) = amrplotfile_read_time(trim(filename)//c_null_char)
                  end do
               end if
            end if
         end block centering_names
      end if

      ! Broadcast ntime and time array to all ranks
      call MPI_BCAST(this%ntime, 1, MPI_INTEGER, 0, this%amr%comm, ierr)
      if (this%ntime .gt. 0) then
         if (.not.this%amr%amRoot) allocate(this%time(this%ntime))
         call MPI_BCAST(this%time, this%ntime, MPI_REAL_WP, 0, this%amr%comm, ierr)
      end if
   end subroutine initialize


   !> Register a scalar field for output
   !> recenter: if .true. (default), interpolate to cell centers; if .false., use native centering
   subroutine add_scalar(this, data, comp, name, recenter)
      use messager, only: die
      implicit none
      class(amrviz), intent(inout) :: this
      type(amrdata), target, intent(in) :: data
      integer, intent(in) :: comp
      character(len=*), intent(in), optional :: name
      logical, intent(in), optional :: recenter
      type(scl), pointer :: new_scl

      ! Check that the component is meaningful
      if (comp.le.0) call die('[amrviz add_scalar] comp must be at least one')
      if (comp.gt.data%ncomp) call die('[amrviz add_scalar] comp is too large for provided amrdata')

      ! Create new scalar node
      allocate(new_scl)
      if (present(name)) then
         new_scl%name = trim(adjustl(name))
      else
         new_scl%name = trim(adjustl(data%name))
      end if
      new_scl%ptr => data
      new_scl%comp = comp
      if (present(recenter)) then
         new_scl%recenter = recenter
      else
         new_scl%recenter = .true.
      end if

      ! Insert at front of list
      new_scl%next => this%first_scl
      this%first_scl => new_scl
   end subroutine add_scalar


   !> Register a vector field for output (3 components from 3 amrdata objects)
   !> For MAC velocity: each component is recentered from its native face to cell centers
   subroutine add_vector(this, datax, compx, datay, compy, dataz, compz, name)
      use messager, only: die
      implicit none
      class(amrviz), intent(inout) :: this
      type(amrdata), target, intent(in) :: datax, datay, dataz
      integer, intent(in) :: compx, compy, compz
      character(len=*), intent(in), optional :: name
      type(vct), pointer :: new_vct

      ! Check that the components are meaningful
      if (compx.le.0) call die('[amrviz add_vector] compx must be at least one')
      if (compx.gt.datax%ncomp) call die('[amrviz add_vector] compx is too large for provided datax')
      if (compy.le.0) call die('[amrviz add_vector] compy must be at least one')
      if (compy.gt.datay%ncomp) call die('[amrviz add_vector] compy is too large for provided datay')
      if (compz.le.0) call die('[amrviz add_vector] compz must be at least one')
      if (compz.gt.dataz%ncomp) call die('[amrviz add_vector] compz is too large for provided dataz')

      ! Create new vector node
      allocate(new_vct)
      if (present(name)) then
         new_vct%name = trim(adjustl(name))
      else
         new_vct%name = 'vector'
      end if
      new_vct%ptrx => datax
      new_vct%ptry => datay
      new_vct%ptrz => dataz
      new_vct%compx = compx
      new_vct%compy = compy
      new_vct%compz = compz
      new_vct%recenter = .true.  ! Vectors always recenter to cell

      ! Insert at front of list
      new_vct%next => this%first_vct
      this%first_vct => new_vct
   end subroutine add_vector


   !> Write all registered fields to HDF5 plotfiles
   !> Fields are grouped by centering type - one file per centering
   !> File pattern: amrviz/<name>/nga2.<centering>.<ntime>.h5
   !> Centerings: cell, xface, yface, zface, xyedge, xzedge, yzedge, node
   subroutine write(this, time)
      implicit none
      class(amrviz), intent(inout) :: this
      real(WP), intent(in) :: time

      integer :: nlev, lev, nscl_grp, nvct_grp, ncomp, icomp, i, n, ngroups
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      real(WP), dimension(:), allocatable :: temp_time
      character(len=str_medium) :: filename
      character(len=8) :: suffix

      ! Centering groups: store unique nodal patterns found
      ! Max 8 possible centerings (2^3), but typically only 1-4 used
      integer, parameter :: MAX_GROUPS = 8
      logical :: group_nodal(3, MAX_GROUPS)
      integer :: group_count(MAX_GROUPS)
      logical :: found, eff_nodal(3), needs_interp, native_nodal(3)
      integer :: ii, jj, kk, di, dj, dk

      ! Working arrays for current group
      type(amrex_multifab), allocatable :: combined(:)
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: dst, src
      type(c_ptr), allocatable :: mf_ptrs(:), geom_ptrs(:)
      integer(c_int), allocatable :: level_steps(:), ref_ratios(:)
      character(len=64), target, allocatable :: varnames(:)
      type(c_ptr), allocatable :: varname_ptrs(:)

      nlev = this%amr%clvl() + 1
      ngroups = 0
      group_count = 0

      ! --- Phase 1: Collect unique centering types and count fields per group ---
      ! If recenter=.true., effective centering is cell-centered
      my_scl => this%first_scl
      do while (associated(my_scl))
         ! Determine effective centering for output
         if (my_scl%recenter) then
            eff_nodal = [.false., .false., .false.]  ! Recenter to cell
         else
            eff_nodal = my_scl%ptr%nodal             ! Native centering
         end if

         ! Check if this centering already exists in groups
         found = .false.
         do i = 1, ngroups
            if (all(eff_nodal .eqv. group_nodal(:,i))) then
               group_count(i) = group_count(i) + 1
               found = .true.
               exit
            end if
         end do
         ! New centering type
         if (.not.found) then
            ngroups = ngroups + 1
            group_nodal(:,ngroups) = eff_nodal
            group_count(ngroups) = 1
         end if
         my_scl => my_scl%next
      end do

      ! Count vectors: vectors always recenter to cell-centered
      my_vct => this%first_vct
      do while (associated(my_vct))
         eff_nodal = [.false., .false., .false.]  ! Vectors always go to cell-centered

         ! Find or create cell-centered group
         found = .false.
         do i = 1, ngroups
            if (all(eff_nodal .eqv. group_nodal(:,i))) then
               group_count(i) = group_count(i) + 3  ! Vector has 3 components
               found = .true.
               exit
            end if
         end do
         if (.not.found) then
            ngroups = ngroups + 1
            group_nodal(:,ngroups) = eff_nodal
            group_count(ngroups) = 3
         end if
         my_vct => my_vct%next
      end do

      if (ngroups == 0) return  ! Nothing to write

      ! --- Time tracking (shared across all groups) ---
      if (this%ntime.eq.0) then
         this%ntime=1
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(this%ntime))
         this%time(1)=time
      else
         n=1
         rewind: do i=this%ntime,1,-1
            if (this%time(i).lt.time-1.0e-6_WP) then
               n=i+1; exit rewind
            end if
         end do rewind
         this%ntime=n; allocate(temp_time(1:this%ntime))
         temp_time=[this%time(1:this%ntime-1),time]
         call move_alloc(temp_time,this%time)
      end if

      ! --- Phase 2: Write one file per centering group ---
      do i = 1, ngroups
         nscl_grp = group_count(i)
         ncomp = nscl_grp  ! For now, scalars only per group

         ! Get centering suffix
         suffix = centering_suffix(group_nodal(:,i))

         ! Allocate varnames for this group
         allocate(varnames(ncomp))
         allocate(varname_ptrs(ncomp))
         icomp = 0

         ! Collect varnames for matching scalars (using effective centering)
         my_scl => this%first_scl
         do while (associated(my_scl))
            if (my_scl%recenter) then
               eff_nodal = [.false., .false., .false.]
            else
               eff_nodal = my_scl%ptr%nodal
            end if
            if (all(eff_nodal .eqv. group_nodal(:,i))) then
               icomp = icomp + 1
               varnames(icomp) = trim(my_scl%name)//c_null_char
               varname_ptrs(icomp) = c_loc(varnames(icomp))
            end if
            my_scl => my_scl%next
         end do

         ! Collect varnames for vectors (cell-centered group only)
         if (all(.not.group_nodal(:,i))) then  ! Cell-centered group
            my_vct => this%first_vct
            do while (associated(my_vct))
               varnames(icomp+1) = trim(my_vct%name)//'_x'//c_null_char
               varnames(icomp+2) = trim(my_vct%name)//'_y'//c_null_char
               varnames(icomp+3) = trim(my_vct%name)//'_z'//c_null_char
               varname_ptrs(icomp+1) = c_loc(varnames(icomp+1))
               varname_ptrs(icomp+2) = c_loc(varnames(icomp+2))
               varname_ptrs(icomp+3) = c_loc(varnames(icomp+3))
               icomp = icomp + 3
               my_vct => my_vct%next
            end do
         end if
         allocate(combined(0:this%amr%clvl()))

         ! Build combined MF and copy matching scalar data
         do lev = 0, this%amr%clvl()
            call this%amr%mfab_build(lvl=lev, mfab=combined(lev), ncomp=ncomp, nover=0, atface=group_nodal(:,i))

            icomp = 0
            my_scl => this%first_scl
            do while (associated(my_scl))
               ! Compute effective nodal for matching
               if (my_scl%recenter) then
                  eff_nodal = [.false., .false., .false.]
               else
                  eff_nodal = my_scl%ptr%nodal
               end if

               if (all(eff_nodal .eqv. group_nodal(:,i))) then
                  icomp = icomp + 1
                  ! Check if interpolation is needed (native != effective)
                  needs_interp = my_scl%recenter .and. any(my_scl%ptr%nodal)

                  call amrex_mfiter_build(mfi, combined(lev), tiling=this%amr%default_tiling)
                  do while (mfi%next())
                     bx = mfi%tilebox()
                     src => my_scl%ptr%mf(lev)%dataptr(mfi)
                     dst => combined(lev)%dataptr(mfi)

                     if (needs_interp) then
                        ! Interpolate nodal/face/edge data to cell centers
                        ! Average 2^n points where n = number of nodal directions
                        native_nodal = my_scl%ptr%nodal
                        do kk = bx%lo(3), bx%hi(3)
                           do jj = bx%lo(2), bx%hi(2)
                              do ii = bx%lo(1), bx%hi(1)
                                 ! Sum all corner contributions based on nodal directions
                                 dst(ii,jj,kk,icomp) = 0.0_WP
                                 do dk = 0, merge(1,0,native_nodal(3))
                                    do dj = 0, merge(1,0,native_nodal(2))
                                       do di = 0, merge(1,0,native_nodal(1))
                                          dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) + src(ii+di, jj+dj, kk+dk, my_scl%comp)
                                       end do
                                    end do
                                 end do
                                 ! Divide by number of points (2^nodal_count)
                                 dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) / real(2**count(native_nodal), WP)
                              end do
                           end do
                        end do
                     else
                        ! Direct copy (same centering)
                        dst(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), icomp) = &
                           src(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), my_scl%comp)
                     end if
                  end do
                  call amrex_mfiter_destroy(mfi)
               end if
               my_scl => my_scl%next
            end do

            ! Copy vector components (cell-centered group only)
            if (all(.not.group_nodal(:,i))) then
               my_vct => this%first_vct
               do while (associated(my_vct))
                  ! Process X component
                  icomp = icomp + 1
                  call amrex_mfiter_build(mfi, combined(lev), tiling=this%amr%default_tiling)
                  do while (mfi%next())
                     bx = mfi%tilebox()
                     src => my_vct%ptrx%mf(lev)%dataptr(mfi)
                     dst => combined(lev)%dataptr(mfi)
                     native_nodal = my_vct%ptrx%nodal
                     do kk = bx%lo(3), bx%hi(3)
                        do jj = bx%lo(2), bx%hi(2)
                           do ii = bx%lo(1), bx%hi(1)
                              dst(ii,jj,kk,icomp) = 0.0_WP
                              do dk = 0, merge(1,0,native_nodal(3))
                                 do dj = 0, merge(1,0,native_nodal(2))
                                    do di = 0, merge(1,0,native_nodal(1))
                                       dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) + src(ii+di, jj+dj, kk+dk, my_vct%compx)
                                    end do
                                 end do
                              end do
                              dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) / real(2**count(native_nodal), WP)
                           end do
                        end do
                     end do
                  end do
                  call amrex_mfiter_destroy(mfi)

                  ! Process Y component
                  icomp = icomp + 1
                  call amrex_mfiter_build(mfi, combined(lev), tiling=this%amr%default_tiling)
                  do while (mfi%next())
                     bx = mfi%tilebox()
                     src => my_vct%ptry%mf(lev)%dataptr(mfi)
                     dst => combined(lev)%dataptr(mfi)
                     native_nodal = my_vct%ptry%nodal
                     do kk = bx%lo(3), bx%hi(3)
                        do jj = bx%lo(2), bx%hi(2)
                           do ii = bx%lo(1), bx%hi(1)
                              dst(ii,jj,kk,icomp) = 0.0_WP
                              do dk = 0, merge(1,0,native_nodal(3))
                                 do dj = 0, merge(1,0,native_nodal(2))
                                    do di = 0, merge(1,0,native_nodal(1))
                                       dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) + src(ii+di, jj+dj, kk+dk, my_vct%compy)
                                    end do
                                 end do
                              end do
                              dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) / real(2**count(native_nodal), WP)
                           end do
                        end do
                     end do
                  end do
                  call amrex_mfiter_destroy(mfi)

                  ! Process Z component
                  icomp = icomp + 1
                  call amrex_mfiter_build(mfi, combined(lev), tiling=this%amr%default_tiling)
                  do while (mfi%next())
                     bx = mfi%tilebox()
                     src => my_vct%ptrz%mf(lev)%dataptr(mfi)
                     dst => combined(lev)%dataptr(mfi)
                     native_nodal = my_vct%ptrz%nodal
                     do kk = bx%lo(3), bx%hi(3)
                        do jj = bx%lo(2), bx%hi(2)
                           do ii = bx%lo(1), bx%hi(1)
                              dst(ii,jj,kk,icomp) = 0.0_WP
                              do dk = 0, merge(1,0,native_nodal(3))
                                 do dj = 0, merge(1,0,native_nodal(2))
                                    do di = 0, merge(1,0,native_nodal(1))
                                       dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) + src(ii+di, jj+dj, kk+dk, my_vct%compz)
                                    end do
                                 end do
                              end do
                              dst(ii,jj,kk,icomp) = dst(ii,jj,kk,icomp) / real(2**count(native_nodal), WP)
                           end do
                        end do
                     end do
                  end do
                  call amrex_mfiter_destroy(mfi)

                  my_vct => my_vct%next
               end do
            end if
         end do
         filename = 'amrviz/'//trim(this%name)//'/nga2.'//trim(suffix)//'.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime

         ! Prepare pointers for HDF5 writer
         allocate(mf_ptrs(0:this%amr%clvl()))
         allocate(geom_ptrs(0:this%amr%clvl()))
         allocate(level_steps(0:this%amr%clvl()))
         allocate(ref_ratios(0:max(this%amr%clvl()-1,0)))

         do lev = 0, this%amr%clvl()
            mf_ptrs(lev) = combined(lev)%p
            geom_ptrs(lev) = this%amr%geom(lev)%p
            level_steps(lev) = this%ntime
         end do
         do lev = 0, this%amr%clvl()-1
            ref_ratios(lev) = this%amr%rref(lev)
         end do

         ! Write HDF5 file for this centering group
         call amrplotfile_write_hdf5(trim(filename)//c_null_char, nlev, mf_ptrs, &
            varname_ptrs, ncomp, geom_ptrs, real(time,c_double), level_steps, &
            ref_ratios, c_null_char)

         ! Cleanup for this group
         do lev = 0, this%amr%clvl()
            call amrex_multifab_destroy(combined(lev))
         end do
         deallocate(combined)
         deallocate(mf_ptrs, geom_ptrs, level_steps, ref_ratios)
         deallocate(varnames, varname_ptrs)
      end do

   contains

      !> Return filename suffix based on centering
      function centering_suffix(nodal) result(suf)
         logical, intent(in) :: nodal(3)
         character(len=8) :: suf
         if (.not.any(nodal)) then
            suf = 'cell'     ! Cell-centered
         else if (all(nodal)) then
            suf = 'node'     ! Node-centered
         else if (nodal(1) .and. .not.nodal(2) .and. .not.nodal(3)) then
            suf = 'xface'    ! X-face
         else if (.not.nodal(1) .and. nodal(2) .and. .not.nodal(3)) then
            suf = 'yface'    ! Y-face
         else if (.not.nodal(1) .and. .not.nodal(2) .and. nodal(3)) then
            suf = 'zface'    ! Z-face
         else if (nodal(1) .and. nodal(2) .and. .not.nodal(3)) then
            suf = 'xyedge'   ! XY-edge
         else if (nodal(1) .and. .not.nodal(2) .and. nodal(3)) then
            suf = 'xzedge'   ! XZ-edge
         else
            suf = 'yzedge'   ! YZ-edge
         end if
      end function centering_suffix

   end subroutine write


   !> Finalize: clean up registered field lists
   subroutine finalize(this)
      implicit none
      class(amrviz), intent(inout) :: this
      type(scl), pointer :: cur_scl, next_scl
      type(vct), pointer :: cur_vct, next_vct

      ! Clean up scalar list
      cur_scl => this%first_scl
      do while (associated(cur_scl))
         next_scl => cur_scl%next
         deallocate(cur_scl)
         cur_scl => next_scl
      end do
      nullify(this%first_scl)

      ! Clean up vector list
      cur_vct => this%first_vct
      do while (associated(cur_vct))
         next_vct => cur_vct%next
         deallocate(cur_vct)
         cur_vct => next_vct
      end do
      nullify(this%first_vct)

      ! Reset state
      nullify(this%amr)
      this%ntime = 0
      if (allocated(this%time)) deallocate(this%time)
   end subroutine finalize

end module amrviz_class
