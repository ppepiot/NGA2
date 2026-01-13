!> AMR I/O class: Checkpoint and restart using AMReX VisMF
!> Provides directory-based checkpoint format with Header + Level_N/ structure
module amrio_class
   use precision,      only: WP
   use string,         only: str_medium
   use iso_c_binding,  only: c_ptr,c_char,c_int,c_null_char
   use amrgrid_class,  only: amrgrid
   use amrdata_class,  only: amrdata
   use amrex_interface, only: amrmfab_vismf_write,amrmfab_vismf_read, &
      amrcheckpoint_prebuild_dirs,amrcheckpoint_mfab_prefix
   implicit none
   private

   public :: amrio

   !> AMR I/O handler for checkpoints
   type :: amrio
      class(amrgrid), pointer :: amr => null()  !< Pointer to AMR grid
   contains
      procedure :: initialize                   !< Initialize with grid
      procedure :: write_checkpoint             !< Write checkpoint to disk
      procedure :: read_checkpoint              !< Read checkpoint from disk
   end type amrio

contains

   !> Initialize I/O handler with AMR grid
   subroutine initialize(this, amr)
      implicit none
      class(amrio), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      this%amr => amr
   end subroutine initialize


   !> Write checkpoint to disk
   !> Creates directory structure: dirname/Header + dirname/Level_N/
   subroutine write_checkpoint(this, dirname, time, step, data, dataname)
      use messager, only: log
      use parallel, only: MPI_REAL_WP
      use mpi_f08,  only: MPI_BCAST,MPI_INTEGER
      implicit none
      class(amrio), intent(inout) :: this
      character(len=*), intent(in) :: dirname      !< Checkpoint directory name
      real(WP), intent(in) :: time                 !< Simulation time
      integer, intent(in) :: step                  !< Step number
      type(amrdata), intent(in) :: data            !< Data to save
      character(len=*), intent(in) :: dataname     !< Name for the data (e.g., 'SC')

      integer :: lev, ierr, iunit
      character(len=256) :: mfab_path
      character(len=str_medium) :: header_file

      ! Create directory hierarchy
      call amrcheckpoint_prebuild_dirs(trim(dirname)//c_null_char, &
         'Level_'//c_null_char, this%amr%nlevels)

      ! Write Header file (root only)
      if (this%amr%amRoot) then
         header_file = trim(dirname)//'/Header'
         open(newunit=iunit, file=trim(header_file), status='replace', action='write')
         write(iunit,'(A)') 'NGA2 AMR Checkpoint'
         write(iunit,'(I0)') this%amr%nlevels - 1  ! finest_level (0-indexed)
         write(iunit,'(I0)') step
         write(iunit,'(ES23.16)') time
         ! Write BoxArray info per level
         do lev = 0, this%amr%maxlvl
            write(iunit,'(I0,A,I0)') lev, ' ', this%amr%rref(min(lev, this%amr%maxlvl-1))
         end do
         close(iunit)
         call log('Wrote checkpoint Header: '//trim(header_file))
      end if

      ! Write MultiFab data for each level
      do lev = 0, this%amr%maxlvl
         ! Get the path prefix for this level's MultiFab
         call amrcheckpoint_mfab_prefix(mfab_path, len(mfab_path), lev, &
            trim(dirname)//c_null_char, 'Level_'//c_null_char, &
            trim(dataname)//c_null_char)
         ! Write the MultiFab
         call amrmfab_vismf_write(data%mf(lev)%p, trim(mfab_path)//c_null_char)
      end do

      if (this%amr%amRoot) call log('Wrote checkpoint: '//trim(dirname))
   end subroutine write_checkpoint


   !> Read checkpoint from disk
   !> Reads Header and MultiFab data, returns time and step
   subroutine read_checkpoint(this, dirname, time, step, data, dataname)
      use messager, only: log,die
      implicit none
      class(amrio), intent(inout) :: this
      character(len=*), intent(in) :: dirname      !< Checkpoint directory name
      real(WP), intent(out) :: time                !< Simulation time (output)
      integer, intent(out) :: step                 !< Step number (output)
      type(amrdata), intent(inout) :: data         !< Data to read into
      character(len=*), intent(in) :: dataname     !< Name for the data (e.g., 'SC')

      integer :: lev, ierr, iunit, finest_level
      character(len=256) :: mfab_path, line
      character(len=str_medium) :: header_file

      ! Read Header file (root only, then broadcast)
      header_file = trim(dirname)//'/Header'
      if (this%amr%amRoot) then
         open(newunit=iunit, file=trim(header_file), status='old', action='read', iostat=ierr)
         if (ierr /= 0) call die('Cannot open checkpoint Header: '//trim(header_file))
         read(iunit,'(A)') line  ! Title
         read(iunit,*) finest_level
         read(iunit,*) step
         read(iunit,*) time
         close(iunit)
         call log('Read checkpoint Header: '//trim(header_file))
      end if

      ! Broadcast to all ranks
      block
         use mpi_f08, only: MPI_BCAST,MPI_INTEGER
         use parallel, only: MPI_REAL_WP
         call MPI_BCAST(step, 1, MPI_INTEGER, 0, this%amr%comm, ierr)
         call MPI_BCAST(time, 1, MPI_REAL_WP, 0, this%amr%comm, ierr)
      end block

      ! Read MultiFab data for each level
      ! NOTE: The amrdata must already be allocated with correct BoxArray/DistroMap
      do lev = 0, this%amr%maxlvl
         call amrcheckpoint_mfab_prefix(mfab_path, len(mfab_path), lev, &
            trim(dirname)//c_null_char, 'Level_'//c_null_char, &
            trim(dataname)//c_null_char)
         call amrmfab_vismf_read(data%mf(lev)%p, trim(mfab_path)//c_null_char)
      end do

      if (this%amr%amRoot) call log('Read checkpoint: '//trim(dirname))
   end subroutine read_checkpoint

end module amrio_class
