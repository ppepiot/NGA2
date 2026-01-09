!> AMR config object is defined based on amrcore AMREX object
!> Amrconfig differs quite a bit from other configs
module amrgrid_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_char,c_int,c_funptr,c_funloc
   use precision,        only: WP
   use string,           only: str_medium
   use amrex_amr_module, only: amrex_geometry
   use mpi_f08,          only: MPI_Comm
   use amrdata_class,    only: amrdata
   use amrflux_class,    only: amrflux
   implicit none
   private


   ! Expose type/constructor/methods
   public :: amrgrid

   ! Wrappers for polymorphism in registry arrays
   type :: data_ptr_wrapper
      class(amrdata), pointer :: p=>null()
   end type data_ptr_wrapper
   type :: flux_ptr_wrapper
      class(amrflux), pointer :: p=>null()
   end type flux_ptr_wrapper

   !> Abstract interface for user-provided tagging callback
   abstract interface
      subroutine tagging_callback(lvl, tags, time)
         use iso_c_binding, only: c_ptr
         use precision, only: WP
         implicit none
         integer, intent(in) :: lvl
         type(c_ptr), intent(in) :: tags  !< amrex_tagboxarray C pointer
         real(WP), intent(in) :: time
      end subroutine tagging_callback
   end interface

   !> Amrgrid object definition based on AMReX's amrcore
   type :: amrgrid
      ! Name of amrgrid
      character(len=str_medium) :: name='UNNAMED_AMRGRID'
      ! Pointer to AMReX's amrcore object
      type(c_ptr) :: amrcore=c_null_ptr
      ! Coordinate system
      integer :: coordsys=0
      ! Domain periodicity
      logical :: xper,yper,zper
      ! Domain extent
      real(WP) :: xlo,ylo,zlo
      real(WP) :: xhi,yhi,zhi
      ! Level 0 mesh dimensions
      integer :: nx,ny,nz
      ! Number of refinement levels
      integer :: nlvl
      ! Max grid size
      integer :: nmax=32
      ! Blocking factor
      integer :: nbloc=8
      ! Refinement ratio
      integer, dimension(:), allocatable :: rref
      ! Geometry object at each level
      type(amrex_geometry), dimension(:), allocatable :: geom
      ! Shortcut for domain volume
      real(WP) :: vol
      ! Shortcut to cell size per level
      real(WP), dimension(:), allocatable :: dx,dy,dz
      real(WP) :: min_meshsize
      ! Parallel info
      type(MPI_Comm) :: comm            !< Communicator for our group
      integer        :: nproc           !< Number of processors
      integer        :: rank            !< Processor grid rank
      logical        :: amRoot          !< Am I the root?
      ! Monitoring info - call get_info() to update
      integer  :: nlevels=-1            !< Current total number of levels
      integer  :: nboxes =-1            !< Current total number of boxes
      real(WP) :: ncells =-1.0_WP       !< Current total number of cells (real!)
      real(WP) :: compression=-1.0_WP   !< Current compression ratio (ncells/total cells with uniform mesh)
      ! Registries
      type(data_ptr_wrapper), dimension(:), allocatable :: data
      type(flux_ptr_wrapper), dimension(:), allocatable :: flux
      ! User-provided tagging callback
      procedure(tagging_callback), pointer, nopass :: tagger => null()
   contains
      procedure :: initialize                !< Initialization of amrgrid object
      procedure :: finalize                  !< Finalization of amrgrid object
      procedure :: initialize_grid           !< Initialize data on armgrid according to registered function
      procedure :: regrid                    !< Perform regriding operation on level baselvl
      procedure :: get_info                  !< Calculate various information on our amrgrid object
      procedure :: print                     !< Print out grid info
      ! Generic register/unregister for both data and flux
      generic   :: register  =>register_data,register_flux
      generic   :: unregister=>unregister_data,unregister_flux
      procedure, private :: register_data
      procedure, private :: unregister_data
      procedure, private :: register_flux
      procedure, private :: unregister_flux
      procedure :: register_internal_callbacks !< Register static callbacks

      procedure, private :: get_boxarray     !< Obtain box array at a given level
      procedure, private :: get_distromap    !< Obtain distromap at a given level

      procedure :: mfiter_build              !< Build mfiter at a given level
      procedure :: mfiter_destroy            !< Destroy mfiter
      procedure :: clvl                      !< Return current finest level
      procedure :: average_down              !< Average down a given multifab throughout all levels
      procedure :: average_downto            !< Average down a given multifab to level lvl
      procedure :: fill                      !< Fill ghost cells and coarse-fine boundaries
      procedure :: mfab_build                !< Build multifab at a given level
      procedure :: mfab_destroy              !< Destroy multifab
      procedure :: set_tagging               !< Set user-provided tagging callback
   end type amrgrid

   ! Singleton for static callback routing
   type(amrgrid),  pointer :: current_amrgrid=>null()
   class(amrdata), pointer :: current_amrdata=>null()

   ! Instance counter for automated AMReX lifecycle management
   integer :: amr_instance_count=0

contains


   !> Initialization of an amrgrid object
   subroutine initialize(this,name)
      use messager, only: die
      implicit none
      class(amrgrid), target, intent(inout) :: this
      character(len=*), intent(in), optional :: name
      ! First of all, confirm that nga2 and amrex reals are compatible
      check_real: block
         use messager,         only: die
         use precision,        only: WP
         use amrex_amr_module, only: amrex_real
         if (amrex_real.ne.WP) call die('Incompatible real type between nga2 and amrex!')
      end block check_real
      ! Check if AMReX has been initialized - if not, initialize it here
      check_if_init: block
         use amrex_amr_module, only: amrex_initialized,amrex_init
         use parallel, only: comm
         ! Increment instance counter
         amr_instance_count=amr_instance_count+1
         ! Initialize AMReX if it's the first instance or not yet initialized
         if (.not.amrex_initialized()) call amrex_init(comm=comm%MPI_VAL,arg_parmparse=.false.)
      end block check_if_init
      ! Set parameters for amrcore and geometry
      set_params: block
         use amrex_amr_module, only: amrex_parmparse,amrex_parmparse_build,amrex_parmparse,amrex_parmparse_destroy
         type(amrex_parmparse) :: pp
         integer, dimension(3) :: per
         call amrex_parmparse_build(pp,'amr')
         call pp%addarr('n_cell'         ,[this%nx,this%ny,this%nz])
         if (this%nlvl.lt.0) call die('[amrgrid initialize] nlvl must be >= 0')
         call pp%add   ('max_level'      ,this%nlvl)
         call pp%add   ('blocking_factor',this%nbloc)
         call pp%add   ('max_grid_size'  ,this%nmax)
         if (.not.allocated(this%rref)) this%rref=[2]
         call pp%addarr('ref_ratio'      ,this%rref)
         call amrex_parmparse_destroy(pp)
         call amrex_parmparse_build(pp,'geometry')
         call pp%add   ('coord_sys'      ,this%coordsys)
         per=0
         if (this%xper) per(1)=1
         if (this%yper) per(2)=1
         if (this%zper) per(3)=1
         call pp%addarr('is_periodic',per)
         call pp%addarr('prob_lo'    ,[this%xlo,this%ylo,this%zlo])
         call pp%addarr('prob_hi'    ,[this%xhi,this%yhi,this%zhi])
         call amrex_parmparse_destroy(pp)
      end block set_params
      ! Create an amrcore object
      create_amrcore_obj: block
         interface
            subroutine amrex_fi_new_amrcore(core) bind(c)
               import :: c_ptr
               implicit none
               type(c_ptr) :: core
            end subroutine amrex_fi_new_amrcore
         end interface
         call amrex_fi_new_amrcore(this%amrcore)
         if (present(name)) this%name=trim(adjustl(name))
         ! Setup Singleton for static callbacks
         current_amrgrid => this
      end block create_amrcore_obj
      ! Get back geometry objects
      store_geometries: block
         use amrex_amr_module, only: amrex_geometry_init_data
         interface
            subroutine amrex_fi_get_geometry(geom,lvl,core) bind(c)
               import :: c_ptr,c_int
               implicit none
               type(c_ptr), intent(out) :: geom
               integer(c_int), value :: lvl
               type(c_ptr), value :: core
            end subroutine amrex_fi_get_geometry
         end interface
         integer :: n
         allocate(this%geom(0:this%nlvl))
         do n=0,this%nlvl
            call amrex_fi_get_geometry(this%geom(n)%p,n,this%amrcore)
            call amrex_geometry_init_data(this%geom(n))
         end do
      end block store_geometries
      ! Store effective refinement ratio
      store_ref_ratio: block
         interface
            subroutine amrex_fi_get_ref_ratio(ref_ratio,core) bind(c)
               import :: c_ptr
               implicit none
               integer, dimension(*), intent(inout) :: ref_ratio
               type(c_ptr), value :: core
            end subroutine amrex_fi_get_ref_ratio
         end interface
         if (allocated(this%rref)) deallocate(this%rref); allocate(this%rref(0:this%nlvl-1))
         call amrex_fi_get_ref_ratio(this%rref,this%amrcore)
      end block store_ref_ratio
      ! Store parallel info
      store_parallel_info: block
         use parallel, only: comm,nproc,rank,amRoot
         this%comm=comm
         this%nproc=nproc
         this%rank=rank
         this%amRoot=amRoot
      end block store_parallel_info
      ! Compute shortcuts
      compute_shortcuts: block
         integer :: lvl
         ! Mesh size at each level
         allocate(this%dx(0:this%nlvl),this%dy(0:this%nlvl),this%dz(0:this%nlvl))
         do lvl=0,this%nlvl
            this%dx(lvl)=this%geom(lvl)%dx(1)
            this%dy(lvl)=this%geom(lvl)%dx(2)
            this%dz(lvl)=this%geom(lvl)%dx(3)
         end do
         ! Total domain volume
         this%vol=(this%xhi-this%xlo)*(this%yhi-this%ylo)*(this%zhi-this%zlo)
         ! Smallest mesh size
         this%min_meshsize=min(this%dx(this%nlvl),this%dy(this%nlvl),this%dz(this%nlvl))
      end block compute_shortcuts
   end subroutine initialize


   !> Finalization of amrgrid object
   impure elemental subroutine finalize(this)
      use mpi_f08,       only: MPI_COMM_NULL
      use iso_c_binding, only: c_associated
      implicit none
      class(amrgrid), intent(inout) :: this
      interface
         subroutine amrex_fi_delete_amrcore(core) bind(c)
            import :: c_ptr
            implicit none
            type(c_ptr), value :: core
         end subroutine amrex_fi_delete_amrcore
      end interface
      ! Delete amrcore object if it exists
      if (c_associated(this%amrcore)) then
         call amrex_fi_delete_amrcore(this%amrcore)
         this%amrcore=c_null_ptr
      end if
      ! Deallocate allocatable arrays
      if (allocated(this%rref)) deallocate(this%rref)
      if (allocated(this%geom)) deallocate(this%geom)
      if (allocated(this%dx))   deallocate(this%dx)
      if (allocated(this%dy))   deallocate(this%dy)
      if (allocated(this%dz))   deallocate(this%dz)
      if (allocated(this%data)) deallocate(this%data)
      if (allocated(this%flux)) deallocate(this%flux)
      ! Do not free comm as it was passed to us
      this%comm=MPI_COMM_NULL
      ! Handle automated AMReX finalization
      auto_finalize: block
         use amrex_amr_module, only: amrex_initialized,amrex_finalize
         amr_instance_count=amr_instance_count-1
         if (amr_instance_count.eq.0.and.amrex_initialized()) call amrex_finalize()
      end block auto_finalize
   end subroutine finalize


   !> Register an amrdata object to be managed by the amrgrid
   !> Optional: name, ncomp, ng to configure the data at registration
   subroutine register_data(this,data,name,ncomp,ng)
      use amrex_amr_module, only: amrex_bc_int_dir
      implicit none
      class(amrgrid), intent(inout) :: this
      class(amrdata), target, intent(inout) :: data
      character(len=*), intent(in), optional :: name
      integer, intent(in), optional :: ncomp,ng
      type(data_ptr_wrapper), dimension(:), allocatable :: temp
      integer :: n
      ! Configure data metadata
      if (present(name))  data%name=trim(adjustl(name))
      if (present(ncomp)) data%ncomp=ncomp
      if (present(ng))    data%ng=ng
      ! Allocate BC arrays - default to int_dir (interior/do nothing)
      ! User sets specific BCs after registration: data%lo_bc(1,:)=amrex_bc_foextrap
      if (allocated(data%lo_bc)) deallocate(data%lo_bc)
      if (allocated(data%hi_bc)) deallocate(data%hi_bc)
      allocate(data%lo_bc(3,data%ncomp),data%hi_bc(3,data%ncomp))
      data%lo_bc=amrex_bc_int_dir
      data%hi_bc=amrex_bc_int_dir
      ! Allocate/resize registry
      if (.not.allocated(this%data)) then
         allocate(this%data(1))
      else
         n=size(this%data)
         allocate(temp(n+1))
         temp(1:n)=this%data
         call move_alloc(temp,this%data)
      end if
      ! Add to registry
      n=size(this%data)
      this%data(n)%p=>data
      ! Allocate data levels to match the grid
      if (.not.allocated(data%mf)) allocate(data%mf(0:this%nlvl))
   end subroutine register_data

   !> Unregister an amrdata object from the registry
   subroutine unregister_data(this,data)
      implicit none
      class(amrgrid), intent(inout) :: this
      class(amrdata), target, intent(in) :: data
      type(data_ptr_wrapper), dimension(:), allocatable :: temp
      integer :: i,n,pos
      if (.not.allocated(this%data)) return
      n=size(this%data)
      pos=0
      ! Find position
      do i=1,n
         if (associated(this%data(i)%p,data)) then
            pos=i
            exit
         end if
      end do
      if (pos.gt.0) then
         if (n.eq.1) then
            deallocate(this%data)
         else
            allocate(temp(n-1))
            if (pos.gt.1) temp(1:pos-1)=this%data(1:pos-1)
            if (pos.lt.n) temp(pos:)   =this%data(pos+1:)
            call move_alloc(temp,this%data)
         end if
      end if
   end subroutine unregister_data

   !> Register an amrflux object to be managed by the amrgrid
   !> Optional: name, ncomp to configure the flux at registration
   subroutine register_flux(this,flux,name,ncomp)
      implicit none
      class(amrgrid), intent(inout) :: this
      class(amrflux), target, intent(inout) :: flux
      character(len=*), intent(in), optional :: name
      integer, intent(in), optional :: ncomp
      type(flux_ptr_wrapper), dimension(:), allocatable :: temp
      integer :: n
      ! Configure flux metadata
      if (present(name))  flux%name=trim(adjustl(name))
      if (present(ncomp)) flux%ncomp=ncomp
      ! Allocate/resize registry
      if (.not.allocated(this%flux)) then
         allocate(this%flux(1))
      else
         n=size(this%flux)
         allocate(temp(n+1))
         temp(1:n)=this%flux
         call move_alloc(temp,this%flux)
      end if
      ! Add to registry
      n=size(this%flux)
      this%flux(n)%p=>flux
      ! Allocate flux register array (levels 1:nlvl, no level 0)
      if (.not.allocated(flux%fr).and.this%nlvl.ge.1) allocate(flux%fr(1:this%nlvl))
   end subroutine register_flux

   !> Unregister an amrflux object from the registry
   subroutine unregister_flux(this,flux)
      implicit none
      class(amrgrid), intent(inout) :: this
      class(amrflux), target, intent(in) :: flux
      type(flux_ptr_wrapper), dimension(:), allocatable :: temp
      integer :: i,n,pos
      if (.not.allocated(this%flux)) return
      n=size(this%flux)
      pos=0
      ! Find position
      do i=1,n
         if (associated(this%flux(i)%p,flux)) then
            pos=i
            exit
         end if
      end do
      if (pos.gt.0) then
         if (n.eq.1) then
            deallocate(this%flux)
         else
            allocate(temp(n-1))
            if (pos.gt.1) temp(1:pos-1)=this%flux(1:pos-1)
            if (pos.lt.n) temp(pos:)   =this%flux(pos+1:)
            call move_alloc(temp,this%flux)
         end if
      end if
   end subroutine unregister_flux

   !> Initialize grid on an amrgrid object
   subroutine initialize_grid(this,time)
      implicit none
      class(amrgrid), target, intent(inout) :: this
      real(WP), intent(in) :: time
      interface
         subroutine amrex_fi_init_from_scratch(t,core) bind(c)
            import :: c_ptr,WP
            implicit none
            real(WP), value :: t
            type(c_ptr), value :: core
         end subroutine amrex_fi_init_from_scratch
      end interface
      ! Setup the singleton for static callbacks
      current_amrgrid=>this
      ! Register our internal static dispatchers with C++
      call this%register_internal_callbacks()
      ! Generate grid and allocate data
      call amrex_fi_init_from_scratch(time,this%amrcore)
      ! Generate info about grid
      call this%get_info()
   end subroutine initialize_grid


   !> Perform regriding operation on baselvl
   subroutine regrid(this,baselvl,time)
      implicit none
      class(amrgrid), target, intent(inout) :: this
      integer,  intent(in) :: baselvl
      real(WP), intent(in) :: time
      interface
         subroutine amrex_fi_regrid(blvl,t,core) bind(c)
            import :: c_ptr,c_int,WP
            implicit none
            integer(c_int), value :: blvl
            real(WP), value :: t
            type(c_ptr), value :: core
         end subroutine amrex_fi_regrid
      end interface
      ! Ensure Singleton matches this instance
      current_amrgrid=>this
      ! Regenerate grid and resize/transfer data
      call amrex_fi_regrid(baselvl,time,this%amrcore)
      ! Generate info about grid
      call this%get_info()
   end subroutine regrid

   !> Internal: Register the static dispatchers with AmrCore
   subroutine register_internal_callbacks(this)
      implicit none
      class(amrgrid), intent(inout) :: this
      interface
         subroutine amrex_fi_init_virtual_functions(maklvl_init,maklvl_crse,maklvl_remk,clrlvl,errest,core) bind(c)
            import :: c_ptr,c_funptr
            implicit none
            type(c_funptr), value :: maklvl_init,maklvl_crse,maklvl_remk,clrlvl,errest
            type(c_ptr), value :: core
         end subroutine amrex_fi_init_virtual_functions
      end interface
      call amrex_fi_init_virtual_functions(c_funloc(dispatch_mak_lvl_init), &
      &                                    c_funloc(dispatch_mak_lvl_crse), &
      &                                    c_funloc(dispatch_mak_lvl_remk), &
      &                                    c_funloc(dispatch_clr_lvl), &
      &                                    c_funloc(dispatch_err_est), &
      &                                    this%amrcore)
   end subroutine register_internal_callbacks

   !> Get info on amrgrid object
   subroutine get_info(this)
      use amrex_amr_module, only: amrex_boxarray,amrex_box
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_boxarray) :: ba
      type(amrex_box)      :: bx
      integer :: n,m,nb
      integer, dimension(3) :: lo,hi
      ! Update number of levels
      this%nlevels=this%clvl()+1
      ! Update number of boxes and cells
      this%nboxes=0
      this%ncells=0.0_WP
      do n=0,this%clvl()
         ! Get boxarray
         ba=this%get_boxarray(lvl=n)
         ! Increment boxes
         nb=int(ba%nboxes())
         this%nboxes=this%nboxes+nb
         ! Traverse boxes and count cells
         do m=1,nb
            ! Get box
            bx=ba%get_box(m-1)
            ! Increment cells
            this%ncells=this%ncells+real((bx%hi(1)-bx%lo(1)+1),WP)*&
            &                       real((bx%hi(2)-bx%lo(2)+1),WP)*&
            &                       real((bx%hi(3)-bx%lo(3)+1),WP)
         end do
      end do
      ! Update compression ratio
      this%compression=this%geom(this%clvl())%dx(1)*&
      &                this%geom(this%clvl())%dx(2)*&
      &                this%geom(this%clvl())%dx(3)*&
      &                this%ncells/((this%xhi-this%xlo)*(this%yhi-this%ylo)*(this%zhi-this%zlo))
   end subroutine get_info


   !> Print amrgrid object
   subroutine print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use amrex_amr_module, only: amrex_boxarray,amrex_box
      use parallel, only: amRoot
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_boxarray) :: ba
      type(amrex_box)      :: bx
      integer :: n,m,nb
      integer, dimension(3) :: lo,hi
      if (this%amRoot) then
         write(output_unit,'("AMR Cartesian grid [",a,"]")') trim(this%name)
         write(output_unit,'(" > amr level = ",i2)') this%clvl()
         write(output_unit,'(" > max level = ",i2)') this%nlvl
         if (allocated(this%rref)) write(output_unit,'(" > ref ratio = ",100(" ",i0))') this%rref
         write(output_unit,'(" >    extent = [",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]")') this%xlo,this%xhi,this%ylo,this%yhi,this%zlo,this%zhi
         write(output_unit,'(" >  periodic = ",l1,"x",l1,"x",l1)') this%xper,this%yper,this%zper
         ! Loop over levels
         do n=0,this%clvl()
            ! Get boxarray at that level
            ba=this%get_boxarray(lvl=n)
            write(output_unit,'(" >>> Level ",i2,"/",i2," with [dx =",es12.5,",dy =",es12.5,",dz =",es12.5,"]")') n,this%clvl(),this%geom(n)%dx(1),this%geom(n)%dx(2),this%geom(n)%dx(3)
            ! Loop over boxes
            nb=int(ba%nboxes())
            do m=1,nb
               bx=ba%get_box(m-1)
               write(output_unit,'(" >>> Box ",i0,"/",i0," from [",i3,",",i3,",",i3,"] to [",i3,",",i3,",",i3,"]")') m,nb,bx%lo,bx%hi
            end do
         end do
      end if
   end subroutine print


   !> Obtain box array at a level
   function get_boxarray(this,lvl) result(ba)
      use amrex_amr_module, only: amrex_boxarray
      implicit none
      class(amrgrid), intent(inout) :: this
      integer, intent(in)  :: lvl
      type(amrex_boxarray) :: ba
      interface
         subroutine amrex_fi_get_boxarray(barray,lev,core) bind(c)
            import :: c_ptr,c_int
            implicit none
            type(c_ptr), intent(out) :: barray
            integer(c_int), value :: lev
            type(c_ptr), value :: core
         end subroutine amrex_fi_get_boxarray
      end interface
      call amrex_fi_get_boxarray(ba%p,lvl,this%amrcore)
   end function get_boxarray


   !> Obtain distromap at a level
   function get_distromap(this,lvl) result(dm)
      use amrex_amr_module, only: amrex_distromap
      implicit none
      class(amrgrid), intent(inout) :: this
      integer, intent(in)  :: lvl
      type(amrex_distromap) :: dm
      interface
         subroutine amrex_fi_get_distromap(dmap,lev,core) bind(c)
            import :: c_ptr,c_int
            implicit none
            type(c_ptr), intent(out) :: dmap
            integer(c_int), value :: lev
            type(c_ptr), value :: core
         end subroutine amrex_fi_get_distromap
      end interface
      call amrex_fi_get_distromap(dm%p,lvl,this%amrcore)
   end function get_distromap


   !> Build mfiter at a level
   subroutine mfiter_build(this,lvl,mfi)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_mfiter_build
      implicit none
      class(amrgrid), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(out) :: mfi
      type(amrex_boxarray) :: ba
      type(amrex_distromap) :: dm
      ba=this%get_boxarray (lvl)
      dm=this%get_distromap(lvl)
      call amrex_mfiter_build(mfi,ba,dm)
   end subroutine mfiter_build


   !> Destroy mfiter
   subroutine mfiter_destroy(this,mfi)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_destroy
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_mfiter), intent(inout) :: mfi
      call amrex_mfiter_destroy(mfi)
   end subroutine mfiter_destroy

   ! --- STATIC DISPATCHERS (bind(c)) ---
   ! These receive calls from C++ and route them to current_amrgrid%data
   subroutine dispatch_mak_lvl_init(lvl, time, ba, dm) bind(c)
      use iso_c_binding, only: c_ptr, c_int, c_double
      use amrex_amr_module, only: amrex_boxarray, amrex_distromap
      implicit none
      integer(c_int), value, intent(in) :: lvl
      real(c_double), value, intent(in) :: time
      type(c_ptr), value, intent(in) :: ba, dm ! Raw pointers
      type(amrex_boxarray) :: ba_obj
      type(amrex_distromap) :: dm_obj
      integer :: i
      ! Wrap pointers
      ba_obj%p=ba
      dm_obj%p=dm
      ! Route to data
      if (associated(current_amrgrid)) then
         if (allocated(current_amrgrid%data)) then
            do i=1,size(current_amrgrid%data)
               call current_amrgrid%data(i)%p%define(int(lvl),ba_obj,dm_obj)
            end do
         end if
         ! Build flux registers for fine levels
         if (int(lvl).ge.1.and.allocated(current_amrgrid%flux)) then
            do i=1,size(current_amrgrid%flux)
               call current_amrgrid%flux(i)%p%build_level(int(lvl),ba_obj,dm_obj,&
               &   current_amrgrid%rref(int(lvl)-1))
            end do
         end if
      end if
   end subroutine dispatch_mak_lvl_init

   subroutine dispatch_mak_lvl_crse(lvl,time,ba,dm) bind(c)
      use iso_c_binding, only: c_ptr,c_int,c_double
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      integer(c_int), value, intent(in) :: lvl
      real(c_double), value, intent(in) :: time
      type(c_ptr), value, intent(in) :: ba,dm
      type(amrex_boxarray) :: ba_obj
      type(amrex_distromap) :: dm_obj
      integer :: i
      ba_obj%p=ba
      dm_obj%p=dm
      if (associated(current_amrgrid)) then
         ! Build data
         if (allocated(current_amrgrid%data)) then
            do i=1,size(current_amrgrid%data)
               call current_amrgrid%data(i)%p%define(int(lvl),ba_obj,dm_obj)
            end do
         end if
         ! Build flux registers for fine levels
         if (int(lvl).ge.1.and.allocated(current_amrgrid%flux)) then
            do i=1,size(current_amrgrid%flux)
               call current_amrgrid%flux(i)%p%build_level(int(lvl),ba_obj,dm_obj,&
               &   current_amrgrid%rref(int(lvl)-1))
            end do
         end if
      end if
   end subroutine dispatch_mak_lvl_crse

   subroutine dispatch_mak_lvl_remk(lvl,time,ba,dm) bind(c)
      use iso_c_binding, only: c_ptr,c_int,c_double
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      integer(c_int), value, intent(in) :: lvl
      real(c_double), value, intent(in) :: time
      type(c_ptr), value, intent(in) :: ba,dm
      type(amrex_boxarray) :: ba_obj
      type(amrex_distromap) :: dm_obj
      integer :: i
      ba_obj%p=ba
      dm_obj%p=dm
      if (associated(current_amrgrid)) then
         ! Build data
         if (allocated(current_amrgrid%data)) then
            do i=1,size(current_amrgrid%data)
               call current_amrgrid%data(i)%p%define(int(lvl),ba_obj,dm_obj)
            end do
         end if
         ! Build flux registers for fine levels
         if (int(lvl).ge.1.and.allocated(current_amrgrid%flux)) then
            do i=1,size(current_amrgrid%flux)
               call current_amrgrid%flux(i)%p%build_level(int(lvl),ba_obj,dm_obj,&
               &   current_amrgrid%rref(int(lvl)-1))
            end do
         end if
      end if
   end subroutine dispatch_mak_lvl_remk

   subroutine dispatch_clr_lvl(lvl) bind(c)
      use iso_c_binding, only: c_int
      implicit none
      integer(c_int), value, intent(in) :: lvl
      integer :: i
      if (associated(current_amrgrid)) then
         ! Clear data at this level
         if (allocated(current_amrgrid%data)) then
            do i=1,size(current_amrgrid%data)
               call current_amrgrid%data(i)%p%clear_level(int(lvl))
            end do
         end if
         ! Clear flux registers at this level
         if (int(lvl).ge.1.and.allocated(current_amrgrid%flux)) then
            do i=1,size(current_amrgrid%flux)
               call current_amrgrid%flux(i)%p%destroy_level(int(lvl))
            end do
         end if
      end if
   end subroutine dispatch_clr_lvl


   !> Error estimation callback - calls user's tagger if set
   subroutine dispatch_err_est(lvl,tags,time,tagval,clrval) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int, c_double
      implicit none
      integer(c_int), value, intent(in) :: lvl
      type(c_ptr), value, intent(in) :: tags
      real(c_double), value, intent(in) :: time
      character(kind=c_char), value, intent(in) :: tagval,clrval
      ! Call user's tagging callback if set
      if (associated(current_amrgrid%tagger)) then
         call current_amrgrid%tagger(int(lvl),tags,real(time,WP))
      end if
   end subroutine dispatch_err_est


   !> Set the user-provided tagging callback
   subroutine set_tagging(this,tagger)
      implicit none
      class(amrgrid), intent(inout) :: this
      procedure(tagging_callback), optional :: tagger
      if (present(tagger)) then
         this%tagger=>tagger
      else
         this%tagger=>null()
      end if
   end subroutine set_tagging


   !> Return current finest level
   function clvl(this) result(cl)
      use amrex_amr_module, only: amrex_boxarray
      implicit none
      class(amrgrid), intent(inout) :: this
      integer :: cl
      interface
         integer(c_int) function amrex_fi_get_finest_level(core) bind(c)
            import :: c_int,c_ptr
            implicit none
            type(c_ptr), value :: core
         end function amrex_fi_get_finest_level
      end interface
      cl=amrex_fi_get_finest_level(this%amrcore)
   end function clvl


   !> Average down entire multifab array
   subroutine average_down(this,mfaba)
      use amrex_amr_module, only: amrex_multifab,amrex_average_down
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_multifab), dimension(0:) :: mfaba
      integer :: n
      do n=this%clvl()-1,0,-1
         call amrex_average_down(mfaba(n+1),mfaba(n),this%geom(n+1),this%geom(n),1,mfaba(0)%nc,this%rref(n))
      end do
   end subroutine average_down


   !> Average entire multifab array down to level lvl
   subroutine average_downto(this,mfaba,lvl)
      use amrex_amr_module, only: amrex_multifab,amrex_average_down
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_multifab), dimension(0:) :: mfaba
      integer, intent(in) :: lvl
      call amrex_average_down(mfaba(lvl+1),mfaba(lvl),this%geom(lvl+1),this%geom(lvl),1,mfaba(0)%nc,this%rref(lvl))
   end subroutine average_downto

   !> Fill ghost cells and coarse-fine boundary data
   !> For lvl=0: fills ghost cells using boundary conditions
   !> For lvl>0: interpolates from coarser level and fills ghosts
   !> Optional: data_old for time interpolation (subcycling)
   subroutine fill(this,data,lvl,time,data_old,t_old,t_new)
      use iso_c_binding, only: c_ptr,c_funloc
      use amrex_amr_module, only: amrex_fillpatch
      implicit none
      class(amrgrid), intent(inout) :: this
      class(amrdata), target, intent(inout) :: data
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      class(amrdata), target, intent(in), optional :: data_old
      real(WP), intent(in), optional :: t_old,t_new
      real(WP) :: time_old,time_new
      ! Set up times (if no old data, use single time)
      if (present(data_old).and.present(t_old).and.present(t_new)) then
         time_old=t_old
         time_new=t_new
      else
         time_old=time-1.0e200_WP
         time_new=time
      end if
      ! Set singleton for callback
      current_amrdata=>data
      ! Call fillpatch
      if (lvl.eq.0) then
         ! Level 0: just fill ghost cells with BC
         call amrex_fillpatch(data%mf(0),time_old,data%mf(0),time_new,data%mf(0),&
         &   this%geom(0),fillbc,time,1,1,data%ncomp)
      else
         ! Fine levels: interpolate from coarse + fill ghosts
         if (present(data_old)) then
            call amrex_fillpatch(data%mf(lvl),&
            &   time_old,data_old%mf(lvl-1),time_new,data%mf(lvl-1),this%geom(lvl-1),fillbc,&
            &   time_old,data_old%mf(lvl  ),time_new,data%mf(lvl  ),this%geom(lvl  ),fillbc,&
            &   time,1,1,data%ncomp,this%rref(lvl-1),data%interp,data%lo_bc,data%hi_bc)
         else
            call amrex_fillpatch(data%mf(lvl),&
            &   time_old,data%mf(lvl-1),time_new,data%mf(lvl-1),this%geom(lvl-1),fillbc,&
            &   time_old,data%mf(lvl  ),time_new,data%mf(lvl  ),this%geom(lvl  ),fillbc,&
            &   time,1,1,data%ncomp,this%rref(lvl-1),data%interp,data%lo_bc,data%hi_bc)
         end if
      end if
   contains
      subroutine fillbc(pmf,scomp,ncomp,t,pgeom) bind(c)
         use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build
         use iso_c_binding, only: c_ptr,c_int
         type(c_ptr), value :: pmf,pgeom
         integer(c_int), value :: scomp,ncomp
         real(WP), value :: t
         type(amrex_geometry) :: geom
         type(amrex_multifab) :: mf
         type(amrex_mfiter) :: mfi
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
         integer, dimension(4) :: plo,phi
         ! Skip if fully periodic
         if (this%xper.and.this%yper.and.this%zper) return
         ! Convert pointers
         geom=pgeom; mf=pmf
         ! Loop over boxes
         call amrex_mfiter_build(mfi,mf)
         do while(mfi%next())
            p=>mf%dataptr(mfi)
            ! Check if part of box is outside the domain
            if (.not.geom%domain%contains(p)) then
               plo=lbound(p); phi=ubound(p)
               call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,&
               &   geom%get_physical_location(plo),current_amrdata%lo_bc,current_amrdata%hi_bc)
            end if
         end do
      end subroutine fillbc
   end subroutine fill

   !> Build multifab at a level
   subroutine mfab_build(this,lvl,mfab,ncomp,nover,atface)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_multifab,amrex_multifab_build
      implicit none
      class(amrgrid), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: mfab
      integer, intent(in) :: ncomp, nover
      logical, intent(in), optional :: atface(3)
      type(amrex_boxarray) :: ba
      type(amrex_distromap) :: dm
      ba=this%get_boxarray(lvl)
      dm=this%get_distromap(lvl)
      call amrex_multifab_build(mfab,ba,dm,ncomp,nover,atface)
   end subroutine mfab_build

   !> Destroy multifab
   subroutine mfab_destroy(this,mfab)
      use amrex_amr_module, only: amrex_multifab,amrex_multifab_destroy
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mfab
      call amrex_multifab_destroy(mfab)
   end subroutine mfab_destroy


   !> Finalization of amrex
   subroutine finalize_amrex()
      use amrex_amr_module, only: amrex_finalize
      implicit none
      call amrex_finalize()
   end subroutine finalize_amrex


end module amrgrid_class
