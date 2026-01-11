!> Centralized C interfaces for all AMReX wrappers in NGA2
!> These wrap the C++ functions in amrex_interface.cpp
!> Provides: AmrCore lifecycle, dispatchers, queries, and FillPatch operations
module amrex_interface
   use iso_c_binding, only: c_ptr,c_funptr,c_double,c_int
   use precision !< Force dep.py to track this module
   implicit none
   private

   !=====================================================================
   ! AmrCore Lifecycle
   !=====================================================================
   public :: amrcore_create
   public :: amrcore_destroy
   public :: amrcore_set_owner

   !=====================================================================
   ! AmrCore Dispatcher Setters
   !=====================================================================
   public :: amrcore_set_on_init_dispatch
   public :: amrcore_set_on_coarse_dispatch
   public :: amrcore_set_on_remake_dispatch
   public :: amrcore_set_on_clear_dispatch
   public :: amrcore_set_on_tag_dispatch
   public :: amrcore_set_on_postregrid_dispatch

   !=====================================================================
   ! AmrCore Grid Operations
   !=====================================================================
   public :: amrcore_init_from_scratch
   public :: amrcore_regrid

   !=====================================================================
   ! AmrCore Queries
   !=====================================================================
   public :: amrcore_finest_level
   public :: amrcore_get_geometry
   public :: amrcore_get_ref_ratio
   public :: amrcore_get_boxarray
   public :: amrcore_get_distromap

   !=====================================================================
   ! FillPatch Operations
   !=====================================================================
   public :: amrmfab_fillpatch_single
   public :: amrmfab_fillpatch_two
   public :: amrmfab_fillcoarsepatch

   interface

      !------------------------------------------------------------------
      ! AmrCore Lifecycle
      !------------------------------------------------------------------

      !> Create an NGA2AmrCore object
      type(c_ptr) function amrcore_create() bind(c)
         import :: c_ptr
      end function amrcore_create

      !> Destroy an NGA2AmrCore object
      subroutine amrcore_destroy(core) bind(c)
         import :: c_ptr
         type(c_ptr), value :: core
      end subroutine amrcore_destroy

      !> Set owner pointer for callbacks
      subroutine amrcore_set_owner(core,owner) bind(c)
         import :: c_ptr
         type(c_ptr), value :: core,owner
      end subroutine amrcore_set_owner

      !------------------------------------------------------------------
      ! AmrCore Dispatcher Setters
      !------------------------------------------------------------------

      !> Set MakeNewLevelFromScratch dispatcher
      subroutine amrcore_set_on_init_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_init_dispatch

      !> Set MakeNewLevelFromCoarse dispatcher
      subroutine amrcore_set_on_coarse_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_coarse_dispatch

      !> Set RemakeLevel dispatcher
      subroutine amrcore_set_on_remake_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_remake_dispatch

      !> Set ClearLevel dispatcher
      subroutine amrcore_set_on_clear_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_clear_dispatch

      !> Set ErrorEst (tagging) dispatcher
      subroutine amrcore_set_on_tag_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_tag_dispatch

      !> Set post-regrid dispatcher
      subroutine amrcore_set_on_postregrid_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_postregrid_dispatch

      !------------------------------------------------------------------
      ! AmrCore Grid Operations
      !------------------------------------------------------------------

      !> Initialize grid from scratch
      subroutine amrcore_init_from_scratch(core,t) bind(c)
         import :: c_ptr,c_double
         type(c_ptr), value :: core
         real(c_double), value :: t
      end subroutine amrcore_init_from_scratch

      !> Perform regrid operation
      subroutine amrcore_regrid(core,blvl,t) bind(c)
         import :: c_ptr,c_int,c_double
         type(c_ptr), value :: core
         integer(c_int), value :: blvl
         real(c_double), value :: t
      end subroutine amrcore_regrid

      !------------------------------------------------------------------
      ! AmrCore Queries
      !------------------------------------------------------------------

      !> Get current finest level
      integer(c_int) function amrcore_finest_level(core) bind(c)
         import :: c_int,c_ptr
         type(c_ptr), value :: core
      end function amrcore_finest_level

      !> Get geometry at a level
      subroutine amrcore_get_geometry(geom,lvl,core) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), intent(out) :: geom
         integer(c_int), value :: lvl
         type(c_ptr), value :: core
      end subroutine amrcore_get_geometry

      !> Get refinement ratios
      subroutine amrcore_get_ref_ratio(ref_ratio,core) bind(c)
         import :: c_ptr
         integer, dimension(*), intent(inout) :: ref_ratio
         type(c_ptr), value :: core
      end subroutine amrcore_get_ref_ratio

      !> Get BoxArray at a level
      subroutine amrcore_get_boxarray(barray,lev,core) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), intent(out) :: barray
         integer(c_int), value :: lev
         type(c_ptr), value :: core
      end subroutine amrcore_get_boxarray

      !> Get DistributionMapping at a level
      subroutine amrcore_get_distromap(dmap,lev,core) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), intent(out) :: dmap
         integer(c_int), value :: lev
         type(c_ptr), value :: core
      end subroutine amrcore_get_distromap

      !------------------------------------------------------------------
      ! FillPatch Operations
      !------------------------------------------------------------------

      !> FillPatch for level 0 (single level, physical BCs only)
      subroutine amrmfab_fillpatch_single(mf,t_old,mf_old,t_new,mf_new, &
      &   geom,solver_ctx,bc_dispatch,time,scomp,dcomp,ncomp) bind(c)
         import :: c_ptr,c_funptr,c_double,c_int
         type(c_ptr), value :: mf,mf_old,mf_new,geom,solver_ctx
         type(c_funptr), value :: bc_dispatch
         real(c_double), value :: t_old,t_new,time
         integer(c_int), value :: scomp,dcomp,ncomp
      end subroutine amrmfab_fillpatch_single

      !> FillPatch for fine levels (two-level interpolation + BCs)
      subroutine amrmfab_fillpatch_two(mf,t_old_c,mf_old_c,t_new_c,mf_new_c,geom_c, &
      &   t_old_f,mf_old_f,t_new_f,mf_new_f,geom_f,solver_ctx,bc_dispatch, &
      &   time,scomp,dcomp,ncomp,ref_ratio,interp_type,lo_bc,hi_bc,nbc) bind(c)
         import :: c_ptr,c_funptr,c_double,c_int
         type(c_ptr), value :: mf,mf_old_c,mf_new_c,geom_c
         type(c_ptr), value :: mf_old_f,mf_new_f,geom_f,solver_ctx
         type(c_funptr), value :: bc_dispatch
         real(c_double), value :: t_old_c,t_new_c,t_old_f,t_new_f,time
         integer(c_int), value :: scomp,dcomp,ncomp,ref_ratio,interp_type,nbc
         integer(c_int), intent(in) :: lo_bc(*),hi_bc(*)
      end subroutine amrmfab_fillpatch_two

      !> FillCoarsePatch - fill fine level from coarse only (for new levels)
      subroutine amrmfab_fillcoarsepatch(mf_f,time,mf_c,geom_c,geom_f,solver_ctx, &
      &   bc_dispatch,scomp,dcomp,ncomp,ref_ratio,interp_type,lo_bc,hi_bc,nbc) bind(c)
         import :: c_ptr,c_funptr,c_double,c_int
         type(c_ptr), value :: mf_f,mf_c,geom_c,geom_f,solver_ctx
         type(c_funptr), value :: bc_dispatch
         real(c_double), value :: time
         integer(c_int), value :: scomp,dcomp,ncomp,ref_ratio,interp_type,nbc
         integer(c_int), intent(in) :: lo_bc(*),hi_bc(*)
      end subroutine amrmfab_fillcoarsepatch

   end interface

end module amrex_interface
