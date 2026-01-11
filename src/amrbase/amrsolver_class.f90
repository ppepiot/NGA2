!> Abstract base class for AMR solvers
!> Provides context for coupled BC callbacks and regrid handling
module amrsolver_class
   use iso_c_binding, only: c_ptr,c_null_ptr,c_loc,c_funloc,c_f_pointer
   use precision,     only: WP
   use string,        only: str_medium
   implicit none
   private

   public :: amrsolver

   !> Abstract solver base - extended by concrete solvers (amrscalar, mpcomp, etc.)
   type, abstract :: amrsolver
      ! Name for logging/debugging
      character(len=str_medium) :: name='UNNAMED_SOLVER'
   contains
      ! Deferred procedures that concrete solvers must implement
      procedure(fillbc_interface), deferred :: fillbc      !< Apply physical BCs
      procedure(regrid_interface), deferred :: on_regrid   !< Handle regrid event
   end type amrsolver

   !> Interface for fillbc - receives solver context, can access all fields
   abstract interface
      subroutine fillbc_interface(this,mf_ptr,geom_ptr,time,scomp,ncomp)
         import :: amrsolver,c_ptr,WP
         class(amrsolver), intent(inout) :: this
         type(c_ptr), value :: mf_ptr     !< MultiFab pointer
         type(c_ptr), value :: geom_ptr   !< Geometry pointer
         real(WP), value :: time
         integer, value :: scomp,ncomp
      end subroutine fillbc_interface

      subroutine regrid_interface(this,lvl,time)
         import :: amrsolver,WP
         class(amrsolver), intent(inout) :: this
         integer, value :: lvl
         real(WP), value :: time
      end subroutine regrid_interface
   end interface

end module amrsolver_class
