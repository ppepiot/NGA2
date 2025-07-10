!> Chemical state
module chem_state_class
    use precision,      only: WP
    use chem_sys_class, only: chem_sys
    implicit none
    private
   
   ! Expose type/constructor/methods
   public :: chem_state

   !> Chemical statete object definition
   type :: chem_state

      ! This is our chemical system
      class(chem_sys), pointer :: sys   !< This is the chemical system the solver is build for

      ! Thermochemical quantities
      real(WP)          :: p            !< Pressure
      real(WP)          :: T            !< Temperature
      real(WP)          :: h            !< Enthalpy
      real(WP), pointer :: zd(:)        !< Moles of determined species (nsd)
      real(WP), pointer :: zu(:)        !< Moles of undetermined species (nsu)
      real(WP), pointer :: lam(:)       !< Lagrange multipliers (nrc)
      real(WP), pointer :: cr(:)        !< Reduced constraint vector
      real(WP), pointer :: Q(:)         !< Consistency vector Q' * Xu = 1 (nsu)
      real(WP), pointer :: gu(:)        !< Gibbs functions of undetermined species (nsu)

      ! Numerical parameters
      integer :: temp_its               !< Number of temperature iterations performed
      integer :: time_steps             !< Number of pseudo time steps
      integer :: newt_calls             !< Number of Newton solves attempted
      integer :: newt_its               !< Number of Newton iterations performed
      integer :: nyeval                 !< Number of y evaluations

   contains
      procedure :: initialize
   end type chem_state

   contains


      !> Chemical state initializer
      subroutine initialize(this,sys)
         use messager, only: die
         implicit none
         class(chem_state), intent(inout) :: this
         class(chem_sys), target, intent(in) :: sys

         ! Point to chem_sys object
         this%sys=>sys
      end subroutine initialize


end module chem_state_class