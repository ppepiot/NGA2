!> Abstract linear solver concept is defined here:
!> given a config, it provides parallel solvers
module linsol_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_medium
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: linsol
   
   
   !> linsol object definition
   type, abstract :: linsol
      
      ! A linear solver works for a config
      type(config), pointer :: cfg                                    !< Config for linsol
      
      ! An linsol has a name
      character(len=str_medium) :: name                               !< Name of solver
      
      ! Specific method to use
      integer  :: method                                              !< Solution method (we assume all solvers will have some sort of options)
      
      ! An linsol has a stencil size
      integer  :: nst                                                 !< Stencil size in 3D
      integer, dimension(:,:),   allocatable :: stc                   !< Stencil map in 3D (from stencil entry to i/j/k shift)
      integer, dimension(:,:,:), allocatable :: stmap                 !< Inverse stencil map in 3D (from i/j/k shift to stencil entry)
      
      ! An linsol stores the linear operator and rhs
      real(WP), dimension(:,:,:,:), allocatable :: opr                !< Linear operator
      real(WP), dimension(:,:,:),   allocatable :: rhs                !< RHS
      real(WP), dimension(:,:,:),   allocatable :: sol                !< Solution
      
      ! Is the solver setup?
      logical :: setup_done                                           !< Check whether the solver has been setup
      
      ! Convergence criteria
      integer  :: maxit                                               !< Maximum number of iterations allowed
      real(WP) :: rcvg                                                !< Desired relative convergence criterion
      real(WP) :: acvg                                                !< Desired absolue convergence criterion
      
      ! Current convergence info
      integer  :: it                                                  !< Current number of iterations
      real(WP) :: rerr                                                !< Current relative error
      real(WP) :: aerr                                                !< Current absolute error
      
   contains
      procedure(in_noarg_interface   ), deferred :: print_short             !< One-line printing of solver status
      procedure(in_noarg_interface   ), deferred :: print                   !< Long-form printing of solver status
      procedure(in_noarg_interface   ), deferred :: log                     !< Long-form logging of solver status
      procedure(inout_noarg_interface), deferred :: init                    !< Grid and stencil initialization - done once for the grid and stencil
      procedure(inout_noarg_interface), deferred :: setup                   !< Solver setup (every time the operator changes)
      procedure(inout_noarg_interface), deferred :: solve                   !< Execute solver (assumes new RHS and initial guess at every call)
      procedure(inout_noarg_interface), deferred :: destroy                 !< Solver destruction (every time the operator changes)
      procedure                                  :: apply_bcond             !< Modify the RHS of the system for boundary conditions
   end type linsol
   
   
   !> Interface
   abstract interface
      subroutine in_noarg_interface(this)
         import linsol
         class(linsol), intent(in) :: this
      end subroutine in_noarg_interface
      subroutine inout_noarg_interface(this)
         import linsol
         class(linsol), intent(inout) :: this
      end subroutine inout_noarg_interface
   end interface


contains


   !> Impose boundary conditions to the RHS
   subroutine apply_bcond(this)
      implicit none
      class(linsol), intent(inout) :: this
      integer :: i,j,k,st,stm,ist
      ! X-direction
      if (this%cfg%nx.gt.1.and.this%cfg%xper.eqv..false.) then
         ! xm boundary
         if (this%cfg%imin_.eq.this%cfg%imin) then
            ! Get the minimum stencil shift
            stm=minval(this%stc(:,1))
            ! Loop over the cells
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do ist=0,-stm-1,+1
                     i=this%cfg%imin_+ist
                     do st=-1-ist,stm,-1
                        this%rhs(i,j,k)=this%rhs(i,j,k)-this%opr(this%stmap(st,0,0),i,j,k)*this%sol(i+st,j,k)
                     end do
                  end do
               end do
            end do
         end if
         ! xp boundary
         if (this%cfg%imax_.eq.this%cfg%imax) then
            ! Get the maximum stencil shift
            stm=maxval(this%stc(:,1))
            ! Loop over the cells
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do ist=0,-stm+1,-1
                     i=this%cfg%imax_+ist
                     do st=+1-ist,stm,+1
                        this%rhs(i,j,k)=this%rhs(i,j,k)-this%opr(this%stmap(st,0,0),i,j,k)*this%sol(i+st,j,k)
                     end do
                  end do
               end do
            end do
         end if
      end if
      ! Y-direction
      if (this%cfg%ny.gt.1.and.this%cfg%yper.eqv..false.) then
         ! ym boundary
         if (this%cfg%jmin_.eq.this%cfg%jmin) then
            ! Get the minimum stencil shift
            stm=minval(this%stc(:,2))
            ! Loop over the cells
            do k=this%cfg%kmin_,this%cfg%kmax_
               do ist=0,-stm-1,+1
                  j=this%cfg%jmin_+ist
                  do st=-1-ist,stm,-1
                     do i=this%cfg%imin_,this%cfg%imax_
                        this%rhs(i,j,k)=this%rhs(i,j,k)-this%opr(this%stmap(0,st,0),i,j,k)*this%sol(i,j+st,k)
                     end do
                  end do
               end do
            end do
         end if
         ! yp boundary
         if (this%cfg%jmax_.eq.this%cfg%jmax) then
            ! Get the maximum stencil shift
            stm=maxval(this%stc(:,2))
            ! Loop over the cells
            do k=this%cfg%kmin_,this%cfg%kmax_
               do ist=0,-stm+1,-1
                  j=this%cfg%jmax_+ist
                  do st=+1-ist,stm,+1
                     do i=this%cfg%imin_,this%cfg%imax_
                        this%rhs(i,j,k)=this%rhs(i,j,k)-this%opr(this%stmap(0,st,0),i,j,k)*this%sol(i,j+st,k)
                     end do
                  end do
               end do
            end do
         end if
      end if
      ! Z-direction
      if (this%cfg%nz.gt.1.and.this%cfg%zper.eqv..false.) then
         ! zm boundary
         if (this%cfg%kmin_.eq.this%cfg%kmin) then
            ! Get the minimum stencil shift
            stm=minval(this%stc(:,3))
            ! Loop over the cells
            do ist=0,-stm-1,+1
               k=this%cfg%kmin_+ist
               do st=-1-ist,stm,-1
                  do j=this%cfg%jmin_,this%cfg%jmax_
                     do i=this%cfg%imin_,this%cfg%imax_
                        this%rhs(i,j,k)=this%rhs(i,j,k)-this%opr(this%stmap(0,0,st),i,j,k)*this%sol(i,j,k+st)
                     end do
                  end do
               end do
            end do
         end if
         ! zp boundary
         if (this%cfg%kmax_.eq.this%cfg%kmax) then
            ! Get the maximum stencil shift
            stm=maxval(this%stc(:,3))
            ! Loop over the cells
            do ist=0,-stm+1,-1
               k=this%cfg%kmax_+ist
               do st=+1-ist,stm,+1
                  do j=this%cfg%jmin_,this%cfg%jmax_
                     do i=this%cfg%imin_,this%cfg%imax_
                        this%rhs(i,j,k)=this%rhs(i,j,k)-this%opr(this%stmap(0,0,st),i,j,k)*this%sol(i,j,k+st)
                     end do
                  end do
               end do
            end do
         end if
      end if
   end subroutine apply_bcond
   
   
end module linsol_class
