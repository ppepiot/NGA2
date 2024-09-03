module rotorDisk_class
   use blade_class, only: blade
   use precision, only: WP
   use config_class, only: config
   implicit none
   private

   public :: rotorDisk

   ! In global cartesian coordinates, x-/y-/z- are the same as config mesh
   ! In rotor disk cartesian coordinates, the z- dir is the axis
   ! the y- dir is the reference direction, the x- dir is the cross product of the axis and reference direction
   ! arctan(y/x) = theta  (rad)
   ! sqrt(x^2 + y^2) = r
   real(WP), dimension(3,3) :: rotMat        ! Rotation matrix from global cartesian to rotor disk cartesian coordinates
   real(WP), dimension(3,3) :: invRotMat     ! Rotation matrix from rotor disk cartesian to global cartesian coordinates

   type :: rotorDisk
      ! User-defined rotor disk properties

      real(WP) :: minR                            ! Minimum radius of the rotor disk
      real(WP) :: maxR                            ! Maximum radius of the rotor disk
      real(WP) :: omega                           ! Rotor angular velocity
      integer :: nblades                          ! Number of blades

      real(WP) :: tipEffectParam                  ! Tip effect parameter, over it*radius, no lift force

      real(WP), dimension(3) :: center            ! Rotor center position
      real(WP), dimension(3) :: axis              ! Rotor axis vector, must be normalized and in x-/y-/z- direction
      real(WP), dimension(3) :: ref_dir           ! Reference direction for theta calculation

      real(WP) :: rho_ref                         ! Reference density, for output

      ! Internal rotor disk properties

      type(blade) :: bl
      type(config) :: cfg

      real(WP), dimension(:,:,:,:), allocatable :: cylPos  ! Cylindrical position in the rotor disk
      real(WP), dimension(:,:,:), allocatable :: area    ! Blade face area in every cell
      real(WP), dimension(:,:,:), allocatable :: forceX  ! Volumetric force in x-direction
      real(WP), dimension(:,:,:), allocatable :: forceY  ! Volumetric force in y-direction
      real(WP), dimension(:,:,:), allocatable :: forceZ  ! Volumetric force in z-direction

      ! output
      logical :: output
      real(WP) :: dragEff, liftEff, powerEff

   contains

      procedure :: prepareRotorDisk         ! Prepare the rotor disk properties
      procedure :: setTipEffect             ! Set the tip effect correction
      procedure :: calculateForce           ! Calculate the volumetric force in every cell from the rotor disk, the force can be added to the momentum equation

   end type rotorDisk

   interface rotorDisk
      procedure rotorDisk_constructor
   end interface rotorDisk

contains

   function rotorDisk_constructor(bl, cfg) result(self)
      implicit none
      type(rotorDisk) :: self
      type(blade), intent(in) :: bl
      type(config), intent(in) :: cfg

      self%cfg = cfg
      self%bl = bl

      ! initialize rotor disk properties
      self%minR = 0.0_WP
      self%maxR = 0.0_WP
      self%omega = 0.0_WP
      self%nblades = 0
      self%center = 0.0_WP
      self%axis = 0.0_WP

      self%tipEffectParam = 1.0_WP

      self%output = .true.
      self%rho_ref = 1000.0_WP

      allocate(self%cylPos(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,3))
      allocate(self%area(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

      allocate(self%forceX(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(self%forceY(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(self%forceZ(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

   end function rotorDisk_constructor

   subroutine prepareRotorDisk(self)
      use mathtools, only: cross_product
      use mathtools, only: normalize
      use mathtools, only: inverse_matrix
      use mathtools, only: arctan
      use mathtools, only: twoPi, pi
      implicit none
      class(rotorDisk) :: self

      integer :: i, j, k

      real(WP), dimension(3,3) :: Ra, Rb, temp
      real(WP), dimension(3) :: x_

      integer :: mask(1), origin_index
      logical :: is_inside
      real(WP) :: m1, m2, n1, n2
      integer, dimension(4) :: insideCircle

      real(WP), dimension(3) :: globalPos, rotatedPos, movedPos ! positions in global, rotated, moved, and cylindrical coordinates

      ! Prepare the rotor disk properties
      ! check center, axis, and ref_dir
      if (sum(self%axis) .lt. 1.0E-9_WP) then
         write(*,*) '[Rotor Disk] Error: axis vector is not defined'
         stop
      end if

      if (sum(self%ref_dir) .lt. 1.0E-9_WP) then
         write(*,*) '[Rotor Disk] Error: reference direction vector is not defined'
         stop
      end if

      ! Normalize
      self%axis = normalize(self%axis)
      self%ref_dir = normalize(self%ref_dir)


      ! Calculate the rotation matrix from cartesian to cylindrical coordinates
      ! A stands for cartesian coordinates, Columns in Ra are the unit vectors of the cartesian coordinates
      ! B stands for cylindrical coordinates, Columns in Rb are the unit vectors of the cylindrical coordinates

      Ra = reshape([1.0_WP, 0.0_WP, 0.0_WP, &
         0.0_WP, 1.0_WP, 0.0_WP, &
         0.0_WP, 0.0_WP, 1.0_WP], [3,3])

      ! In rotor disk cartesian coordinates, the z- dir is the axis
      ! the y- dir is the reference direction, the x- dir is the cross product of the axis and reference direction
      ! arctan(y/x) = theta
      x_ = cross_product(self%axis, self%ref_dir)

      Rb = reshape([x_(1), self%ref_dir(1), self%axis(1), &
         x_(2), self%ref_dir(2), self%axis(2), &
         x_(3), self%ref_dir(3), self%axis(3)], [3,3])

      ! The rotation matrix from cartesian to cylindrical coordinates (A => B) is Rb * Ra^(-1)
      ! The rotation matrix from cylindrical to cartesian coordinates (B => A) is Ra * Rb^(-1)
      temp = transpose(Ra)
      rotMat = matmul(Rb, temp)

      temp = transpose(Rb)
      invRotMat = matmul(Ra, temp)

      ! calculate the BET area in every cell
      self%area = 0.0_WP

      ! Get the axis direction. If mask(1) is 1, the axis is x-direction,
      ! if mask(1) is 2, the axis is y-direction, if mask(1) is 3, the axis is z-direction
      mask = maxloc(abs(self%axis))
      is_inside = .false.

      ! Based on mask(axis direction), find the index of the origin
      if (mask(1) == 1) then
         ! loop over the x-direction
         do i=self%cfg%imin_,self%cfg%imax_
            if ((self%cfg%x(i) .le. self%center(1)) .and. (self%cfg%x(i+1) .gt. self%center(1))) then
               origin_index = i
               is_inside = .true.
            end if
         end do
      else if (mask(1) == 2) then
         ! loop over the y-direction
         do j=self%cfg%jmin_,self%cfg%jmax_
            if ((self%cfg%y(j) .le. self%center(2)) .and. (self%cfg%y(j+1) .gt. self%center(2))) then
               origin_index = j
               is_inside = .true.
            end if
         end do
      else
         ! loop over the z-direction
         do k=self%cfg%kmin_,self%cfg%kmax_
            if ((self%cfg%z(k) .le. self%center(3)) .and. (self%cfg%z(k+1) .gt. self%center(3))) then
               origin_index = k
               is_inside = .true.
            end if
         end do
      end if

      ! If not is_inside, the rotor disk does not intersect with the domain in this processor
      if (is_inside) then
         ! loop over the domain
         do k=self%cfg%kmin_,self%cfg%kmax_
            do j=self%cfg%jmin_,self%cfg%jmax_
               do i=self%cfg%imin_,self%cfg%imax_
                  if (mask(1)==1) then
                     if (i == origin_index) then
                        ! if the disk is in y-z plane

                        ! calculate the relative position
                        m1 = self%cfg%y(j) - self%center(2) ; m2 = self%cfg%y(j+1) - self%center(2) ! realtive y position
                        n1 = self%cfg%z(k) - self%center(3) ; n2 = self%cfg%z(k+1) - self%center(3) ! realtive z position

                        ! determine number of points inside the circle
                        insideCircle = 0
                        if (((m1**2 + n1**2) .le. self%maxR**2) .and. ((m1**2 + n1**2) .ge. self%minR**2)) insideCircle(1) = 1
                        if (((m2**2 + n1**2) .le. self%maxR**2) .and. ((m2**2 + n1**2) .ge. self%minR**2)) insideCircle(2) = 1
                        if (((m1**2 + n2**2) .le. self%maxR**2) .and. ((m1**2 + n2**2) .ge. self%minR**2)) insideCircle(3) = 1
                        if (((m2**2 + n2**2) .le. self%maxR**2) .and. ((m2**2 + n2**2) .ge. self%minR**2)) insideCircle(4) = 1

                        ! estimate the area of the circle
                        self%area(i,j,k) = sum(insideCircle) * self%cfg%dy(j) * self%cfg%dz(k) / 4.0_WP
                     end if
                  else if (mask(1)==2) then
                     if (j == origin_index) then
                        ! if the disk is in x-z plane

                        ! calculate the relative position
                        m1 = self%cfg%x(i) - self%center(1) ; m2 = self%cfg%x(i+1) - self%center(1) ! realtive x position
                        n1 = self%cfg%z(k) - self%center(3) ; n2 = self%cfg%z(k+1) - self%center(3) ! realtive z position

                        ! determine number of points inside the circle
                        insideCircle = 0
                        if (((m1**2 + n1**2) .le. self%maxR**2) .and. ((m1**2 + n1**2) .ge. self%minR**2)) insideCircle(1) = 1
                        if (((m2**2 + n1**2) .le. self%maxR**2) .and. ((m2**2 + n1**2) .ge. self%minR**2)) insideCircle(2) = 1
                        if (((m1**2 + n2**2) .le. self%maxR**2) .and. ((m1**2 + n2**2) .ge. self%minR**2)) insideCircle(3) = 1
                        if (((m2**2 + n2**2) .le. self%maxR**2) .and. ((m2**2 + n2**2) .ge. self%minR**2)) insideCircle(4) = 1

                        ! estimate the area of the circle
                        self%area(i,j,k) = sum(insideCircle) * self%cfg%dx(i) * self%cfg%dz(k) / 4.0_WP
                     end if
                  else
                     if (k == origin_index) then
                        ! if the disk is in x-y plane

                        ! calculate the relative position
                        m1 = self%cfg%x(i) - self%center(1) ; m2 = self%cfg%x(i+1) - self%center(1) ! realtive x position
                        n1 = self%cfg%y(j) - self%center(2) ; n2 = self%cfg%y(j+1) - self%center(2) ! realtive y position

                        ! determine number of points inside the circle
                        insideCircle = 0
                        if (((m1**2 + n1**2) .le. self%maxR**2) .and. ((m1**2 + n1**2) .ge. self%minR**2)) insideCircle(1) = 1
                        if (((m2**2 + n1**2) .le. self%maxR**2) .and. ((m2**2 + n1**2) .ge. self%minR**2)) insideCircle(2) = 1
                        if (((m1**2 + n2**2) .le. self%maxR**2) .and. ((m1**2 + n2**2) .ge. self%minR**2)) insideCircle(3) = 1
                        if (((m2**2 + n2**2) .le. self%maxR**2) .and. ((m2**2 + n2**2) .ge. self%minR**2)) insideCircle(4) = 1

                        ! estimate the area of the circle
                        self%area(i,j,k) = sum(insideCircle) * self%cfg%dx(i) * self%cfg%dy(j) / 4.0_WP
                     end if
                  end if
               end do
            end do
         end do
      end if

      call self%cfg%sync(self%area)

      ! Prepare the relative position in the rotor disk

      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               globalPos = [self%cfg%xm(i), self%cfg%ym(j), self%cfg%zm(k)]
               rotatedPos = matmul(rotMat, globalPos)       ! rotate the position to disk cartesian direction
               movedPos = rotatedPos + matmul(rotMat, self%center)         ! move the origin position
               self%cylPos(i,j,k,1) = sqrt(movedPos(1)**2 + movedPos(2)**2)     ! r
               self%cylPos(i,j,k,2) = arctan(movedPos(1), movedPos(2))          ! theta
               self%cylPos(i,j,k,3) = movedPos(3)                               ! z
            end do
         end do
      end do
      
      ! Prepare the angular velocity in rad/s
      self%omega = self%omega * 2.0_WP * pi / 60.0_WP 

   end subroutine prepareRotorDisk

   subroutine setTipEffect(self, tf)
      implicit none
      class(rotorDisk) :: self
      real(WP), intent(in) :: tf

      self%tipEffectParam = tf

   end subroutine setTipEffect


   subroutine calculateForce(self, rho, U, V, W)
      use mathtools, only: arctan
      use mathtools, only: cross_product
      use mathtools, only: twoPi, pi
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(rotorDisk) :: self
      real(WP), dimension(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_), intent(in) :: rho, U, V, W
  

      real(WP), dimension(3) :: cylPos
      real(WP), dimension(3) :: globalVel       ! U, V, W
      real(WP), dimension(3) :: CylVel          ! radial, tangential, axial
      real(WP), dimension(3,3) :: rotMat_vec
      integer :: i, j, k
      real(WP) :: radius, chord, twist, cl, cd
      real(WP) :: alphaGeom, alphaEff, alphaFlow              ! Geometric and effective angle of attack
      real(WP) :: tipFactor                        ! Tip effect factor
      real(WP) :: fz, ftheta, f, flift, fdrag
      real(WP), dimension(3) :: localforce
      real(WP) :: signOmega

      real(WP) :: mydragEff, myliftEff, mypowerEff
      integer :: ierr

      self%forceX = 0.0_WP
      self%forceY = 0.0_WP
      self%forceZ = 0.0_WP

      mydragEff = 0.0_WP
      myliftEff = 0.0_WP
      mypowerEff = 0.0_WP

      signOmega = sign(1.0_WP, self%omega)

      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               if (self%area(i,j,k) .gt. 1.0E-32_WP) then
                  ! Get cylindrical position
                  cylPos = self%cylPos(i,j,k,:)

                  radius = cylPos(1)

                  ! prepare vector rotation matrix, disk-cartesian to radial-tangential-axial
                  rotMat_vec = reshape([cos(cylPos(2)), -sin(cylPos(2)), 0.0_WP, &
                     sin(cylPos(2)), cos(cylPos(2)), 0.0_WP, &
                     0.0_WP, 0.0_WP, 1.0_WP], [3,3])
                  ! Transform the velocity from ori cartesian to rotor disk cylindrical coordinates
                  globalVel = [U(i,j,k), V(i,j,k), W(i,j,k)]
                  globalVel = matmul(rotMat, globalVel)        ! rotate the velocity to disk cartesian direction
                  CylVel = matmul(rotMat_vec, globalVel)       ! rotate again to radial-tangential-axial


                  ! Set Radical component to be zero
                  CylVel(1) = 0.0_WP

                  ! Set blade normal component (theta) of velocity
                  CylVel(2) = self%omega * radius - CylVel(2)

                  ! Determine blade data for this radius
                  call self%bl%interpolate_r(radius, chord, twist)

                  ! Effective angle of attack
                  alphaGeom = twist                               ! Twist angle in degrees, no trim angle
                  alphaFlow = arctan(CylVel(2), CylVel(3))/pi*180.0_WP

                  if (alphaFlow .gt. 180.0_WP) alphaFlow = alphaFlow - 360.0_WP

                  alphaEff = alphaGeom - alphaFlow


                  ! Calculate the lift and drag coefficients
                  call self%bl%interpolate_a(alphaEff, cl, cd)

                  ! Tip effect correction
                  tipFactor = 1.0_WP
                  if (radius .gt. (self%tipEffectParam*self%maxR)) tipFactor = 0.0_WP

                  ! Calculate forces
                  alphaFlow = alphaFlow / 180.0_WP * pi

                  f = 0.5*rho(i,j,k)*sum(CylVel**2)*chord*self%area(i,j,k)*self%nblades/radius/twoPi
                  fdrag = -1.0_WP*f*signOmega*cd
                  flift = tipFactor*f*cl

                  ftheta = flift*sin(alphaFlow) + fdrag*cos(alphaFlow)
                  fz = flift*cos(alphaFlow) - fdrag*sin(alphaFlow)
                  

                  if (radius .lt. 0.018_WP) then
                     print*, 'alphaEff:', alphaEff, "Cl:", cl, "Cd:", cd   
                     print*, 'fz:', fz
                     print*, CylVel
                  end if


                  ! Calculate the force in the rotor disk cylindrical coordinates
                  ! radial, tangential, axial
                  localforce = [0.0_WP, ftheta , fz]

                  

                  mydragEff = mydragEff + localforce(2)*self%rho_ref 
                  myliftEff = myliftEff + localforce(3)*self%rho_ref
                  mypowerEff = mypowerEff + localforce(2)*radius*signOmega*self%omega*self%rho_ref

                  ! Transform the force from rotor disk cylindrical to global cartesian coordinates
                  ! Note : the force transformation does not consider the trimming rotation of the rotor disk (LRF - RSP, common in helicopter and so on)
                  rotMat_vec = transpose(rotMat_vec)
                  localforce = matmul(rotMat_vec, localforce)
                  localforce = matmul(invRotMat, localforce)

                  ! Transform the force into volumetric force
                  localforce = localforce / self%cfg%vol(i,j,k)

                  self%forceX(i,j,k) = localforce(1)
                  ! self%forceY(i,j,k) = localforce(2)
                  ! self%forceZ(i,j,k) = localforce(3)

               end if
            end do
         end do
      end do

      call self%cfg%sync(self%forceX)
      call self%cfg%sync(self%forceY)
      call self%cfg%sync(self%forceZ)

      ! Synchronize the output
      call MPI_ALLREDUCE(mydragEff, self%dragEff, 1, MPI_REAL_WP, MPI_SUM, self%cfg%comm, ierr)
      call MPI_ALLREDUCE(myliftEff, self%liftEff, 1, MPI_REAL_WP, MPI_SUM, self%cfg%comm, ierr)
      call MPI_ALLREDUCE(mypowerEff, self%powerEff, 1, MPI_REAL_WP, MPI_SUM, self%cfg%comm, ierr)

   end subroutine calculateForce


end module rotorDisk_class
