module rotorDisk_class
   use blade_class, only: blade
   use precision, only: WP
   use config_class, only: config
   implicit none
   private

   public :: rotorDisk, cart2cyl_vec, cyl2cart_vec

   ! In ori cartesian coordinates, x-/y-/z- are the same as config mesh
   ! In rotor disk cartesian coordinates, the z- dir is the axis
   ! the y- dir is the reference direction, the x- dir is the cross product of the axis and reference direction
   ! arctan(y/x) = theta  (rad)
   ! sqrt(x^2 + y^2) = r
   real(WP), dimension(3,3) :: rotMat        ! Rotation matrix from ori cartesian to rotor disk cartesian coordinates
   real(WP), dimension(3,3) :: invRotMat     ! Rotation matrix from rotor disk cartesian to ori cartesian coordinates

   type :: rotorDisk
      ! User-defined rotor disk properties

      real(WP) :: minR                            ! Minimum radius of the rotor disk
      real(WP) :: maxR                            ! Maximum radius of the rotor disk
      real(WP) :: omega                           ! Rotor angular velocity
      integer :: nblades                          ! Number of blades

      real(WP), dimension(3) :: center            ! Rotor center position
      real(WP), dimension(3) :: axis              ! Rotor axis vector, must be normalized and in x-/y-/z- direction
      real(WP), dimension(3) :: ref_dir           ! Reference direction for theta calculation

      ! Internal rotor disk properties

      type(blade) :: bl
      type(config) :: cfg

   contains

      procedure :: prepareRotorDisk         ! Prepare the rotor disk properties
      procedure :: setFaceArea              ! Set the blade face area in every cell
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
   end function rotorDisk_constructor

   subroutine prepareRotorDisk(self)
      use mathtools, only: cross_product
      use mathtools, only: normalize
      use mathtools, only: inverse_matrix
      implicit none
      class(rotorDisk) :: self

      real(WP), dimension(3,3) :: Ra, Rb, temp
      real(WP), dimension(3) :: x_

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


   end subroutine prepareRotorDisk

   ! Note: this function is to transform a vector from ori cartesian (Mesh coordinates) to rotor disk cylindric coordinates
   function cart2cyl_vec(carVec) result(CylVec)
      use mathtools, only: cross_product
      use mathtools, only: arctan
      implicit none
      real(WP), dimension(3), intent(in) :: carVec
      real(WP), dimension(3) :: CylVec     ! r, theta, z

      real(WP), dimension(3) :: temp

      ! Rotate the vector from ori cartesian to rotor disk cartesian coordinates
      temp = matmul(rotMat, carVec)

      CylVec(1) = sqrt(temp(1)**2 + temp(2)**2)
      CylVec(2) = arctan(temp(1), temp(2))
      CylVec(3) = temp(3)
   end function cart2cyl_vec

   ! Note: this function is to transform a vector from rotor disk cylindric to ori cartesian coordinates (Mesh coordinates)
   function cyl2cart_vec(CylVec) result(carVec)
      use mathtools, only: cross_product
      use mathtools, only: arctan
      implicit none
      real(WP), dimension(3), intent(in) :: CylVec    ! r, theta, z
      real(WP), dimension(3) :: carVec                ! x, y, z

      real(WP), dimension(3) :: temp                  ! x, y, z

      ! transform the vector from cylindrical to cartesian coordinates
      temp(1) = CylVec(1) * cos(CylVec(2))
      temp(2) = CylVec(2) * sin(CylVec(2))
      temp(3) = CylVec(3)

      ! Rotate the vector from rotor disk cartesian to ori cartesian coordinates
      carVec = matmul(invRotMat, temp)
   end function cyl2cart_vec

   subroutine setFaceArea(self, area_)
      implicit none
      class(rotorDisk) :: self
      real(WP), dimension(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_), intent(inout) :: area_


      integer :: i, j, k
      integer :: mask(1), origin_index
      logical :: is_inside
      real(WP) :: m1, m2, n1, n2
      integer, dimension(4) :: insideCircle

      ! calculate the BET area in every cell
      area_ = 0.0_WP

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
                        if (((m1**2 + n1**2) .lt. self%maxR**2) .and. ((m1**2 + n1**2) .gt. self%minR**2)) insideCircle(1) = 1
                        if (((m2**2 + n1**2) .lt. self%maxR**2) .and. ((m2**2 + n1**2) .gt. self%minR**2)) insideCircle(2) = 1
                        if (((m1**2 + n2**2) .lt. self%maxR**2) .and. ((m1**2 + n2**2) .gt. self%minR**2)) insideCircle(3) = 1
                        if (((m2**2 + n2**2) .lt. self%maxR**2) .and. ((m2**2 + n2**2) .gt. self%minR**2)) insideCircle(4) = 1

                        ! estimate the area of the circle
                        area_(i,j,k) = sum(insideCircle) * self%cfg%dy(j) * self%cfg%dz(k) / 4.0_WP
                     end if
                  else if (mask(1)==2) then
                     if (j == origin_index) then
                        ! if the disk is in x-z plane

                        ! calculate the relative position
                        m1 = self%cfg%x(i) - self%center(1) ; m2 = self%cfg%x(i+1) - self%center(1) ! realtive x position
                        n1 = self%cfg%z(k) - self%center(3) ; n2 = self%cfg%z(k+1) - self%center(3) ! realtive z position

                        ! determine number of points inside the circle
                        insideCircle = 0
                        if (((m1**2 + n1**2) .lt. self%maxR**2) .and. ((m1**2 + n1**2) .gt. self%minR**2)) insideCircle(1) = 1
                        if (((m2**2 + n1**2) .lt. self%maxR**2) .and. ((m2**2 + n1**2) .gt. self%minR**2)) insideCircle(2) = 1
                        if (((m1**2 + n2**2) .lt. self%maxR**2) .and. ((m1**2 + n2**2) .gt. self%minR**2)) insideCircle(3) = 1
                        if (((m2**2 + n2**2) .lt. self%maxR**2) .and. ((m2**2 + n2**2) .gt. self%minR**2)) insideCircle(4) = 1

                        ! estimate the area of the circle
                        area_(i,j,k) = sum(insideCircle) * self%cfg%dx(i) * self%cfg%dz(k) / 4.0_WP
                     end if
                  else
                     if (k == origin_index) then
                        ! if the disk is in x-y plane

                        ! calculate the relative position
                        m1 = self%cfg%x(i) - self%center(1) ; m2 = self%cfg%x(i+1) - self%center(1) ! realtive x position
                        n1 = self%cfg%y(j) - self%center(2) ; n2 = self%cfg%y(j+1) - self%center(2) ! realtive y position

                        ! determine number of points inside the circle
                        insideCircle = 0
                        if (((m1**2 + n1**2) .lt. self%maxR**2) .and. ((m1**2 + n1**2) .gt. self%minR**2)) insideCircle(1) = 1
                        if (((m2**2 + n1**2) .lt. self%maxR**2) .and. ((m2**2 + n1**2) .gt. self%minR**2)) insideCircle(2) = 1
                        if (((m1**2 + n2**2) .lt. self%maxR**2) .and. ((m1**2 + n2**2) .gt. self%minR**2)) insideCircle(3) = 1
                        if (((m2**2 + n2**2) .lt. self%maxR**2) .and. ((m2**2 + n2**2) .gt. self%minR**2)) insideCircle(4) = 1

                        ! estimate the area of the circle
                        area_(i,j,k) = sum(insideCircle) * self%cfg%dx(i) * self%cfg%dy(j) / 4.0_WP
                     end if
                  end if
               end do
            end do
         end do
      end if

   end subroutine setFaceArea

   subroutine calculateForce(self, area, rho, U, V, W, forceX, forceY, forceZ)
      use mathtools, only: arctan
      use mathtools, only: twoPi
      implicit none
      class(rotorDisk) :: self
      real(WP), dimension(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_), intent(in) :: area, rho, U, V, W
      real(WP), dimension(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_), intent(out) :: forceX, forceY, forceZ

      real(WP) :: relaPos(3), carVelocity(3), cylVelocity(3)
      integer :: i, j, k
      integer :: mask(1)
      real(WP) :: radius, chord, twist, cl, cd
      real(WP) :: alphaGeom, alphaEff              ! Geometric and effective angle of attack
      real(WP) :: pDyn, f
      real(WP), dimension(3) :: localforce

      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               if (area(i,j,k) .gt. 0.0_WP) then
                  ! Get the radius, from cell center to the rotor disk center
                  mask = maxloc(abs(self%axis))
                  relaPos = [self%cfg%xm(i) - self%center(1),self%cfg%ym(j) - self%center(2),self%cfg%zm(k) - self%center(3)]
                  relaPos(mask(1)) = 0.0_WP
                  radius = sqrt(sum(relaPos**2))

                  ! Transform the velocity from ori cartesian to rotor disk cylindrical coordinates
                  carVelocity = [U(i,j,k), V(i,j,k), W(i,j,k)]
               
                  cylVelocity = cart2cyl_vec(carVelocity)


                  ! Set Radical component to be zero
                  cylVelocity(1) = 0.0_WP

                  ! Set blade normal component (theta) of velocity
                  cylVelocity(2) = self%omega * radius - cylVelocity(2)

                  ! Determine blade data for this radius
                  call self%bl%interpolate_r(radius, chord, twist)

                  ! Effective angle of attack
                  alphaGeom = twist                               ! Twist angle in degrees
                  alphaEff = alphaGeom - arctan(sign(1.0_WP, self%omega)*cylVelocity(2), cylVelocity(3))*360.0_WP/twoPi    ! Effective angle of attack in degrees

                  ! Calculate the lift and drag coefficients
                  call self%bl%interpolate_a(alphaEff, cl, cd)

                  ! Tip effect correction, ignore for now

                  ! Calculate forces perpendicular to blade
                  pDyn = 0.5_WP * rho(i,j,k) * sum(cylVelocity**2)
                  f = pDyn * chord * self%nblades * area(i,j,k) / radius / twoPi

                  localforce = [0.0_WP, sign(1.0_WP, self%omega)*-1.0_WP*f*cd, f*cl]      ! r, theta, z

                  ! Transform the force from rotor disk cylindrical to ori cartesian coordinates
                  ! Note : the force transformation does not consider the trimming rotation of the rotor disk (LRF - RSP, common in helicopter and so on)
                  localforce = cyl2cart_vec(localforce)
                  
                  ! Transform the force into volumetric force
                  localforce = localforce / self%cfg%vol(i,j,k)

                  forceX(i,j,k) = localforce(1)
                  forceY(i,j,k) = localforce(2)
                  forceZ(i,j,k) = localforce(3)
                  
               end if
            end do
         end do
      end do

   end subroutine calculateForce


end module rotorDisk_class
