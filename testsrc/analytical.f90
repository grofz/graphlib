  module analytical_mod
    use iso_fortran_env, only : dp => real64
    implicit none (type, external)
    private
    public analytical_1d_temperature

    real(dp), parameter :: PI = 3.141592653589793_dp
    integer, parameter  :: MAX_N = 200

  contains

    pure function analytical_1d_temperature(y, t, alfa, H, temp_bc) result(temp)
      real(dp), intent(in) :: y, t, alfa, H, temp_bc(2)
      real(dp) :: temp
!
! One dimensional heat conduction with Dirichlet boundary conditions on both
! sides, return temperature T(y, t)
!
! IN
!   y          - position,    0 <= y <= H
!   t          - time,        t > 0
!   alfa       - thermal diffusivity (i.e. conductivity / volumetric capacity)
!   H          - domain size
!   temp_bc    - the initial and boundary conditions
!                  T(y, t=0) = temp_bc(1)
!                  T(y=0, t) = temp_bc(1)
!                  T(y=H, t) = temp_bc(2)
!
      real(dp) :: delta_temp
      integer :: i

      if (t < 0.0_dp) error stop &
        'analytical_1d_temperature - time must have a positive value'
      if (y < 0.0_dp .or. y > H) error stop &
        'analytical_1d_temperature - y must be between 0 and H'

      delta_temp = temp_bc(2) - temp_bc(1)
      temp = temp_bc(1) + delta_temp * (1.0_dp - y/H)
      do i=1, MAX_N
        temp = temp - 2.0_dp/PI*delta_temp * sin(real(i,dp)*PI*y/H)/real(i,dp) &
            * exp(-real(i,dp)**2 * PI**2 * alfa * t / H**2)
      end do
    end function analytical_1d_temperature

  end module analytical_mod
