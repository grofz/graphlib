  module analytical_mod
    use iso_fortran_env, only : dp => real64
    implicit none (type, external)
    private
    public conduction_temp

    real(dp), parameter :: PI = 3.141592653589793_dp

  contains

    pure function conduction_temp(y, t, alfa, H, temp_bc) result(temp)
      real(dp), intent(in) :: y, t, alfa, H, temp_bc(2)
      real(dp) :: temp
!
! One dimensional heat conduction, return temperature T(y, t)
! IN
!   y - position from 0 to H
!   t - time
!   alfa - heat diffusivity (conductivity / volumetric capacity)
!   H - domain size
!   temp_bc(1) - Dirichlet B.C for y=0, and the initial temperature
!   temp_bc(2) - Dirichlet B.C for y=H
!
      real(dp) :: temp_dif
      integer :: i
      integer, parameter :: MAX_N = 200

      temp_dif = temp_bc(2) - temp_bc(1)
      temp = temp_bc(1) + temp_dif * (1.0_dp - y/H)
      do i=1, MAX_N
        temp = temp - 2.0_dp/pi*temp_dif * sin(real(i,dp)*PI*y/H)/real(i,dp) &
            * exp(-real(i,dp)**2 * PI**2 * alfa * t / H**2)
      end do
    end function conduction_temp
  end module analytical_mod
