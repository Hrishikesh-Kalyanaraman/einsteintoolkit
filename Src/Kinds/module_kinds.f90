module kinds
!! Definition of the basic kind values and some basic constants.

  implicit none

  integer, parameter :: sp = selected_real_kind(5,30)
  !! Single precision floating point.
  integer, parameter :: dp = selected_real_kind(9,99)
  !! Double precision floating point.
  integer, parameter :: qp = selected_real_kind(20,199)
  !! Quad precision floating point.

  integer, parameter :: wp = dp
  !! The working precsion.
  integer, parameter :: ip = selected_int_kind(8)
  !! 32 bit integers.

! These empty arrays are used to initialize variables to either the min or
! max possible number of kind wp or integer.
  real(wp), dimension(2:1) :: empty
  !! Empty floating point array used to initialize variables to either the min
  !! or max value of kind wp.
  integer(ip), dimension(2:1) :: iempty
  !! Empty integer array used to initialize variables to either the min
  !! or max value of kind ip.

  real(wp), parameter :: rzero = 0.0_wp
  !! Zero real type constant.
  complex(wp), parameter :: czero = cmplx(0.0_wp,0.0_wp,wp)
  !! Zero complex type constant.
  integer(ip), parameter :: izero = 0_ip
  !! Zero integer type constant.
  complex(wp), parameter :: zi = cmplx(0.0_wp,1.0_wp,wp)
  !! The imaginary unit.
end module kinds
