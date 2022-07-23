submodule(hermite) hermite3_implementation

contains

  module procedure H3_interpolate

    implicit none

    real(wp) :: fp0, fp1, xval, dx, dxinv

    if ( present(dydx) ) then
      fp0 = dydx(0)
      fp1 = dydx(1)
    else
      fp0 = calc_deriv1 ( x(-3:0), y(-3:0) )
      fp1 = calc_deriv1 ( x(-2:1), y(-2:1) )
    end if
    
    dx = x(1) - x(0)
    dxinv = 1.0_wp / dx
    xval = ( xv - x(0) ) * dxinv

    res(1)  = H3_P0 ( xval ) * y(0) &
             +H3_P1 ( xval ) * dx * fp0 &
             +H3_P0 ( 1.0_wp - xval ) * y(1) &
             -H3_P1 ( 1.0_wp - xval ) * dx * fp1
    res(2)  = ( H3_P01 ( xval ) * y(0) &
               +H3_P11 ( xval ) * dx * fp0 &
               -H3_P01 ( 1.0_wp - xval ) * y(1) &
               +H3_P11 ( 1.0_wp - xval ) * dx * fp1 ) * dxinv
    res(3)  = ( H3_P02 ( xval ) * y(0) &
               +H3_P12 ( xval ) * dx * fp0 &
               +H3_P02 ( 1.0_wp - xval ) * y(1) &
               -H3_P12 ( 1.0_wp - xval ) * dx * fp1 ) * dxinv**2

  end procedure H3_interpolate
 

  function calc_deriv1 ( x, y )

    implicit none

    real(wp), dimension(4) :: x
    real(wp), dimension(4) :: y
    real(wp) :: calc_deriv1
    real(wp) :: dx12, dx13, dx14, dx23, dx24, dx34

    dx12 = x(1) - x(2)
    dx13 = x(1) - x(3)
    dx14 = x(1) - x(4)
    dx23 = x(2) - x(3)
    dx24 = x(2) - x(4)
    dx34 = x(3) - x(4)

    calc_deriv1 = dx24*dx34/(dx12*dx13*dx14)*y(1) &
                 -dx14*dx34/(dx12*dx23*dx24)*y(2) &
                 +dx14*dx24/(dx13*dx23*dx34)*y(3) &
                 -(1.0_wp/dx14+1.0_wp/dx24+1.0_wp/dx34)*y(4)
    
  end function calc_deriv1


  function H3_P0 ( x )

    implicit none

    real(wp) :: H3_P0
    real(wp), intent(in) :: x

    ! P0 = 2*x^3-3*x^2+1
    H3_P0 = 1.0_wp + x**2*(-3.0_wp + 2.0_wp*x)

  end function H3_P0


  function H3_P01 ( x )

    implicit none

    real(wp) :: H3_P01
    real(wp), intent(in) :: x

    ! P01 = 6*x^2-6*x
    H3_P01 = 6.0_wp*x*(-1.0_wp+x)

  end function H3_P01


  function H3_P02 ( x )

    implicit none

    real(wp) :: H3_P02
    real(wp), intent(in) :: x

    ! P02 = 12*x-6
    H3_P02 = 6_wp*(-1.0_wp+2.0*x)

  end function H3_P02


  function H3_P1 ( x )

    implicit none

    real(wp) :: H3_P1
    real(wp), intent(in) :: x

    ! P1 = x^3-2*x^2+x
    H3_P1= x*(1.0_wp + x*(-2.0_wp + x))

  end function H3_P1


  function H3_P11 ( x )

    implicit none

    real(wp) :: H3_P11
    real(wp), intent(in) :: x

    ! P11 = 3*x^2-4*x+1
    H3_P11= 1.0_wp + x*(-4.0_wp + 3.0*x)

  end function H3_P11


  function H3_P12 ( x )

    implicit none

    real(wp) :: H3_P12
    real(wp), intent(in) :: x

    ! P12 = 6*x-4
    H3_P12= -4.0_wp + 6.0_wp*x

  end function H3_P12

end submodule
