submodule(hermite) hermite5_implementation

contains

  module procedure H5_interpolate

    implicit none

    real(wp) :: fp0, fp1, fpp0, fpp1, xval, dx, dx2, dxinv
    real(wp), dimension(2) :: derivs

    if ( present(dydx) .and. present(d2ydx2) ) then
      fp0 = dydx(0)
      fp1 = dydx(1)
      fpp0 = d2ydx2(0)
      fpp1 = d2ydx2(1)
    else
      derivs = calc_deriv12 ( x(-5:0), y(-5:0) )
      fp0 = derivs(1)
      fpp0 = derivs(2)
      derivs = calc_deriv12 ( x(-4:1), y(-4:1) )
      fp1 = derivs(1)
      fpp1 = derivs(2)
    end if
!    derivs = calc_deriv12 ( x(-5:0), y(-5:0) )
!    print*, 'x = ', x
!    print*, 'y = ', y
!    print*,'der1 = ', fp0, fpp0, derivs(1), derivs(2)
!    derivs = calc_deriv12 ( x(-4:1), y(-4:1) )
!    print*,'der2 = ', fp1, fpp1, derivs(1), derivs(2)

!    stop
    dx = x(1) - x(0)
    dx2 = dx*dx
    dxinv = 1.0_wp / dx
    xval = ( xv - x(0) ) * dxinv

    res(1) = H5_P0 ( xval ) * y(0) &
           + H5_P1 ( xval ) * dx * fp0 &
           + H5_P2 ( xval ) * dx2 * fpp0 &
           + H5_P0 ( 1.0_wp - xval ) * y(1) &
           - H5_P1 ( 1.0_wp - xval ) * dx * fp1 &
           + H5_P2 ( 1.0_wp - xval ) * dx2 * fpp1
    res(2) = ( H5_P01 ( xval ) * y(0) &
             + H5_P11 ( xval ) * dx * fp0 &
             + H5_P21 ( xval ) * dx2 * fpp0 &
             - H5_P01 ( 1.0_wp - xval ) * y(1) &
             + H5_P11 ( 1.0_wp - xval ) * dx * fp1 &
             - H5_P21 ( 1.0_wp - xval ) * dx2 * fpp1 ) * dxinv
    res(3) = ( H5_P02 ( xval ) * y(0) &
             + H5_P12 ( xval ) * dx * fp0 &
             + H5_P22 ( xval ) * dx2 * fpp0 &
             + H5_P02 ( 1.0_wp - xval ) * y(1) &
             - H5_P12 ( 1.0_wp - xval ) * dx * fp1 &
             + H5_P22 ( 1.0_wp - xval ) * dx2 * fpp1 ) * dxinv**2

  end procedure H5_interpolate


  function calc_deriv12 ( x, y )

    implicit none

    real(wp), dimension(6) :: x
    real(wp), dimension(6) :: y
    real(wp), dimension(2) :: calc_deriv12
    real(wp) :: dx12, dx13, dx14, dx15, dx16, dx23, dx24, dx25, dx26, &
                dx34, dx35, dx36, dx45, dx46, dx56

    dx12 = x(1) - x(2)
    dx13 = x(1) - x(3)
    dx14 = x(1) - x(4)
    dx15 = x(1) - x(5)
    dx16 = x(1) - x(6)
    dx23 = x(2) - x(3)
    dx24 = x(2) - x(4)
    dx25 = x(2) - x(5)
    dx26 = x(2) - x(6)
    dx34 = x(3) - x(4)
    dx35 = x(3) - x(5)
    dx36 = x(3) - x(6)
    dx45 = x(4) - x(5)
    dx46 = x(4) - x(6)
    dx56 = x(5) - x(6)

    calc_deriv12(1) = dx26*dx36*dx46*dx56/(dx12*dx13*dx14*dx15*dx16)*y(1) &
                     -dx16*dx36*dx46*dx56/(dx12*dx23*dx24*dx25*dx26)*y(2) &
                     +dx16*dx26*dx46*dx56/(dx13*dx23*dx34*dx35*dx36)*y(3) &
                     -dx16*dx26*dx36*dx56/(dx14*dx24*dx34*dx45*dx46)*y(4) &
                     +dx16*dx26*dx36*dx46/(dx15*dx25*dx35*dx45*dx56)*y(5) &
                     -(1.0_wp/dx16+1.0_wp/dx26+1.0_wp/dx36+1.0_wp/dx46 &
                      +1.0_wp/dx56)*y(6)

    calc_deriv12(2) = ( - 2.0_wp*x(3)*x(4)*x(5) - 2.0_wp*x(2)*( x(4)*x(5) &
                           + x(3)*( x(4) + x(5) ) ) &
                        + 4.0_wp*( x(4)*x(5) + x(3)*( x(4) + x(5) ) &
                           + x(2)*( x(3) + x(4) + x(5) ) )*x(6) &
                        - 6.0_wp*( x(2) + x(3) + x(4) + x(5) )*x(6)**2 &
                        + 8*x(6)**3 ) &
                    / ( dx12*dx13*dx14*dx15*dx16 ) * y(1)
    calc_deriv12(2) = calc_deriv12(2) &
                    + ( 2.0_wp*( x(1)*x(3)*x(4) + x(1)*x(3)*x(5) &
                               + x(1)*x(4)*x(5) + x(3)*x(4)*x(5) &
                               - 2.0_wp*( x(4)*x(5) + x(3)*( x(4) + x(5) ) &
                                        + x(1)*( x(3) + x(4) + x(5) ) )*x(6) &
                               + 3.0_wp*( x(1) + x(3) + x(4) + x(5) )*x(6)**2 &
                               - 4.0_wp*x(6)**3 ) ) &
                    / ( dx12*dx23*dx24*dx25*dx26 ) * y(2)
    calc_deriv12(2) = calc_deriv12(2) &
                    - ( 2.0_wp*( x(1)*x(2)*x(4) + x(1)*x(2)*x(5) &
                               + x(1)*x(4)*x(5) + x(2)*x(4)*x(5) &
                               - 2.0_wp*( x(4)*x(5) + x(2)*( x(4) + x(5) ) &
                                        + x(1)*( x(2) + x(4) + x(5) ) )*x(6) &
                               + 3.0_wp*( x(1) + x(2) + x(4) + x(5) )*x(6)**2 &
                               - 4.0_wp*x(6)**3 ) ) &
                    / ( dx13*dx23*dx34*dx35*dx36 ) * y(3)
    calc_deriv12(2) = calc_deriv12(2) &
                    + ( 2.0_wp*( x(1)*x(2)*x(3) + x(1)*x(2)*x(5) &
                               + x(1)*x(3)*x(5) + x(2)*x(3)*x(5) &
                               - 2.0_wp*( x(3)*x(5) + x(2)*( x(3) + x(5) ) &
                                        + x(1)*( x(2) + x(3) + x(5) ) )*x(6) &
                               + 3.0_wp*( x(1) + x(2) + x(3) + x(5) )*x(6)**2 &
                               - 4.0_wp*x(6)**3 ) ) &
                    / ( dx14*dx24*dx34*dx45*dx46 ) * y(4)
    calc_deriv12(2) = calc_deriv12(2) &
                    - (2.0_wp*( x(1)*x(2)*x(3) + x(1)*x(2)*x(4) &
                              + x(1)*x(3)*x(4) + x(2)*x(3)*x(4) &
                              - 2.0_wp*( x(3)*x(4) + x(2)* (x(3) + x(4) ) &
                                       + x(1)*( x(2) + x(3) + x(4) ) )*x(6) &
                              + 3.0_wp*( x(1) + x(2) + x(3) + x(4) )*x(6)**2 &
                              - 4.0_wp*x(6)**3 ) ) &
                    / ( dx15*dx25*dx35*dx45*dx56 ) * y(5)
    calc_deriv12(2) = calc_deriv12(2) &
                    + ( 2.0_wp*( x(2)*x(3)*x(4) + x(2)*x(3)*x(5) &
                               + x(2)*x(4)*x(5) + x(3)*x(4)*x(5) &
                               - 3.0_wp*( x(3)*x(4) + ( x(3) + x(4) )*x(5) &
                                        + x(2)*( x(3) + x(4) + x(5) ) )*x(6) &
                               + 6.0_wp*( x(2) + x(3) + x(4) + x(5) )*x(6)**2 &
                               - 10.0_wp*x(6)**3 + x(1)* ( x(3)*x(4) &
                                                         + x(3)*x(5) &
                                                         + x(4)*x(5) &
                                                         + x(2)*( x(3) + x(4) &
                                                                + x(5) &
                                                                - 3.0_wp*x(6) &
                                                                ) &
                                                    - 3.0_wp*( x(3) + x(4) &
                                                             + x(5) )*x(6) &
                                                    + 6.0_wp*x(6)**2 ) ) ) &
                    / ( dx16*dx26*dx36*dx46*dx56 ) * y(6)

  end function calc_deriv12


  function H5_P0 ( x )

    implicit none

    real(wp) :: H5_P0
    real(wp), intent(in) :: x

    ! P0 = -(-1+x)^3*(1+3*x+6*x^2)
    H5_P0 = -(-1.0_wp + x)**3*(1.0_wp + x*(3.0_wp + 6.0_wp*x))

  end function H5_P0


  function H5_P01 ( x )

    implicit none

    real(wp) :: H5_P01
    real(wp), intent(in) :: x

    ! P01 = -30*(-1+x^2)*x^2
    H5_P01 = -30.0_wp*(-1.0_wp + x)**2*x**2

  end function H5_P01


  function H5_P02 ( x )

    implicit none

    real(wp) :: H5_P02
    real(wp), intent(in) :: x

    ! P02 = -60*x*(1-3*x+2*x^2)
    H5_P02 = -60.0_wp*x*(1.0_wp + x*(-3.0_wp + 2.0_wp*x))

  end function H5_P02


  function H5_P1 ( x )

    implicit none

    real(wp) :: H5_P1
    real(wp), intent(in) :: x

    ! P1 = -(-1+x)^3*x*(1+3*x)
    H5_P1 = -(-1.0_wp + x)**3*x*(1.0_wp + 3.0_wp*x)

  end function H5_P1


  function H5_P11 ( x )

    implicit none

    real(wp) :: H5_P11
    real(wp), intent(in) :: x

    ! P11 = -(-1+x)^2*(-1-2*x+15*x^2)
    H5_P11 = -(-1.0_wp + x)**2*(-1.0_wp + x*(-2.0_wp + 15.0_wp*x))

  end function H5_P11


  function H5_P12 ( x )

    implicit none

    real(wp) :: H5_P12
    real(wp), intent(in) :: x

    ! P12 = -12*x*(3-8*x+5*x^2)
    H5_P12 = -12.0_wp*x*(3.0_wp + x*(-8.0_wp + 5.0_wp*x))

  end function H5_P12

  function H5_P2 ( x )

    implicit none

    real(wp) :: H5_P2
    real(wp), intent(in) :: x

    ! P2 = -1/2*(-1+x)^3*x^2
    H5_P2 = -0.5_wp*(-1.0_wp + x)**3*x**2

  end function H5_P2


  function H5_P21 ( x )

    implicit none

    real(wp) :: H5_P21
    real(wp), intent(in) :: x

    ! P21 = -1/2*(-1+x)**2*x*(-2+5*x)
    H5_P21 = -0.5_wp*(-1.0_wp + x)**2*x*(-2.0_wp+5.0_wp*x)

  end function H5_P21


  function H5_P22 ( x )

    implicit none

    real(wp) :: H5_P22
    real(wp), intent(in) :: x

    ! P22 = 1-9*x+18*x^2-10*x^3
    H5_P22 = 1.0_wp + x*(-9.0_wp + x*(18.0_wp - 10.0_wp*x))

  end function H5_P22

end submodule hermite5_implementation
