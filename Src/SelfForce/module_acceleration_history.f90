module acceleration_history
!! Module that keeps track of the acceleration history in order to try
!! to get smoother data at intermediate Runge-Kutta substeps.
!!
!! The implementation is found in
!! [[submodule_acceleration_history_implementation.f90]].

  use kinds

  type :: accel_history
  !! A type to keep track of the acceleration history.
    integer(ip) :: nlevels
    !! The number of timelevels of accelerations to keep.
    integer(ip) :: extrap_order
    !! The order to use for extrapolation.
    integer(ip) :: nfilled
    !! The number of timelevels that has been stored already.
    real(wp), dimension(:), allocatable :: at
    !! A 1d real array of size nlevels to keep track
    !! of the history of \(a^t\).
    real(wp), dimension(:), allocatable :: ar
    !! A 1d real array of size nlevels to keep track
    !! of the history of \(a^r\).
    real(wp), dimension(:), allocatable :: aphi
    !! A 1d real array of size nlevels to keep track
    !! of the history of \(a^{\phi}\).
    real(wp), dimension(:), allocatable :: dt
    !! A 1d real array of size nlevels to keep track
    !! of the time difference \(\Delta t\) between timelevels.
    real(wp), dimension(:), allocatable :: t_rel
    !! A 1d real array of sliding times with the current time always zero.
    real(wp), dimension(:,:), allocatable, private :: xm
    !! Local work array for the least square fitting.
    real(wp), dimension(:), allocatable, private :: yv
    !! Local work array for the least square fitting.
    real(wp), dimension(:), allocatable, private :: cv
    !! Local work array for the least square fitting.
    real(wp), dimension(:,:), allocatable, private :: covm
    !! Local work array for the least square fitting.
  contains
    procedure :: cycle_timelevels
    !! Routine to cycle previous timelevels and storing data in the last
    !! timelevel.
    procedure :: extrapolate
    !! Extrapolate from stored timelevels by a given \(\Delta t\).
    procedure :: calc_smooth_derivs
    !! Calculate smooth derivatives using least square fitting.
    final :: deallocate_timelevels
    !! The finalizer that will deallocate the timelevel arrays.
  end type accel_history

  interface accel_history
  !! The constructor for the acceleration history class.
    module procedure init_timelevels
  end interface accel_history

  interface
    module function init_timelevels ( nlevels, extrap_order )
    !! The interface for the constructor for the acceleration history class.
      type(accel_history) :: init_timelevels
      !! The return object has to be of [[accel_history]] class.
      integer(ip), intent(in) :: nlevels
      !! The number of timelevels to store.
      integer(ip), intent(in) :: extrap_order
      !! The order of extrapolation to use.
    end function

    module subroutine deallocate_timelevels ( this )
    !! The interface for the finalizer for the [[accel_history]] class.
      type(accel_history), intent(inout) :: this
      !! The argument has to be of the [[accel_history]] class.
    end subroutine deallocate_timelevels

    module subroutine cycle_timelevels ( this, accel, dtime )
    !! The interface for the [[cycle_timelevels]] routine.
      class(accel_history), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), dimension(4), intent(in) :: accel
      !! The 4-acceleration to store in the last timelevel.
      real(wp), intent(in) :: dtime
      !! The \(\Delta t\) for the last timestep.
    end subroutine cycle_timelevels

    module function extrapolate ( this, dtime )
    !! The interface for the [[extrapolate]] function.
      class(accel_history), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: dtime
      !! The \(\Delta t\) since the last stored timelevel.
      real(wp), dimension(4,3) :: extrapolate
      !! The return value is the 4-acceleration and it's first and second
      !! time derivatives extrapolated \(\Delta t\) into the future from the
      !! last stored timelevel.
    end function extrapolate

    module function calc_smooth_derivs ( this, x, y, dx, order )
    !! The interface for the [[calc_smooth_derivs]] function.
      class(accel_history), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), dimension(:), intent(in) :: x
      !! The x-values of tabulated data.
      real(wp), dimension(:), intent(in) :: y
      !! The y-values of tabulated data.
      real(wp), intent(in) :: dx
      !! The \(\Delta \) value (from the current zero) at which to evaluate
      !! The derivatives.
      integer(ip), intent(in) :: order
      real(wp), dimension(3) :: calc_smooth_derivs
      !! The return value is the function and it's first 2 derivatives
      !! evaluated at \(\Delta x\).
    end function calc_smooth_derivs

  end interface

  type(accel_history) :: ah

!  contains
!
!  function H3_interpolate ( x, y, xv )
!
!    implicit none
!
!    real(wp), dimension(3) :: H3_interpolate
!    real(wp), dimension(-2:1), intent(in) :: x, y
!    real(wp), intent(in) :: xv
!    real(wp) :: dxm10, dxm20, dxm2m1, dx01, dxm11
!    real(wp) :: fp0, fp1
!    real(wp) :: xval
!
!    dxm10 = x(-1) - x(0)
!    dxm20 = x(-2) - x(0)
!    dxm2m1 = x(-2) - x(-1)
!
!    dx01 = x(0) - x(1)
!    dxm11 = x(-1) - x(1)
!
!    fp0 = - dxm10 * y(-2) / ( dxm2m1 * dxm20 ) &
!          + dxm20 * y(-1) / ( dxm2m1 * dxm10 ) &
!          + ( 2.0_wp*x(0) - x(-2) - x(-1) ) * y(0) / ( dxm10 * dxm20 )
!
!    fp1 = - dx01 * y(-1) / ( dxm10 * dxm11 ) &
!          + dxm11 * y(0) / ( dxm10 * dx01 ) &
!          + ( 2.0_wp*x(1) - x(-1) - x(0) ) * y(1) / ( dx01 * dxm11 )
!
!    xval = ( x(0) - xv ) / dx01
!
!    H3_interpolate(1)  = H3_P0 ( xval ) * y(0) - H3_P1 ( xval ) * dx01 * fp0 &
!                       + H3_P0 ( 1.0_wp - xval ) * y(1) &
!                       + H3_P1 ( 1.0_wp - xval ) * dx01 * fp1
!    H3_interpolate(2)  = ( H3_P01 ( xval ) * y(0) &
!                         - H3_P11 ( xval ) * dx01 * fp0 &
!                         - H3_P01 ( 1.0_wp - xval ) * y(1) &
!                         - H3_P11 ( 1.0_wp - xval ) * dx01 * fp1 ) / ( -dx01 )
!    H3_interpolate(3)  = ( H3_P02 ( xval ) * y(0) &
!                         - H3_P12 ( xval ) * dx01 * fp0 &
!                         + H3_P02 ( 1.0_wp - xval ) * y(1) &
!                         + H3_P12 ( 1.0_wp - xval ) * dx01 * fp1 ) / dx01**2 
!
!  end function H3_interpolate
!
!
!  function H3_P0 ( x )
!
!    implicit none
!
!    real(wp) :: H3_P0
!    real(wp), intent(in) :: x
!
!    ! P0 = 2*x^3-3*x^2+1
!    H3_P0 = 1.0_wp + x**2*(-3.0_wp + 2.0_wp*x)
!
!  end function H3_P0
!
!
!  function H3_P01 ( x )
!
!    implicit none
!
!    real(wp) :: H3_P01
!    real(wp), intent(in) :: x
!
!    ! P01 = 6*x^2-6*x
!    H3_P01 = 6.0_wp*x*(-1.0_wp+x)
!
!  end function H3_P01
!
!
!  function H3_P02 ( x )
!
!    implicit none
!
!    real(wp) :: H3_P02
!    real(wp), intent(in) :: x
!
!    ! P02 = 12*x-6
!    H3_P02 = 6_wp*(-1.0_wp+2.0*x)
!
!  end function H3_P02
!  function H3_P1 ( x )
!
!    implicit none
!
!    real(wp) :: H3_P1
!    real(wp), intent(in) :: x
!
!    ! P1 = x^3-2*x^2+x
!    H3_P1= x*(1.0_wp + x*(-2.0_wp + x))
!
!  end function H3_P1
!
!
!  function H3_P11 ( x )
!
!    implicit none
!
!    real(wp) :: H3_P11
!    real(wp), intent(in) :: x
!
!    ! P11 = 3*x^2-4*x+1
!    H3_P11= 1.0_wp + x*(-4.0_wp + 3.0*x)
!
!  end function H3_P11
!
!
!  function H3_P12 ( x )
!
!    implicit none
!
!    real(wp) :: H3_P12
!    real(wp), intent(in) :: x
!
!    ! P12 = 6*x-4
!    H3_P12= -4.0_wp + 6.0_wp*x
!
!  end function H3_P12
!
!  function H5_interpolate ( x, y, xv )
!
!    implicit none
!
!    real(wp), dimension(3) :: H5_interpolate
!    real(wp), dimension(-4:1), intent(in) :: x, y
!    real(wp), intent(in) :: xv
!    real(wp) :: dx0m1, dx0m2, dx0m3, dx0m4, dxm1m4, dxm2m4, dxm3m4, &
!                dxm1m3, dxm2m3, dxm1m2, dx10, dx1m1, dx1m2, dx1m3
!    real(wp) :: fp0, fpp0, fp1, fpp1
!    real(wp) :: xval
!
!    dx0m1 = x(0) - x(-1)
!    dx0m2 = x(0) - x(-2)
!    dx0m3 = x(0) - x(-3)
!    dx0m4 = x(0) - x(-4)
!    dxm1m4 = x(-1) - x(-4)
!    dxm2m4 = x(-2) - x(-4)
!    dxm3m4 = x(-3) - x(-4)
!    dxm1m3 = x(-1) - x(-3)
!    dxm2m3 = x(-2) - x(-3)
!    dxm1m2 = x(-1) - x(-2)
!    dx10 = x(1) - x(0)
!    dx1m1 = x(1) - x(-1)
!    dx1m2 = x(1) - x(-2)
!    dx1m3 = x(1) - x(-3)
!
!
!    fp0 = dx0m1*dx0m2*dx0m3/(dx0m4*dxm1m4*dxm2m4*dxm3m4)*y(-4) &
!         -dx0m1*dx0m2*dx0m4/(dx0m3*dxm1m3*dxm2m3*dxm3m4)*y(-3) &
!         +dx0m1*dx0m3*dx0m4/(dx0m2*dxm1m2*dxm2m3*dxm2m4)*y(-2) &
!         -dx0m2*dx0m3*dx0m4/(dx0m1*dxm1m2*dxm1m3*dxm1m4)*y(-1) &
!         -( ( -4.0_wp*x(0)**3 + x(-2)*x(-3)*x(-4) &
!             + 3.0_wp*x(0)**2*( x(-1) + x(-2) + x(-3) + x(-4) ) &
!             + x(-1)*( x(-3)*x(-4) + x(-2)*( x(-3) + x(-4) ) ) &
!             - 2.0_wp*x(0)*( x(-3)*x(-4) + x(-2)*( x(-3) + x(-4) ) &
!                           + x(-1)*( x(-2) + x(-3) + x(-4) ) ) ) &
!            / ( dx0m1*dx0m2*dx0m3*dx0m4 ) )*y(0)
!
!    fpp0 = 2.0_wp*( 3.0_wp*x(0)**2 + x(-2)*x(-3) + x(-1)* ( x(-2) + x(-3) ) &
!                   -2.0_wp*x(0)*( x(-1) + x(-2) + x(-3) ) ) &
!          /( dx0m4*dxm1m4*dxm2m4*dxm3m4 ) * y(-4) &
!          -2.0_wp*( 3.0_wp*x(0)**2 + x(-2)*x(-4) + x(-1)* ( x(-2) + x(-4) ) &
!                   -2.0_wp*x(0)*( x(-1) + x(-2) + x(-4) ) ) &
!          /( dx0m3*dxm1m3*dxm2m3*dxm3m4 ) * y(-3) &
!          +2.0_wp*( 3.0_wp*x(0)**2 + x(-3)*x(-4) + x(-1)* ( x(-3) + x(-4) ) &
!                   -2.0_wp*x(0)*( x(-1) + x(-3) + x(-4) ) ) &
!          /( dx0m2*dxm1m2*dxm2m3*dxm2m4 ) * y(-2) &
!          -2.0_wp*( 3.0_wp*x(0)**2 + x(-3)*x(-4) + x(-2)* ( x(-3) + x(-4) ) &
!                   -2.0_wp*x(0)*( x(-2) + x(-3) + x(-4) ) ) &
!          /( dx0m1*dxm1m2*dxm1m3*dxm1m4 ) * y(-1) &
!          +2.0_wp*(6.0_wp*x(0)**2 + x(-2)*x(-3) + x(-2)*x(-4) + x(-3)*x(-4) &
!                  +x(-1)*( x(-2) + x(-3) + x(-4) ) &
!                  -3.0_wp*x(0)*( x(-1)+x(-2)+x(-3)+x(-4) ) ) &
!          /( dx0m1*dx0m2*dx0m3*dx0m4 ) * y(0)
!
!    fp1 = dx10*dx1m1*dx1m2/(dx1m3*dx0m3*dxm1m3*dxm2m3)*y(-3) &
!         -dx10*dx1m1*dx1m3/(dx1m2*dx0m2*dxm1m2*dxm2m3)*y(-2) &
!         +dx10*dx1m2*dx1m3/(dx1m1*dx0m1*dxm1m2*dxm1m3)*y(-1) &
!         -dx1m1*dx1m2*dx1m3/(dx10*dx0m1*dx0m2*dx0m3)*y(0) &
!         +( ( 4.0_wp*x(1)**3 - x(-1)*x(-2)*x(-3) &
!             - 3.0_wp*x(1)**2*( x(0) + x(-1) + x(-2) + x(-3) ) &
!             - x(0)*( x(-2)*x(-3) + x(-1)*( x(-2) + x(-3) ) ) &
!             + 2.0_wp*x(1)*( x(-2)*x(-3) + x(-1)*( x(-2) + x(-3) ) &
!                           + x(0)*( x(-1) + x(-2) + x(-3) ) ) ) &
!            / ( dx10*dx1m1*dx1m2*dx1m3 ) )*y(1)
!
!    fpp1 = 2.0_wp*( 3.0_wp*x(1)**2 + x(-1)*x(-2) + x(0)* ( x(-1) + x(-2) ) &
!                   -2.0_wp*x(1)*( x(0) + x(-1) + x(-2) ) ) &
!          /( dx1m3*dx0m3*dxm1m3*dxm2m3 ) * y(-3) &
!          -2.0_wp*( 3.0_wp*x(1)**2 + x(-1)*x(-3) + x(0)* ( x(-1) + x(-3) ) &
!                   -2.0_wp*x(1)*( x(0) + x(-1) + x(-3) ) ) &
!          /( dx1m2*dx0m2*dxm1m2*dxm2m3 ) * y(-2) &
!          +2.0_wp*( 3.0_wp*x(1)**2 + x(-2)*x(-3) + x(0)* ( x(-2) + x(-3) ) &
!                   -2.0_wp*x(1)*( x(0) + x(-2) + x(-3) ) ) &
!          /( dx1m1*dx0m1*dxm1m2*dxm1m3 ) * y(-1) &
!          -2.0_wp*( 3.0_wp*x(1)**2 + x(-2)*x(-3) + x(-1)* ( x(-2) + x(-3) ) &
!                   -2.0_wp*x(1)*( x(-1) + x(-2) + x(-3) ) ) &
!          /( dx10*dx0m1*dx0m2*dx0m3 ) * y(0) &
!          +2.0_wp*(6.0_wp*x(1)**2 + x(-1)*x(-2) + x(-1)*x(-3) + x(-2)*x(-3) &
!                  +x(0)*( x(-1) + x(-2) + x(-3) ) &
!                  -3.0_wp*x(1)*( x(0)+x(-1)+x(-2)+x(-3) ) ) &
!          /( dx10*dx1m1*dx1m2*dx1m3 ) * y(1)
!
!    xval = ( xv - x(0) ) / dx10
!
!!      print*,'coeff1 = ', 2.0_wp*( 3.0_wp*x(1)**2 + x(-1)*x(-2) + x(0)* ( x(-1) + x(-2) ) &
!!                   -2.0_wp*x(1)*( x(0) + x(-1) + x(-2) ) ) &
!!          /( dx1m3*dx0m3*dxm1m3*dxm2m3 )
!!      print*,'coeff2 = ', -2.0_wp*( 3.0_wp*x(1)**2 + x(-1)*x(-3) + x(0)* ( x(-1) + x(-3) ) &
!!                   -2.0_wp*x(1)*( x(0) + x(-1) + x(-3) ) ) &
!!          /( dx1m2*dx0m2*dxm1m2*dxm2m3 )
!!      print*,'coeff3 = ', +2.0_wp*( 3.0_wp*x(1)**2 + x(-2)*x(-3) + x(0)* ( x(-2) + x(-3) ) &
!!                   -2.0_wp*x(1)*( x(0) + x(-2) + x(-3) ) ) &
!!          /( dx1m1*dx0m1*dxm1m2*dxm1m3 )
!!      print*,'coeff4 = ', -2.0_wp*( 3.0_wp*x(1)**2 + x(-2)*x(-3) + x(-1)* ( x(-2) + x(-3) ) &
!!                   -2.0_wp*x(1)*( x(-1) + x(-2) + x(-3) ) ) &
!!          /( dx10*dx0m1*dx0m2*dx0m3 )
!!      print*,'coeff5 = ', +2.0_wp*(6.0_wp*x(1)**2 + x(-1)*x(-2) + x(-1)*x(-3) + x(-2)*x(-3) &
!!                  +x(0)*( x(-1) + x(-2) + x(-3) ) &
!!                  -3.0_wp*x(1)*( x(0)+x(-1)+x(-2)+x(-3) ) ) &
!!          /( dx10*dx1m1*dx1m2*dx1m3 )
!!
!    H5_interpolate(1)  = H5_P0 ( xval ) * y(0) &
!                       + H5_P1 ( xval ) * dx10 * fp0 &
!                       + H5_P2 ( xval ) * dx10**2 * fpp0 &
!                       + H5_P0 ( 1.0_wp - xval ) * y(1) &
!                       - H5_P1 ( 1.0_wp - xval ) * dx10 * fp1 &
!                       + H5_P2 ( 1.0_wp - xval ) * dx10**2 * fpp1
!    H5_interpolate(2)  = ( H5_P01 ( xval ) * y(0) &
!                         + H5_P11 ( xval ) * dx10 * fp0 &
!                         + H5_P21 ( xval ) * dx10**2 * fpp0 &
!                         - H5_P01 ( 1.0_wp - xval ) * y(1) &
!                         + H5_P11 ( 1.0_wp - xval ) * dx10 * fp1 &
!                         - H5_P21 ( 1.0_wp - xval ) * dx10**2 * fpp1 ) / dx10
!    H5_interpolate(3)  = ( H5_P02 ( xval ) * y(0) &
!                         + H5_P12 ( xval ) * dx10 * fp0 &
!                         + H5_P22 ( xval ) * dx10**2 * fpp0 &
!                         + H5_P02 ( 1.0_wp - xval ) * y(1) &
!                         - H5_P12 ( 1.0_wp - xval ) * dx10 * fp1 &
!                         + H5_P22 ( 1.0_wp - xval ) * dx10**2 * fpp1 ) / dx10**2
!!    if (xv == 0.5_wp) then
!!      print*,'x = ', x
!!      print*,'y = ', y
!!      print*,'yp = ', 5.0_wp*x**4-4.0_wp*x**3+3.0_wp*x**2-2.0_wp*x+1.0_wp
!!      print*,'ypp = ', 20.0_wp*x**3-12.0_wp*x**2+6.0_wp*x-2.0_wp
!!      print*,'fp0 = ', fp0
!!      print*,'fp1 = ', fp1
!!      print*,'fpp0 = ', fpp0
!!      print*,'fpp1 = ', fpp1
!!      print*,'P = ', H5_interpolate(1)
!!      print*,'Pp = ', H5_interpolate(2)
!!      print*,'Ppp = ', H5_interpolate(3)
!!      print*,'xval = ', xval
!!      print*,'H5_P0(xval) = ', H5_P0 ( xval )
!!      print*,'H5_P1(xval) = ', H5_P1 ( xval )
!!      print*,'H5_P2(xval) = ', H5_P2 ( xval )
!!      print*,'H5_P0(1-xval) = ', H5_P0 ( 1.0_wp - xval )
!!      print*,'H5_P1(1-xval) = ', H5_P1 ( 1.0_wp - xval )
!!      print*,'H5_P2(1-xval) = ', H5_P2 ( 1.0_wp - xval )
!!      print*,'H5_P01(xval) = ', H5_P01 ( xval )
!!      print*,'H5_P11(xval) = ', H5_P11 ( xval )
!!      print*,'H5_P21(xval) = ', H5_P21 ( xval )
!!      print*,'H5_P01(1-xval) = ', H5_P01 ( 1.0_wp - xval )
!!      print*,'H5_P11(1-xval) = ', H5_P11 ( 1.0_wp - xval )
!!      print*,'H5_P21(1-xval) = ', H5_P21 ( 1.0_wp - xval )
!!      print*,'H5_P02(xval) = ', H5_P02 ( xval )
!!      print*,'H5_P12(xval) = ', H5_P12 ( xval )
!!      print*,'H5_P22(xval) = ', H5_P22 ( xval )
!!      print*,'H5_P02(1-xval) = ', H5_P02 ( 1.0_wp - xval )
!!      print*,'H5_P12(1-xval) = ', H5_P12 ( 1.0_wp - xval )
!!      print*,'H5_P22(1-xval) = ', H5_P22 ( 1.0_wp - xval )
!!      print*
!!    endif
!
!  end function H5_interpolate
!
!
!  function H5_P0 ( x )
!
!    implicit none
!
!    real(wp) :: H5_P0
!    real(wp), intent(in) :: x
!
!    ! P0 = -(-1+x)^3*(1+3*x+6*x^2)
!    H5_P0 = -(-1.0_wp + x)**3*(1.0_wp + x*(3.0_wp + 6.0_wp*x))
!
!  end function H5_P0
!
!
!  function H5_P01 ( x )
!
!    implicit none
!
!    real(wp) :: H5_P01
!    real(wp), intent(in) :: x
!
!    ! P01 = -30*(-1+x^2)*x^2
!    H5_P01 = -30.0_wp*(-1.0_wp + x)**2*x**2
!
!  end function H5_P01
!
!
!  function H5_P02 ( x )
!
!    implicit none
!
!    real(wp) :: H5_P02
!    real(wp), intent(in) :: x
!
!    ! P02 = -60*x*(1-3*x+2*x^2)
!    H5_P02 = -60.0_wp*x*(1.0_wp + x*(-3.0_wp + 2.0_wp*x))
!
!  end function H5_P02
!
!
!  function H5_P1 ( x )
!
!    implicit none
!
!    real(wp) :: H5_P1
!    real(wp), intent(in) :: x
!
!    ! P1 = -(-1+x)^3*x*(1+3*x)
!    H5_P1 = -(-1.0_wp + x)**3*x*(1.0_wp + 3.0_wp*x)
!
!  end function H5_P1
!
!
!  function H5_P11 ( x )
!
!    implicit none
!
!    real(wp) :: H5_P11
!    real(wp), intent(in) :: x
!
!    ! P11 = -(-1+x)^2*(-1-2*x+15*x^2)
!    H5_P11 = -(-1.0_wp + x)**2*(-1.0_wp + x*(-2.0_wp + 15.0_wp*x))
!
!  end function H5_P11
!
!
!  function H5_P12 ( x )
!
!    implicit none
!
!    real(wp) :: H5_P12
!    real(wp), intent(in) :: x
!
!    ! P12 = -12*x*(3-8*x+5*x^2)
!    H5_P12 = -12.0_wp*x*(3.0_wp + x*(-8.0_wp + 5.0_wp*x**2))
!
!  end function H5_P12
!
!
!  function H5_P2 ( x )
!
!    implicit none
!
!    real(wp) :: H5_P2
!    real(wp), intent(in) :: x
!
!    ! P2 = -1/2*(-1+x)^3*x^2
!    H5_P2 = -0.5_wp*(-1.0_wp + x)**3*x**2
!
!  end function H5_P2
!
!
!  function H5_P21 ( x )
!
!    implicit none
!
!    real(wp) :: H5_P21
!    real(wp), intent(in) :: x
!
!    ! P21 = -1/2*(-1+x)**2*x*(-2+5*x)
!    H5_P21 = -0.5_wp*(-1.0_wp + x)**2*x*(-2.0_wp+5.0_wp*x)
!
!  end function H5_P21
!
!
!  function H5_P22 ( x )
!
!    implicit none
!
!    real(wp) :: H5_P22
!    real(wp), intent(in) :: x
!
!    ! P22 = 1-9*x+18*x^2-10*x^3
!    H5_P22 = 1.0_wp + x*(-9.0_wp + x*(18.0_wp - 10.0_wp*x))
!
!  end function H5_P22

end module acceleration_history
