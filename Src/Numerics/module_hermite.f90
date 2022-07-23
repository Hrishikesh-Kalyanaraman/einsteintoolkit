module hermite

  use kinds

  interface
    module function H3_interpolate ( x, y, xv, dydx ) result ( res )
      real(wp), dimension(-3:1), intent(in) :: x
      real(wp), dimension(-3:1), intent(in) :: y
      real(wp), intent(in) :: xv
      real(wp), dimension(0:1), intent(in), optional :: dydx
      real(wp), dimension(3) :: res
    end function H3_interpolate

    module function H5_interpolate ( x, y, xv, dydx, d2ydx2 ) result ( res )
      real(wp), dimension(-5:1), intent(in) :: x
      real(wp), dimension(-5:1), intent(in) :: y
      real(wp), intent(in) :: xv
      real(wp), dimension(0:1), intent(in), optional :: dydx
      real(wp), dimension(0:1), intent(in), optional :: d2ydx2
      real(wp), dimension(3) :: res
    end function H5_interpolate
  end interface


end module hermite
