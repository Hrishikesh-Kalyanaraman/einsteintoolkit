module numerics
!! Module that contains a number of more or less useful numerical routines.

  use kinds

  implicit none

  real(sp), parameter :: epss = 1.0e-5_sp
  !! A single precision constant with a value used to determine convergence
  !! of iterative methods.
  real(dp), parameter :: epsd = 1.0e-12_dp
  !! A double precision constant with a value used to determine convergence
  !! of iterative methods.
  real(qp), parameter :: epsq = 1.0e-32_qp
  !! A quad precision constant with a value used to determine convergence
  !! of iterative methods.
  integer(ip), parameter :: maxiter = 100
  !! An integer constant with a maximum number of iterations for iterative
  !! methods.

  interface eps
    module procedure eps_prec_s, eps_prec_d, eps_prec_q
    !! The generic name for a functions that returns the convergence criterium
    !! depending on the floating point data type.
  end interface eps

contains

  function eps_prec_s (real_var)
  !! Single precision version of a function that returns the convergence
  !! criterium.
    real(sp) :: real_var
    !! Any single precision variable.
    real(sp) :: eps_prec_s
    !! The single precision convergence criterium.

    eps_prec_s = epss
  end function eps_prec_s


  function eps_prec_d (real_var)
  !! Double precision version of a function that returns the convergence
  !! criterium.
    real(dp) :: real_var
    !! Any double precision variable.
    real(dp) :: eps_prec_d
    !! The double precision convergence criterium.

    eps_prec_d = epsd
  end function eps_prec_d


  function eps_prec_q (real_var)
  !! Quad precision version of a function that returns the convergence
  !! criterium.
    real(qp) :: real_var
    !! Any quad precision variable.
    real(qp) :: eps_prec_q
    !! The quad precision convergence criterium.

    eps_prec_q = epsq
  end function eps_prec_q


  function Lambert ( z )
  !! Function to calculate Lambert's W-function.
  !!
  !! Algorithm found at
  !! http://en.citizendium.org/wiki/Lambert_W_function#Numerical_calculation

    implicit none

    real(wp), intent(in) :: z
    !! The input value.
    real(wp) :: Lambert
    !! The return value Lambert(z).
    real(wp) :: wcurrent, wnew, expw, diff
    integer(ip) :: iter
!    real(wp), parameter :: eps = 1.0e-12_wp

    iter = 0
    wcurrent = 1.0_wp

    loop: do
      iter = iter+1
      expw = exp(wcurrent)
      diff = wcurrent*expw - z
      wnew = wcurrent - diff / ( expw*(wcurrent+1) - &
             ((wcurrent+2)*diff/(2*wcurrent+2)))
      if ( abs(wnew - wcurrent) < eps(1.0_wp) ) exit loop
      if (iter>maxiter) stop 'Function Lambert failed to converge. Aborting'
      wcurrent = wnew
    end do loop
    Lambert = wnew
  end function Lambert


  function rschw( z, mass )
  !! Function to invert the tortoise radius as a function of Schwarzschild
  !! radius.
  !!
  !! Simple Newton root finding algorithm. Have problems converging for
  !! negative rstar as it overshoots and rschw-2M can become negative.

    implicit none

    real(wp), intent(in) :: z
    !! The tortoise radius to invert, \(z\).
    real(wp), intent(in) :: mass
    !! The mass of the black hole, \(M\).
    real(wp) :: rschw
    !! The return value is \(r_{\mathrm{schw}}(z)-2M\).
    real(wp) :: xcurrent, xnew
    integer(ip) :: iter
!   real(wp), parameter :: eps = 1.0e-12_wp

    iter = 0
    xcurrent = 1.0_wp

    loop: do
      iter = iter+1
      xnew = (z - 2.0_wp*mass*log(0.5_wp*xcurrent/mass)) &
             /(1+2.0_wp*mass/xcurrent)
      if ( abs ( xnew - xcurrent) < eps(1.0_wp) ) exit loop
      if (iter>maxiter) stop 'Function rschw failed to converge. Aborting'
      xcurrent = xnew
    end do loop
    rschw = xnew
  end function rschw

 
  function invert_tortoise ( rstar, mass )
  !! Function to invert the tortoise radius as a function of Schwarzschild
  !! radius.
  !!
  !! Uses rschw for rstar>=0 and Lambert for rstar<0. The Lmabert method has
  !! problems for large rstar due to numerical overflow of it's argument.

    implicit none

    real(wp), intent(in) :: rstar
    !! The tortoise radius to invert, \(r_*\).
    real(wp), intent(in) :: mass
    !! The mass of the black hole, \(M\).
    real(wp) :: invert_tortoise
    !! The return value is \(r_{\mathrm{schw}}(r_*) - 2M\).


    if ( rstar >= 0 ) then
      invert_tortoise = rschw ( rstar, mass )
    else
      invert_tortoise = 2.0_wp *mass * Lambert(exp(0.5_wp*rstar/mass - 1.0_wp))
    end if
  end function invert_tortoise


  function rstar_of_r ( r, mass )
  !! Function to calculate the tortoise radius, \(r_*\), as function of the
  !! Schwarzschild radius, \(r_{\mathrm{schw}}\).

    implicit none

    real(wp), intent(in) :: r
    !! The schwarzschild radius, \(r_{\mathrm{schw}}\).
    real(wp), intent(in) :: mass
    !! The mass of the black hole, \(M\).
    real(wp) :: rstar_of_r
    !! The return value is \(r_* = r_{\mathrm{schw}}+2 M\log\left (
    !! \frac{r_{\mathrm{schw}}}{2 M} - 1\right )\).

    rstar_of_r = r+2.0_wp*mass*log(r/(2.0_wp*mass)-1.0_wp)
  end function rstar_of_r


  subroutine transition ( rho, R, S, fT, fTp, fTpp )
  !! Routine to calculate the smooth transition function, \(fT\), and it's first
  !! \(fT'\) and second \(fT''\) derivative with respect to the computational
  !! coordinate, \(\rho\).

    implicit none

    real(wp), intent(in) :: rho
    !! The computational coordinate, \(\rho\), where the transition function
    !! should be calculated. Should be between \(R\) and \(S\).
    real(wp), intent(in) :: R
    !! The value of \(\rho\) where \(fT\) should be \(0\).
    real(wp), intent(in) :: S
    !! The value of \(\rho\) where \(fT\) should be \(1\).
    real(wp), intent(out) :: fT
    !! On return contains \(fT\).
    real(wp), intent(out) :: fTp
    !! On return contains \(fT'\).
    real(wp), intent(out) :: fTpp
    !! On return contains \(fT''\).
    real(wp), parameter :: es = 1.5_wp
    real(wp), parameter :: q2 = 1.3_wp
    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp), parameter :: half = 0.5_wp
    real(wp) :: x, tanx, cotx, fac, tanhfac, cscx2, secx2, sechfac2

    x = max ( half * pi * ( rho - R ) / ( S - R ), 0.0_wp )
    if ( x > 0.0_wp ) then
!      print*,'x = ', x
      tanx = tan ( x )
!      print*,'tanx = ', tanx
      cotx = 1.0_wp / tanx
!      print*,'cotx = ', cotx
      fac = es/pi * ( tanx - q2 * cotx )
!      print*,'fac = ', fac
      tanhfac = tanh ( fac)
!      print*,'tanhfac = ', tanhfac
      fT = half + half * tanhfac
!      print*,'fT = ', fT
      cscx2 = 1.0_wp / sin ( x )**2
      secx2 = 1.0_wp / cos ( x )**2
      sechfac2 = 1.0_wp / cosh ( fac )**2
      fTp = 0.25_wp * es * ( q2 * cscx2 + secx2 ) * sechfac2 / ( S - R )
      fTpp = -0.25_wp * es * sechfac2 * ( &
               pi * ( q2 * cotx * cscx2 - secx2 * tanx ) + &
               es * ( q2 * cscx2 + secx2 )**2 * tanhfac ) / ( S - R )**2
    else
      fT = 0.0_wp
      fTp = 0.0_wp
      fTpp = 0.0_wp
    end if
  end subroutine transition


  subroutine time_window ( time, tsigma, norder, tfac, dtfac_dt, d2tfac_dt2 )
  !! Routine to calculate a smooth "Gaussian" type time window function,
  !! \(f_t\) and it's first \(f_t'\) and second \(f_t''\) time derivative.

    implicit none

    real(wp), intent(in) :: time
    !! The time at which to calculate the time window function, \(t\).
    real(wp), intent(in) :: tsigma
    !! The width of the time window function, \(\sigma_t\).
    integer(ip), intent(in) :: norder
    !! The order of the time window function, \(n\).
    real(wp), intent(out) :: tfac
    !! On output contains \(f_t = 1-e^{\left (\frac{t}{\sigma_t}\right )^n}.\)
    real(wp), intent(out) :: dtfac_dt
    !! On output contains \(f_t'\)
    real(wp), intent(out) :: d2tfac_dt2
    !! On output contains \(f_t''\)
    real(wp) :: tfactor, expfactor

    tfactor = (time/tsigma)**norder
    expfactor = exp(-tfactor)
    tfac = 1.0_wp - expfactor
    dtfac_dt = norder*time**(norder-1)/tsigma**norder*expfactor
    d2tfac_dt2 = -norder*(1.0_wp+norder*(tfactor-1.0_wp)) &
                  *time**(norder-2)/tsigma**norder*expfactor
  end subroutine time_window


  function ldep ( l, order )
  !! Function that calculates the higher order terms, \(c_n(\ell)\) in the
  !! self-force expansion over \(\ell\).

    integer(ip), intent(in) :: l
    !! The \(\ell\)-value.
    integer(ip), intent(in) :: order
    !! The order of the higher order term, \(n\).
    real(wp) :: ldep
    !! The return value is \[c_n(\ell)=\begin{cases}
    !!   \frac{1}{(2\ell-1)(2\ell+3)}, & n=2 \\
    !!   \frac{1}{(2\ell-3)(2\ell-1)(2\ell+3)(2\ell+5)}, & n=3 \\
    !!   \frac{1}{(2\ell-5)(2\ell-3)(2\ell-1)(2\ell+3)(2\ell+5)(2\ell+7)}, & n=4
    !! \end{cases} \].

    real(wp) :: lm

    lm = real(l,wp)
    select case (order)
    case (2)
      ldep = 1.0_wp/((2*lm-1)*(2*lm+3))
    case (4)
      ldep = 1.0_wp/((2*lm-3)*(2*lm-1)*(2*lm+3)*(2*lm+5))
    case (6)
      ldep = 1.0_wp/((2*lm-5)*(2*lm-3)*(2*lm-1)*(2*lm+3)*(2*lm+5)*(2*lm+7))
    case default
      stop "error: ldep called with unsupported order"
    end select
  end function ldep


  function lsum ( lmin, order )
  !! Function that calculates the sum of higher order terms, \(c_n(\ell)\),
  !! provided by [[ldep]] from \(\ell_{\mathrm{min}}\) to \(\infty\).

    integer(ip), intent(in) :: lmin
    !! The value of \(\ell\) at which to start the sum,
    !! \(\ell_{\mathrm{min}}\).
    integer(ip), intent(in) :: order
    !! The order of the higher order term, \(n\).
    real(wp) :: lsum
    !! The return value is \[\sum_{\ell=\ell_{\mathrm{min}}}^{\infty} c_n(\ell) =
    !!   \begin{cases}
    !!     \frac{\ell_{\mathrm{min}}}{4\ell_{\mathrm{min}}^2-1}, & n=2 \\
    !!     \frac{\ell_{\mathrm{min}}}{48\ell_{\mathrm{min}}^4-
    !!       120\ell_{\mathrm{min}}^2+27}, & n=3 \\
    !!     \frac{\ell_{\mathrm{min}}}{5 (64\ell_{\mathrm{min}}^6-
    !!       560\ell_{\mathrm{min}}^4+1036\ell_{\mathrm{min}}^2-225)}, & n=4
    !!   \end{cases} \].
    real(wp) :: lm

    lm = real(lmin,wp)
    select case (order)
    case (2)
      lsum = lm/(4.0_wp*lm**2-1.0_wp)
    case (4)
      lsum = lm/(48.0_wp*lm**4-120.0_wp*lm**2+27.0_wp)
    case (6)
      lsum = lm/(5.0_wp &
                 *(64.0_wp*lm**6-560.0_wp*lm**4+1036.0_wp*lm**2-225.0_wp))
    case default
      stop "error: lsum called with unsupported order"
    end select

  end function lsum


  function correct_for_higher_modes ( ydat, nmodes, &
                                      startfit, endfit, order, npar )
  !! Function that fits higher order terms to data containing amplitude as
  !! function of l-values and corrects the sum over the evolved l-modes with
  !! the contribution from the not evolved l-modes.

    use gsl_interface

    real(wp), dimension(:), intent(in) :: ydat
    !! A 1d-array containing the amplitudes of the \(\ell\)-modes. It is
    !! assumed, for now, that the data starts at \(\ell=0\).
    integer(ip), intent(in) :: nmodes
    !! The number of \(\ell\)-modes.
    integer(ip), intent(in) :: startfit
    !! At what \(\ell\)-value should the fit begin.
    integer(ip), intent(in) :: endfit
    !! At what \(\ell\)-value should the fit end.
    integer(ip), intent(in) :: order
    !! What is the lowest order term that should be fitted.
    integer(ip), intent(in) :: npar
    !! How many terms should be included in the correction. One more term
    !! is included in the fit, but the coefficient for the last fit
    !! coefficient is only used as an error estimator.
    real(wp), dimension(3) :: correct_for_higher_modes
    !! A 1d-array of size 3 that on return contains the sum of evolved modes
    !! (first entry), the corrected sum (second entry) and the error
    !! estimate (third entry).
    integer(ip) :: n
    real(wp), dimension(npar+1,endfit-startfit+1) :: xm
    real(wp), dimension(endfit-startfit+1) :: yv
    real(wp), dimension(npar+1) :: cv
    real(wp), dimension(npar+1,npar+1) :: covm
    real(wp) :: chisq
    integer(ip) :: ierr, i, j, l
    real(wp) :: c1, c2, c3

    n = endfit-startfit+1
    do i = 1, n
      l = startfit+i-1
      yv(i) = ydat(startfit+i)
      do j = 1, npar+1
        xm(j,i) = ldep(l,order+2*(j-1))
      end do
    end do

    call multifit_linear ( n, npar+1, reshape(xm,(/n*(npar+1)/)), yv, cv, &
                           reshape(covm,(/(npar+1)*(npar+1)/)), chisq, ierr )

    correct_for_higher_modes(1) = sum(ydat)
    correct_for_higher_modes(2) = correct_for_higher_modes(1)
    do j = 1, npar
      correct_for_higher_modes(2) = correct_for_higher_modes(2) + &
                                   cv(j)*lsum(nmodes,order+2*(j-1))
    end do
    correct_for_higher_modes(3) = cv(npar+1)*lsum(nmodes,order+2*npar)
  end function correct_for_higher_modes


  recursive function factorial ( n ) result ( fac )
  !! A simple factorial function. Use only for small values of \(n\) as
  !! no consideration of efficiency has been made.

    implicit none

    integer(ip), intent(in) :: n
    !! The value for which the factorial function should be calculated.
    integer(ip) :: fac
    !! The return value, \(n!\).

    if (n<2) then
      fac = 1
    else
      fac = n*factorial(n-1)
    end if

  end function factorial

end module numerics
