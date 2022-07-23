submodule(osculating_schwarzschild) time_derivative_implementation
!! Routines to calculate the time derivative of the 4-acceleration based
!! on the time derivatives of the osculating orbit variables.
!!
!! Currently requires the fifth order Adams-Bashford-Moulton time integrator.

  use iso_c_binding

  implicit none

  interface
    function comp_d2edt2 ( alpha, beta, dalphadt, dbetadt, &
                           d2alphadt2, d2betadt2 ) bind(c, &
                         name='compute_d2edt2')
    !! Fortran interface to the C++ routine 'compute_d2edt2' that computes
    !! the second time derivative of the eccentricity from the time
    !! derivatives of \(\alpha\) and \(\beta\).

      use iso_c_binding
      real(c_double), intent(in), value :: alpha
      !! The current value of \(\alpha\).
      real(c_double), intent(in), value :: beta
      !! The current value of \(\beta\).
      real(c_double), intent(in), value :: dalphadt
      !! The current value of \(\dot{\alpha}\).
      real(c_double), intent(in), value :: dbetadt
      !! The current value of \(\dot{\beta}\).
      real(c_double), intent(in), value :: d2alphadt2
      !! The current value of \(\ddot{\alpha}\).
      real(c_double), intent(in), value :: d2betadt2
      !! The current value of \(\ddot{\beta}\).
      real(c_double) :: comp_d2edt2
      !! The return value is \(\ddot{e}\).
    end function comp_d2edt2

    function comp_d2wdt2 ( alpha, beta, dalphadt, dbetadt, &
                           d2alphadt2, d2betadt2 ) bind(c, &
                         name='compute_d2wdt2')
    !! Fortran interface to the C++ routine 'compute_d2wdt2' that computes
    !! the second time derivative of \(w\) from the time
    !! derivatives of \(\alpha\) and \(\beta\).

      use iso_c_binding
      real(c_double), intent(in), value :: alpha
      !! The current value of \(\alpha\).
      real(c_double), intent(in), value :: beta
      !! The current value of \(\beta\).
      real(c_double), intent(in), value :: dalphadt
      !! The current value of \(\dot{\alpha}\).
      real(c_double), intent(in), value :: dbetadt
      !! The current value of \(\dot{\beta}\).
      real(c_double), intent(in), value :: d2alphadt2
      !! The current value of \(\ddot{\alpha}\).
      real(c_double), intent(in), value :: d2betadt2
      !! The current value of \(\ddot{\beta}\).
      real(c_double) :: comp_d2wdt2
      !! The return value is \(\ddot{w}\).
    end function comp_d2wdt2

    function comp_dedt ( alpha, beta, dalphadt, dbetadt ) bind(c, &
                         name='compute_dedt')
    !! Fortran interface to the C++ routine 'compute_dedt' that computes
    !! the time derivative of the eccentricity from the time
    !! derivatives of \(\alpha\) and \(\beta\).

      use iso_c_binding
      real(c_double), intent(in), value :: alpha
      !! The current value of \(\alpha\).
      real(c_double), intent(in), value :: beta
      !! The current value of \(\beta\).
      real(c_double), intent(in), value :: dalphadt
      !! The current value of \(\dot{\alpha}\).
      real(c_double), intent(in), value :: dbetadt
      !! The current value of \(\dot{\beta}\).
      real(c_double) :: comp_dedt
      !! The return value is \(\dot{e}\).
    end function comp_dedt

    function comp_dwdt ( alpha, beta, dalphadt, dbetadt ) bind(c, &
                         name='compute_dwdt')
    !! Fortran interface to the C++ routine 'compute_dedt' that computes
    !! the time derivative of \(w\) from the time
    !! derivatives of \(\alpha\) and \(\beta\).

      use iso_c_binding
      real(c_double), intent(in), value :: alpha
      !! The current value of \(\alpha\).
      real(c_double), intent(in), value :: beta
      !! The current value of \(\beta\).
      real(c_double), intent(in), value :: dalphadt
      !! The current value of \(\dot{\alpha}\).
      real(c_double), intent(in), value :: dbetadt
      !! The current value of \(\dot{\beta}\).
      real(c_double) :: comp_dwdt
      !! The return value is \(\dot{w}\).
    end function comp_dwdt

    function comp_w ( alpha, beta ) bind(c, name='compute_w')
    !! Fortran interface to the C++ routine 'comp_w' that computes
    !! \(w\) from \(\alpha\) and \(\beta\).

      use iso_c_binding
      real(c_double), intent(in), value :: alpha
      !! The current value of \(\alpha\).
      real(c_double), intent(in), value :: beta
      !! The current value of \(\beta\).
      real(c_double) :: comp_w
      !! The return value is \(w\).
    end function comp_w

    function comp_v ( chi, w ) bind(c, name='compute_v')
    !! Fortran interface to the C++ routine 'comp_v' that computes
    !! \(v\) from \(\chi\) and \(w\).

      use iso_c_binding
      real(c_double), intent(in), value :: chi
      !! The current value of \(\chi\).
      real(c_double), intent(in), value :: w
      !! The current value of \(w\).
      real(c_double) :: comp_v
      !! The return value is \(v\).
    end function comp_v

    function comp_dardt ( d2edt2, d2wdt2, d2chidt2, dpdt, dedt, dwdt, dchidt, &
                          e, p, v, mass ) bind(c, name='compute_dardt')
    !! Fortran interface to the C++ routine 'comp_dardt' that computes
    !! \(\dot{a}^r\) from \(\ddot{e}\), \(\ddot{w}\), \(\ddot{\chi}\),
    !! \(\dot{p}\), \(\dot{e}\), \(\dot{w}\), \(\dot{\chi}\), \(e\), \(p\),
    !! \(v\) and \(M\).

      use iso_c_binding
      real(c_double), intent(in), value :: d2edt2
      !! The current value of \(\ddot{e}\).
      real(c_double), intent(in), value :: d2wdt2
      !! The current value of \(\ddot{w}\).
      real(c_double), intent(in), value :: d2chidt2
      !! The current value of \(\ddot{\chi}\).
      real(c_double), intent(in), value :: dpdt
      !! The current value of \(\dot{p}\).
      real(c_double), intent(in), value :: dedt
      !! The current value of \(\dot{e}\).
      real(c_double), intent(in), value :: dwdt
      !! The current value of \(\dot{w}\).
      real(c_double), intent(in), value :: dchidt
      !! The current value of \(\dot{\chi}\).
      real(c_double), intent(in), value :: e
      !! The current value of \(e\).
      real(c_double), intent(in), value :: p
      !! The current value of \(p\).
      real(c_double), intent(in), value :: v
      !! The current value of \(v\).
      real(c_double), intent(in), value :: mass
      !! The black hole mass.
      real(c_double) :: comp_dardt
      !! The return value is \(\dot{a}^r\).
    end function comp_dardt

    function comp_daphidt ( d2edt2, d2wdt2, d2chidt2, dpdt, dedt, dwdt, &
                            dchidt, e, p, v, mass ) bind(c, &
                            name='compute_daphidt')
    !! Fortran interface to the C++ routine 'comp_daphidt' that computes
    !! \(\dot{a}^{\phi}\) from \(\ddot{e}\), \(\ddot{w}\), \(\ddot{\chi}\),
    !! \(\dot{p}\), \(\dot{e}\), \(\dot{w}\), \(\dot{\chi}\), \(e\), \(p\),
    !! \(v\) and \(M\).

      use iso_c_binding
      real(c_double), intent(in), value :: d2edt2
      !! The current value of \(\ddot{e}\).
      real(c_double), intent(in), value :: d2wdt2
      !! The current value of \(\ddot{w}\).
      real(c_double), intent(in), value :: d2chidt2
      !! The current value of \(\ddot{\chi}\).
      real(c_double), intent(in), value :: dpdt
      !! The current value of \(\dot{p}\).
      real(c_double), intent(in), value :: dedt
      !! The current value of \(\dot{e}\).
      real(c_double), intent(in), value :: dwdt
      !! The current value of \(\dot{w}\).
      real(c_double), intent(in), value :: dchidt
      !! The current value of \(\dot{\chi}\).
      real(c_double), intent(in), value :: e
      !! The current value of \(e\).
      real(c_double), intent(in), value :: p
      !! The current value of \(p\).
      real(c_double), intent(in), value :: v
      !! The current value of \(v\).
      real(c_double), intent(in), value :: mass
      !! The black hole mass.
      real(c_double) :: comp_daphidt
      !! The return value is \(\dot{a}^{\phi}\).
    end function comp_daphidt

  end interface

contains

  module procedure osc_schw_extract_td

    use all_integrators, only : mol_int
    use parameters, only : mass

    implicit none

    real(wp), dimension(2) :: dchi, dp, dalpha, dbeta
    real(wp) :: d2edt2, d2wdt2, dedt, dwdt, w, v, dardt, daphidt

    if (mol_int%iname == 'abmv5') then
      dchi(1:2) = extract_time_derivatives ( this%tmp_data(:,1), 2 )
      dp(1:1) = extract_time_derivatives ( this%tmp_data(:,3), 1 )
      dalpha(1:2) = extract_time_derivatives ( this%tmp_data(:,4), 2 )
      dbeta(1:2) = extract_time_derivatives ( this%tmp_data(:,5), 2 )

      associate ( alpha => this%var_data(4), &
                  beta => this%var_data(5), &
                  dalphadt => dalpha(1), &
                  dbetadt => dbeta(1), &
                  d2alphadt2 => dalpha(2), &
                  d2betadt2 => dbeta(2), &
                  chi => this%var_data(1), &
                  dchidt => dchi(1), &
                  d2chidt2 => dchi(2), &
                  dpdt => dp(1), &
                  e => this%e, &
                  p => this%var_data(3) )

        d2edt2 = comp_d2edt2 ( alpha, beta, dalphadt, &
                               dbetadt, d2alphadt2, d2betadt2 )
        d2wdt2 = comp_d2wdt2 ( alpha, beta, dalphadt, &
                               dbetadt, d2alphadt2, d2betadt2 )
        dedt = comp_dedt ( alpha, beta, dalphadt, dbetadt )
        dwdt = comp_dwdt ( alpha, beta, dalphadt, dbetadt )
        w = comp_w ( alpha, beta )
        v = comp_v ( chi, w )
        dardt = comp_dardt ( d2edt2, d2wdt2, d2chidt2, dpdt, dedt, dwdt, &
                             dchidt, e, p, v, mass )
        daphidt = comp_daphidt ( d2edt2, d2wdt2, d2chidt2, dpdt, dedt, dwdt, &
                                 dchidt, e, p, v, mass )
      end associate
      osc_schw_extract_td = (/ 0.0_wp, dardt, 0.0_wp, daphidt /)
    else
      print*,'Extraction of time derivatives of 4-acceleration only &
             &supported when the time integrator is abmv5'
    end if
  end procedure osc_schw_extract_td


  function extract_time_derivatives ( tmpdata, order )
  !! Function that extracts the time derivatives of evolution variables when
  !! the fifth order Adams-Bashford-Moulton integrator is used.

    use numerics, only : factorial
    use time_info, only : get_current_dtime

    implicit none

    real(wp), dimension(:), intent(in) :: tmpdata
    !! The ABMV5 temporary date for an evolved variable.
    integer(ip), intent(in) :: order
    !! The order of the time derivative to extract.
    real(wp), dimension(order) :: extract_time_derivatives
    integer(ip) :: i
    real(wp) :: dtime

    if ( order>5 .or. order<0 ) then
      print*,'extract_scalar_time_derivatives called with illegal order'
      stop
    end if

    dtime = get_current_dtime ( )

    do i = 1, order
      extract_time_derivatives(i) = tmpdata(i+1)*factorial(order)/dtime**i
    end do
  end function extract_time_derivatives


end submodule time_derivative_implementation
