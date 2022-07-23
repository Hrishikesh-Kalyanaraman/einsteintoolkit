module osculating_schwarzschild
!! Module that defines an equation class for evolving the geodesic equations
!! using the osculating orbits framework with forcing in a Schwarzschild
!! spacetime.
!!
!! The implementation is found in [[submodule_osc_schw_implementation.f90]].

  use kinds
  use ode_equations

  implicit none

  type, extends(ode_equation) :: osc_schw
  !! A class that defines the osculating orbit evolution equations for a
  !! particle !! moving in a Schwarzschild spacetime.
  !!
  !! For this system [[osc_schw:nvars]]=6. The quantities stored in
  !! [[osc_schw:var_data]] are 1: \(\chi\), 2: \(\phi\), 3: \(p\),
  !! 4: \(\alpha\), 5: \(\beta\), 6: \(m_q\), where \(\chi\) is a variable
  !! that varies from \(0\) to \(2\pi\) over a full radial cycle, \(\phi\) is
  !! the azimuthal angle, \(p\) is the semi-latus rectum, \(\alpha=e\sin(w)\),
  !! \(\beta=e\cos(w)\) and \(m_q\) is the mass of the scalar charge. The
  !! quantities \(e\) and \(w\) are defined below. Note that \(m_q\) is
  !! stricly constant for all other cases than a scalar charge
  !! and should stricly not be included in a generic geodesic evolution.
    real(wp) :: En
    !! The orbital energy per unit mass.
    real(wp) :: Lz
    !! The orbital angular momentum per unit mass.
    real(wp) :: e
    !! The orbital eccentricity, \(e=\sqrt{\alpha^2+\beta^2}\).
    real(wp) :: r
    !! The radial coordinate in Schwarzschild coordinates,
    !! \(r=\frac{p M}{1+e\cos(\chi-w)}\).
    real(wp) :: d2chidt2
    !! The second time derivative of \(\chi\).
    real(wp) :: drdt
    !! The time derivative of \(r\).
    real(wp) :: d2rdt2
    !! The second time derivative of \(r\).
    real(wp) :: w
    !! The value of \(\chi\) at periapsis.
    real(wp) :: mass
    !! Local copy of the run time parameter setting the black hole mass, \(M\).
    real(wp) :: ur
    !! The radial component of the 4-velocity, \(u^r\). Needed for the
    !! effective source.
    real(wp), dimension(4) :: force
     !! The self-force or any external force \(f_{\mu}\).
    real(wp), dimension(4) :: accel
    !! The acceleration \(a^{\mu}=
    !!   \frac{q}{m_q}(g^{\mu\nu}+u^{\mu}u^{\nu})f_{\nu}\).
    real(wp) :: udota
    !! The dot-product of the 4-velocity and the force \(q\, u^{\mu}f_{\mu}\).
  contains
    procedure :: init => osc_schw_init
    !! The [[equation:init]] routine is provided by [[osc_schw_init]].
    procedure :: rhs => osc_schw_rhs
    !! The [[equation:rhs]] routine is provided by [[osc_schw_rhs]].
    procedure :: output => osc_schw_output
    !! The [[equation:output]] routine is provided by [[osc_schw_output]].
    procedure :: save_globals_1 => osc_schw_save_globals_1
    !! The [[equation:save_globals_1]] routine is provided by
    !! [[geod_schw_save_globals_1]].
    procedure :: save_globals_2 => osc_schw_save_globals_2
    !! The [[equation:save_globals_2]] routine is provided by
    !! [[geod_schw_save_globals_2]].
    procedure :: load_globals => osc_schw_load_globals
    !! The [[equation:load_globals]] routine is provided by
    !! [[geod_schw_load_globals]].
    procedure :: calc_dependent
    !! A routine for calculating the dependent from the evolved variables.
    procedure :: extract_time_derivatives => osc_schw_extract_td
    !! A routine for extracting the time derivatives of the 4-acceleration.
    !! Currently only works when the Adams-Bashford-Moulton multi-value
    !! integrator is used.
    final :: close_osc_schw
    !! The finalizer.
  end type osc_schw

  interface
    module subroutine osc_schw_init ( this )
    !! The interface for the [[osc_schw]] version of [[equation:init]].
    !! This interface is consistent with [[eq_init_interface]].
      class(osc_schw), target, intent(inout) :: this
      !! The equation that is being initialized.
    end subroutine osc_schw_init

    module subroutine osc_schw_rhs ( this )
    !! The interface for the [[osc_schw]] version of [[equation:rhs]]. This
    !! interface is consistent with [[eq_rhs_interface]].
      class(osc_schw), intent(inout) :: this
      !! The equation for which the RHS is calculated.
    end subroutine osc_schw_rhs

    module subroutine osc_schw_output ( this )
    !! The interface for the [[osc_schw]] version of [[equation:output]].
    !! This interface is consistent with [[eq_output]].
      class(osc_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine osc_schw_output

    module subroutine osc_schw_save_globals_1 ( this )
    !! The interface for the [[osc_schw]] version of
    !! [[equation:save_globals_1]]. This interface is consistent with
    !! [[eq_save_globals_1]].
      class(osc_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine osc_schw_save_globals_1

    module subroutine osc_schw_save_globals_2 ( this )
    !! The interface for the [[osc_schw]] version of
    !! [[equation:save_globals_2]]. This interface is consistent with
    !! [[eq_save_globals_2]].
      class(osc_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine osc_schw_save_globals_2

    module subroutine osc_schw_load_globals ( this )
    !! The interface for the [[osc_schw]] version of
    !! [[equation:load_globals]]. This interface is consistent with
    !! [[eq_load_globals]].
      class(osc_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine osc_schw_load_globals

    module subroutine calc_dependent ( this ) 
    !! The interface for the [[calc_dependent]] routine that calculates
    ! dependent variables from the evolved variables.
      class(osc_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine calc_dependent

    module subroutine close_osc_schw ( this )
    !! The interface for the [[close_osc_schw]] routine that provides a
    !! finalizer.
      type(osc_schw) :: this
      !! The equation that is being finalized.
    end subroutine close_osc_schw

    module function osc_schw_extract_td ( this )
    !! The interface for the routine that extracts the time derivative
    !! of the 4-acceleration.
      class(osc_schw) :: this
      !! The routine is called on this equation.
      real(wp), dimension(4) :: osc_schw_extract_td
      !! The return value is a size (4) real array containing the coordinate
      !! time derivative of the 4-acceleration \[\dot{a}^{\mu}\].
    end function osc_schw_extract_td
  end interface

end module osculating_schwarzschild
