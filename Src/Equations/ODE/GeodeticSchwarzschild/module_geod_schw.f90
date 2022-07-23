module geodesic_schwarzschild
!! Module that defines an equation class for evolving the geodesic equations
!! with forcing in a Schwarzschild spacetime.
!!
!! The implementation is found in [[submodule_geod_schw_implementation.f90]].

  use kinds
  use ode_equations

  implicit none

  type, extends(ode_equation) :: geod_schw
  !! A class that defines the geodesic evolution equations for a particle
  !! moving in a Schwarzschild spacetime.
  !!
  !! For this system [[geod_schw:nvars]]=7. The quantities stored in
  !! [[geod_schw:var_data]] are 1: \(r\), 2: \(\phi\), 3: \(u^t\),
  !! 4: \(u^r\), 5: \(u^{\phi}\), 6: \(m_q\), 7: \(\chi\), where
  !! \(r\) is the radial coordinate in Schwarzschild coordinates, \(\phi\)
  !! is the azimuthal angle,
  !! \(u^t, u^r\) and \(u^{\phi}\) are the time, radial and \(\phi\)
  !! components of the 4-velocity, \(m_q\) is the mass of the scalar charge
  !! and \(\chi\) s a variable defined in the osculating orbits
  !! framework that varies from \(0\) to \(2\pi\) over a full radial cycle.
  !! Note that \(m_q\) is stricly constant for all other cases than a
  !! scalar charge and should stricly not be included in a generic geodesic
  !! evolution. Also \(\chi\) is included for convenience as it is useful to
  !! control the turn on of the back-reaction in terms of timescales related
  !! to the radial cycle.
    real(wp) :: En
    !! The orbital energy per unit mass.
    real(wp) :: Lz
    !! The orbital angular momentum per unit mass.
    real(wp) :: e
    !! The orbital eccentricity.
    real(wp) :: p
    !! The orbital semi-latus rectum.
    real(wp) :: w
    !! The value of \(\chi\) at periapsis.
    real(wp) :: d2rdt2
    !! The second time derivative of the radial coordinate. Needed for the
    !! time dependent coordinate transformation.
    real(wp), dimension(4) :: force
    !! The self-force or any external force \(f_{\mu}\).
    real(wp), dimension(4) :: accel
    !! The acceleration \(a^{\mu}=
    !!   \frac{q}{m_q}(g^{\mu\nu}+u^{\mu}u^{\nu})f_{\nu}\).
    real(wp) :: udota
    !! The dot-product of the 4-velocity and the force \(q\, u^{\mu}f_{\mu}\).
  contains
    procedure :: init => geod_schw_init
    !! The [[equation:init]] routine is provided by [[geod_schw_init]].
    procedure :: rhs => geod_schw_rhs
    !! The [[equation:rhs]] routine is provided by [[geod_schw_rhs]].
    procedure :: output => geod_schw_output
    !! The [[equation:output]] routine is provided by [[geod_schw_output]].
    procedure :: save_globals_1 => geod_schw_save_globals_1
    !! The [[equation:save_globals_1]] routine is provided by
    !! [[geod_schw_save_globals_1]].
    procedure :: save_globals_2 => geod_schw_save_globals_2
    !! The [[equation:save_globals_2]] routine is provided by
    !! [[geod_schw_save_globals_2]].
    procedure :: load_globals => geod_schw_load_globals
    !! The [[equation:load_globals]] routine is provided by
    !! [[geod_schw_load_globals]].
!    procedure :: extract_time_derivatives => geod_schw_extract_td
!    !! Provide a routine [[geod_schw:extract_time_derivatives]] used
!    !! for extracting the time derivatives of the 4-acceleration.
    final :: close_geod_schw
    !! The finalizer.
  end type geod_schw

  interface
    module subroutine geod_schw_init ( this )
    !! The interface for the [[geod_schw]] version of [[equation:init]].
    !! This interface is consistent with [[eq_init_interface]].
      class(geod_schw), target, intent(inout) :: this
      !! The equation that is being initialized.
    end subroutine

    module subroutine geod_schw_rhs ( this )
    !! The interface for the [[geod_schw]] version of [[equation:rhs]]. This
    !! interface is consistent with [[eq_rhs_interface]].
      class(geod_schw), intent(inout) :: this
      !! The equation for which the RHS is calculated.
    end subroutine geod_schw_rhs

    module subroutine geod_schw_save_globals_1 ( this )
    !! The interface for the [[geod_schw]] version of
    !! [[equation:save_globals_1]]. This interface is consistent with
    !! [[eq_save_globals_1]].
      class(geod_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine geod_schw_save_globals_1

    module subroutine geod_schw_save_globals_2 ( this )
    !! The interface for the [[geod_schw]] version of
    !! [[equation:save_globals_2]]. This interface is consistent with
    !! [[eq_save_globals_2]].
      class(geod_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine geod_schw_save_globals_2

    module subroutine geod_schw_load_globals ( this )
    !! The interface for the [[geod_schw]] version of
    !! [[equation:load_globals]]. This interface is consistent with
    !! [[eq_load_globals]].
      class(geod_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine geod_schw_load_globals

    module subroutine geod_schw_output ( this )
    !! The interface for the [[geod_schw]] version of [[equation:output]].
    !! This interface is consistent with [[eq_output]].
      class(geod_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine geod_schw_output

    module subroutine close_geod_schw ( this )
    !! The interface for the [[close_geod_schw]] routine that provides a
    !! finalizer.
      type(geod_schw) :: this
      !! The equation that is being finalized.
    end subroutine close_geod_schw

  end interface

end module geodesic_schwarzschild
