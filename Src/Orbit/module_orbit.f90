module orbit_base
!! Module that defines a global orbit object that can be used to pass
!! information between different equations using the save_globals/load_globals
!! mechanism.

  use kinds

  type :: g_orbit
  !! A global orbit type.
    real(wp), private :: r
    !! The radial coordinate, \(r\) (in Schwarzschild coordinates).
    real(wp), private :: phi
    !! The azimuthal angle, \(\phi\).
    real(wp), private :: ur
    !! The radial component of the 4-velocity, \(u^r\).
    real(wp), private :: En
    !! The energy per unit mass of the orbit, \(E\).
    real(wp), private :: Lz
    !! The angular momentum per unit mass of the orbit, \(L_z\).
    real(wp), private :: chi
    !! The osculaing orbits parameter, \(\chi\) that changes by \(2\pi\) over
    !! a full radial cycle.
  contains
    procedure :: set_orbit
    !! Routine to set the [[g_orbit]] values.
    procedure :: get_orbit
    !! Routine to get the [[g_orbit]] values.
    procedure :: get_chi
    !! Routine to get \(\chi\) alone.
  end type g_orbit

  type :: tdc_orbit
  !! A global type with information needed by the time dependent coordinate
  !! transformation object.
    real(wp), private :: r
    !! The radial coordinate, \(r\) (in Schwarzschild coordinates).
    real(wp), private :: drdt
    !! The time derivative of the radial coordinate, \(\dot{r}\).
    real(wp), private :: d2rdt2
    !! The second time derivative of the radial coordinate, \(\ddot{r}\).
  contains
    procedure :: set_tdc
    !! Routine to set the [[tdc_orbit]] values.
    procedure :: get_tdc
    !! Routine to get the [[tdc_orbit]] values.
  end type tdc_orbit

  type(g_orbit) :: orbit_info
  !! The global [[g_orbit]] object that is available by use association.
  type(tdc_orbit) :: tdc_info
  !! The global [[tdc_orbit]] object that is available by use association.

contains

  subroutine set_orbit ( this, r, phi, ur, En, Lz, chi )
  !! Routine that sets all the orbit variables.
    class(g_orbit), intent(inout) :: this
    !! The routine is called on this [[g_orbit]] object.
    real(wp), intent(in) :: r
    !! The radial coordinate, \(r\) (in Schwarzschild coordinates).
    real(wp), intent(in) :: phi
    !! The azimuthal angle, \(\phi\).
    real(wp), intent(in) :: ur
    !! The radial component of the 4-velocity, \(u^r\).
    real(wp), intent(in) :: En
    !! The energy per unit mass of the orbit, \(E\).
    real(wp), intent(in) :: Lz
    !! The angular momentum per unit mass of the orbit, \(L_z\).
    real(wp), intent(in) :: chi
    !! The osculaing orbits parameter, \(\chi\) that changes by \(2\pi\) over
    !! a full radial cycle.

    this%r = r
    this%phi = phi
    this%ur = ur
    this%En = En
    this%Lz = Lz
    this%chi = chi
  end subroutine set_orbit

  subroutine get_orbit( this, r, phi, ur, En, Lz )
  !! Routine that gets all the orbit variables except for \(\chi\).
    class(g_orbit), intent(inout) :: this
    !! The routine is called on this [[g_orbit]] object.
    real(wp), intent(out) :: r
    !! The radial coordinate, \(r\) (in Schwarzschild coordinates).
    real(wp), intent(out) :: phi
    !! The azimuthal angle, \(\phi\).
    real(wp), intent(out) :: ur
    !! The radial component of the 4-velocity, \(u^r\).
    real(wp), intent(out) :: En
    !! The energy per unit mass of the orbit, \(E\).
    real(wp), intent(out) :: Lz
    !! The angular momentum per unit mass of the orbit, \(L_z\).

    r = this%r
    phi = this%phi
    ur = this%ur 
    En = this%En
    Lz = this%Lz
  end subroutine get_orbit

  subroutine get_chi ( this, chi )
  !! Routine that gets \(\chi\).
    class(g_orbit), intent(inout) :: this
    !! The routine is called on this [[g_orbit]] object.
    real(wp), intent(out) :: chi
    !! The osculaing orbits parameter, \(\chi\) that changes by \(2\pi\) over
    !! a full radial cycle.

    chi = this%chi
  end subroutine get_chi

  subroutine set_tdc ( this, r, drdt, d2rdt2 )
  ! Routine that sets the information needed by the time dependent time object.
    class(tdc_orbit), intent(inout) :: this
    !! The routine is called on this [[tdc_orbit]] object.
    real(wp), intent(in) :: r
    !! The radial coordinate, \(r\) (in Schwarzschild coordinates).
    real(wp), intent(in) :: drdt
    !! The time derivative of the radial coordinate, \(\dot{r}\).
    real(wp), intent(in) :: d2rdt2
    !! The second time derivative of the radial coordinate, \(\ddot{r}\).

    this%r = r
    this%drdt = drdt
    this%d2rdt2 = d2rdt2
  end subroutine set_tdc

  subroutine get_tdc ( this, r, drdt, d2rdt2 )
  ! Routine that gets the information needed by the time dependent time object.
    class(tdc_orbit), intent(inout) :: this
    !! The routine is called on this [[tdc_orbit]] object.

    real(wp), intent(out) :: r
    !! The radial coordinate, \(r\) (in Schwarzschild coordinates).
    real(wp), intent(out) :: drdt
    !! The time derivative of the radial coordinate, \(\dot{r}\).
    real(wp), intent(out) :: d2rdt2
    !! The second time derivative of the radial coordinate, \(\ddot{r}\).

    r = this%r
    drdt = this%drdt
    d2rdt2 = this%d2rdt2
  end subroutine get_tdc


end module orbit_base
