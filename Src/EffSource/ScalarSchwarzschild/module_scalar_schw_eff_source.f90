module scalar_schw_eff
!! Module that defines an effective source class for a generic effective source
!! for a scalar charge in orbit around a Schwarzschild black hole.
!!
!! The implementation in
!! [[submodule_scalar_schw_eff_source_implementation.f90]] is an interface to
!! a C++ class provided by Barry Wardell.
  use kinds
  use effective_source
  use iso_c_binding

  implicit none

  type, extends(eff_source) :: scal_schw_eff
  !! A class that interfaces with a C++ class that provides an effective
  !! source for a scalar charge on a generic accelerated orbit around a
  !! Schwarzschild black hole.
    real(wp), dimension(:,:), allocatable :: sre
    !! A 2d-array of reals to hold the real part of the effective source
    !! for all nodes in a DG element and for all modes. When allocated
    !! the size is ([[element_data:order]]+1:[[eff_source:nmodes]]).
    real(wp), dimension(:,:), allocatable :: sim
    !! A 2d-array of reals to hold the complex part of the effective source
    !! for all nodes in a DG element and for all modes. When allocated
    !! the size is ([[element_data:order]]+1:[[eff_source:nmodes]]).
  contains
    procedure :: init => scal_schw_eff_init
    !! The [[eff_source:init]] routine is provided by [[scal_schw_eff_init]].
    procedure :: set_time_window => scal_schw_eff_set_time_window
    !! The [[eff_source:set_particle_pos]] routine is provided by
    !! [[scal_schw_eff_set_particle_pos]].
    procedure :: set_particle_pos => scal_schw_eff_set_particle_pos
    !! The [[eff_source:set_particle_pos]] routine is provided by
    !! [[scal_schw_eff_set_particle_pos]].
    procedure :: evaluate_source => scal_schw_eff_evaluate_source
    !! The [[eff_source:evaluate_source]] routine is provided by
    !! [[scal_schw_eff_evaluate_source]].
    procedure :: get_singular => scal_schw_eff_get_singular
    !! The [[eff_source:get_singular]] routine is provided by
    !! [[scal_schw_eff_get_singular]].
    procedure :: get_dsingular_dt => scal_schw_eff_get_dsingular_dt
    !! The [[eff_source:get_dsingular_dt]] routine is provided by
    !! [[scal_schw_eff_get_dsingular_dt]].
    procedure :: get_dsingular_dr => scal_schw_eff_get_dsingular_dr
    !! The [[eff_source:get_dsingular_dr]] routine is provided by
    !! [[scal_schw_eff_get_dsingular_dr]].
  end type scal_schw_eff

  interface
    module subroutine scal_schw_eff_init ( this, nmodes, nvars, l, m, mass )
    !! Interface of the constructor compatible with [[eff_source:init]].
      class(scal_schw_eff), intent(inout) :: this
      !! On return, the constructed object.
      integer(ip), intent(in) :: nmodes
      !! The number of modes.
      integer(ip), intent(in) :: nvars
      !! The number of variables. 
      integer(c_int), dimension(nmodes), intent(in) :: l
      !! A 1d array of size [[eff_source:nmodes]] containing the \(\ell\)-values
      !! of the modes.
      integer(c_int), dimension(nmodes), intent(in) :: m
      !! A 1d array of size [[eff_source:nmodes]]  containing the m-values
      !! of the modes.
      real(wp), intent(in) :: mass
      !! The mass of the black hole.
    end subroutine scal_schw_eff_init 

    module subroutine scal_schw_eff_set_time_window ( this, tfac, dtfac_dt, &
                                                      d2tfac_dt2, &
                                                      do_smooth_after_lmax )
    !! Interface to the set_time_window routine compatible with
    !! [[eff_source:set_time_window]].
      class(scal_schw_eff), intent(inout) ::  this
      !! The routine is called on this object.
      real(wp), intent(in) :: tfac
      !! The current value for the time window.
      real(wp), intent(in) :: dtfac_dt
      !! The current value for the time derivative of the time window.
      real(wp), intent(in) :: d2tfac_dt2
      !! The current value for the second time derivative of the time window.
      integer(ip), intent(in) :: do_smooth_after_lmax
      !! If \(l<\)do_smooth_after_lmax use tfac=1, dtfac_dt=0, d2tfac_dt2=0.
      !! This allows for smooth turn on of the effective source for the modes
      !! with \(l>\)do_smooth_after_lmax, while the rest gets turned on
      !! instantaneously (e.g. when external initial data is available).
    end subroutine scal_schw_eff_set_time_window

    module subroutine scal_schw_eff_set_particle_pos ( this, r, phi, ur, En, &
                                                       Lz, ar, aphi, dardt, &
                                                       daphidt, d2ardt2, &
                                                       d2aphidt2 )
    !! Interface of the routine that sets the state of the particle.
      class(scal_schw_eff), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: r
      !! The radial coordinate of the particle (in Schwarzschild coordinates),
      !! \(r\).
      real(wp), intent(in) :: phi
      !! The azimuthal angle of the particle, \(\phi\).
      real(wp), intent(in) :: ur
      !! The radial component of the 4-velocity, \(u^r\).
      real(wp), intent(in) :: En
      !! The energy per unit mass of the particle, \(E\).
      real(wp), intent(in) :: Lz
      !! The angular momentum per unit mass of the particle \(L_z\).
      real(wp), intent(in) :: ar
      !! The radial component of the 4-acceleration, \(a^r\).
      real(wp), intent(in) :: aphi
      !! The \(\phi\)-component of the 4-acceleration, \(a^{\phi}\).
      real(wp), intent(in) :: dardt
      !! The derivative of \(a^r\) with respect to Schwarzschild coordinate
      !! time.
      real(wp), intent(in) :: daphidt
      !! The derivative of \(a^{\phi}\) with respect to Schwarzschild
      !! coordinate time.
      real(wp), intent(in) :: d2ardt2
      !! The second derivative of \(a^r\) with respect to Schwarzschild
      !! coordinate time.
      real(wp), intent(in) :: d2aphidt2
      !! The second derivative of \(a^{\phi}\) with respect to Schwarzschild
      !! coordinate time.
    end subroutine scal_schw_eff_set_particle_pos

    module subroutine scal_schw_eff_evaluate_source ( this, r, wt )
    !! Interface of the routine that evaluates the effective source.
    !!
    !! On return [[scal_schw_eff:source]] contains the evaluated effective
    !! source.
      class(scal_schw_eff), intent(inout) :: this
      !! The routine is called on this object.
      type(rgf), intent(in) :: r
      !! A real values grid function that contain the radial coordinate
      !! (in Schwarzschild coordinates).
      type(wtube), intent(in) :: wt
      !! A world-tube object that ensures that the effective source is only
      !! non-zero inside the world-tube.
    end subroutine scal_schw_eff_evaluate_source

    module subroutine scal_schw_eff_get_singular ( this, r, mode, psi )
    !! Interface of the routine that evaluates the singular field for a given
    !! mode and at a given radial coordinate.
      class(scal_schw_eff), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: r
      !! The radial coordinate (Schwarzschild coordinates).
      integer(ip), intent(in) :: mode
      !! The mode.
      complex(wp), dimension(:), intent(out) :: psi
      !! A 1d-array of size [[scal_schw_eff:nvars]] of complex values that on
      !! return contains the singular field for all variables.
    end subroutine scal_schw_eff_get_singular

    module subroutine scal_schw_eff_get_dsingular_dt ( this, r, sgn, &
                                                       mode, dpsidt )
    !! Interface of the routine that evaluates the time derivative of the
    !! singular field for a given mode and at a given radial coordinate.
      class(scal_schw_eff), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: r
      !! The radial coordinate (Schwarzschild coordinates).
      integer(ip), intent(in) :: sgn
      !! The value should be 1 when the derivative to the right of the
      !! particle is needed and -1 when the derivative to the left of the
      !! particle is needed.
      integer(ip), intent(in) :: mode
      !! The mode.
      complex(wp), dimension(:), intent(out) :: dpsidt
      !! A 1d-array of size [[scal_schw_eff:nvars]] of complex values that on
      !! return contains the time derivative of the singular field for all
      !! variables.
    end subroutine scal_schw_eff_get_dsingular_dt

    module subroutine scal_schw_eff_get_dsingular_dr ( this, r, sgn, &
                                                       mode, dpsidr )
    !! Interface of the routine that evaluates the radial derivative of the
    !! singular field for a given mode and at a given radial coordinate.
      class(scal_schw_eff), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: r
      !! The radial coordinate (Schwarzschild coordinates).
      integer(ip), intent(in) :: sgn
      !! The value should be 1 when the derivative to the right of the
      !! particle is needed and -1 when the derivative to the left of the
      !! particle is needed.
      integer(ip), intent(in) :: mode
      !! The mode.
      complex(wp), dimension(:), intent(out) :: dpsidr
      !! A 1d-array of size [[eff_source:nvars]] of complex values that on
      !! return contains the radial derivative of the singular field for all
      !! variables.
    end subroutine scal_schw_eff_get_dsingular_dr

  end interface

end module scalar_schw_eff
