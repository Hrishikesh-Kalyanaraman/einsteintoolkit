module effective_source
!! Module that defines the abstract interface of an effective source class.

  use kinds
  use grid_function
  use world_tube
  use iso_c_binding

  implicit none

  type, abstract :: eff_source
  !! An abstract effective source interface.
    integer(ip) :: nmodes
    !! The number of modes an effective source is provided for.
    integer(ip) :: nvars
    !! The number of variables an effective source is provided for.
    type(cgf), dimension(:,:), allocatable :: source
    !! A 2d-array of complex grid functions. When allocated the size is
    !! ([[eff_source:nvars]]:[[eff_source:nmodes]]).
  contains
    procedure (eff_source_init), deferred, pass :: init
    !! The constructor.
!    procedure (eff_source_set_window), deferred, pass :: set_window
!    procedure (eff_source_calc_window), deferred, pass :: calc_window    
    procedure (eff_source_set_time_window), deferred, pass :: set_time_window
    !! Routine to set a time window.
    procedure (eff_source_set_particle_pos), deferred, pass :: set_particle_pos
    !! Routine to specify the current state of the particle.
    procedure (eff_source_evaluate_source), deferred, pass :: evaluate_source
    !! Routine to evaluate the effective source.
    procedure (eff_source_get_singular), deferred, pass :: get_singular
    !! Routine to get the singular field at a given radius.
    procedure (eff_source_get_dsingular_dt), deferred, pass :: get_dsingular_dt
    !! Routine to get the time derivative of the singular field at a given
    !! radius.
    procedure (eff_source_get_dsingular_dr), deferred, pass :: get_dsingular_dr
    !! Routine to get the radial derivative of the singular field at a given
    !! radius.
  end type eff_source


  abstract interface
    subroutine eff_source_init ( this, nmodes, nvars, l, m, mass )
    !! The interface of the constructor.
      import :: eff_source, ip, wp, c_int
      class(eff_source), intent(inout) :: this
      !! On return, the constructed object.
      integer(ip), intent(in) :: nmodes
      !! The number of modes.
      integer(ip), intent(in) :: nvars
      !! The number of variables.
      integer(c_int), dimension(nmodes), intent(in) :: l
      !! A 1d array of size [[eff_source:nmodes]] containing the \(\ell\)-values
      !! of the modes.
      integer(c_int), dimension(nmodes), intent(in) :: m
      !! A 1d array of size [[eff_source:nmodes]] containing the m-values
      !! of the modes.
      real(wp), intent(in) :: mass
      !! The mass of the black hole.
    end subroutine eff_source_init

!    subroutine eff_source_set_window ( this, r1, w1, q1, s1, &
!                                       r2, w2, q2, s2 )
!      import :: eff_source, wp
!      class(eff_source(*,*)), intent(in) :: this
!      real(wp), intent(in) :: r1, w1, q1, s1, r2, w2, q2, s2
!    end subroutine eff_source_set_window 
!
!    subroutine eff_source_calc_window ( this, r )
!      import :: eff_source, rgf
!      class(eff_source(*,*)), intent(inout) :: this
!      type(rgf), intent(in) :: r
!    end subroutine eff_source_calc_window
!
    subroutine eff_source_set_time_window ( this, tfac, dtfac_dt, d2tfac_dt2, &
                                            do_smooth_after_lmax )
    !! The interface of the routine that sets the time window for smooth turn
    !! on of the effective source.
      import :: eff_source, ip, wp
      class(eff_source), intent(inout) ::  this
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
    end subroutine eff_source_set_time_window

    subroutine eff_source_set_particle_pos ( this, r, phi, ur, En, Lz, &
                                             ar, aphi, dardt, daphidt, &
                                             d2ardt2, d2aphidt2 )
    !! The interface of the routine that sets the state of the particle.
      import :: eff_source, wp
      class(eff_source), intent(inout) :: this
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
    end subroutine eff_source_set_particle_pos

    subroutine eff_source_evaluate_source ( this, r, wt )
    !! The interface of the routine that evaluates the effective source.
    !!
    !! On return [[eff_source:source]] contains the evaluated effective source.
      import :: eff_source, rgf, wtube
      class(eff_source), intent(inout) :: this
      !! The routine is called on this object.
      type(rgf), intent(in) :: r
      !! A real values grid function that contain the radial coordinate
      !! (in Schwarzschild coordinates).
      type(wtube), intent(in) :: wt
      !! A world-tube object that ensures that the effective source is only
      !! non-zero inside the world-tube.
    end subroutine eff_source_evaluate_source

    subroutine eff_source_get_singular ( this, r, mode, psi )
    !! The interface of the routine that evaluates the singular field for
    !! a given mode and at a given radial coordinate.
      import :: eff_source, wp, ip
      class(eff_source), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: r
      !! The radial coordinate (Schwarzschild coordinates).
      integer(ip), intent(in) :: mode
      !! The mode.
      complex(wp), dimension(:), intent(out) :: psi
      !! A 1d-array of size [[eff_source:nvars]] of complex values that on
      !! return contains the singular field for all variables.
    end subroutine eff_source_get_singular

    subroutine eff_source_get_dsingular_dt ( this, r, sgn, mode, dpsidt )
    !! The interface of the routine that evaluates the time derivative of the
    !! singular field for a given mode and at a given radial coordinate.
      import :: eff_source, wp, ip
      class(eff_source), intent(inout) :: this
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
      !! A 1d-array of size [[eff_source:nvars]] of complex values that on
      !! return contains the time derivative of the singular field for all
      !! variables.
    end subroutine eff_source_get_dsingular_dt

    subroutine eff_source_get_dsingular_dr ( this, r, sgn, mode, dpsidr )
    !! The interface of the routine that evaluates the radial derivative of the
    !! singular field for a given mode and at a given radial coordinate.
      import :: eff_source, wp, ip
      class(eff_source), intent(inout) :: this
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
    end subroutine eff_source_get_dsingular_dr

  end interface

end module effective_source
