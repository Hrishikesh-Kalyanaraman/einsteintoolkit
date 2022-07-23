module self_force_base
!! Module that defines a global self-force object that can be used to pass
!! information between different equations using the save_globals/load_globals
!! mechanism.

  use kinds

  type :: self_force
  !! A global self-force type.
    real(wp), dimension(4), private :: f = 0.0_wp
    !! A 4-vector containing the components of the self-force, \(f_{\mu}\).
    real(wp), dimension(4), private :: a = 0.0_wp
    !! A 4-vector containing the components of the 4-acceleration, \(a^{\mu}\).
    real(wp), dimension(4), private :: da_dt = 0.0_wp
    !! A 4-vector containing the components of the coordinate time derivatives
    !! of the 4-acceleration, \(\dot{a}^{\mu}\).
    real(wp), dimension(4), private :: d2a_dt2 = 0.0_wp
    !! A 4-vector containing the components of the second coordinate time
    !! derivatives of the 4-acceleration, \(\ddot{a}^{\mu}\).
    real(wp), private :: udota = 0.0_wp
    !! A scalar containing the inner product of the 4-velocity and the
    !! self-force \(u^{\mu}f_{\mu}\).
    integer(ip), private :: ioo_id = -1
    !! The file unit number used for output.
  contains
    procedure :: set_force
    !! Routine to set the self-force values.
    procedure :: get_force
    !! Routine to get the self-force values.
    procedure :: set_accel
    !! Routine to set the values of the 4-acceleration.
    procedure :: set_daccel_dt
    !! Routine to set the values of the time derivative of the 4-acceleration.
    procedure :: set_d2accel_dt2
    !! Routine to set the values of the second time derivative of the
    !! 4-acceleration.
    procedure :: get_accel
    !! Routine to get the values of the 4-acceleration.
    procedure :: get_daccel_dt
    !! Routine to get the values of the time derivative of the 4-acceleration.
    procedure :: get_d2accel_dt2
    !! Routine to get the values of the second time derivative of the
    !! 4-acceleration.
    procedure :: output
    !! Routine to perform output.
  end type self_force

  type(self_force) :: sf
  !! The global [[self_force]] object that is available by use association.

contains

  subroutine set_force ( this, ft, fr, ftheta, fphi )
  !! Routine that sets the self-force variable.

    implicit none

    class(self_force), intent(inout) :: this
    !! The routine is called on this [[self_force]] object.
    real(wp), intent(in) :: ft
    !! The \(t\)-component of the self-force, \(f_t\).
    real(wp), intent(in) :: fr
    !! The \(r\)-component of the self-force, \(f_r\).
    real(wp), intent(in) :: ftheta
    !! The \(\theta\)-component of the self-force, \(f_{\theta}\).
    real(wp), intent(in) :: fphi
    !! The \(\phi\)-component of the self-force, \(f_{\phi}\).

    this%f = (/ ft, fr, ftheta, fphi /)
  end subroutine set_force


  subroutine get_force ( this, force )
  !! Routine that gets the self-force variable.
  !!
  !! This also take into account run time parameters that specifies when
  !! the back-reaction should be turned on and will return zero if it is
  !! not yet time to do so.

    use time_info, only : get_current_time
    use orbit_base, only : orbit_info
    use parameters, only : evolve_orbit, turn_on_force_smoothly, use_chi, &
                           evolve_after, force_sigma, torder
    use numerics, only : time_window

    implicit none

    class(self_force), intent(inout) :: this
    !! The routine is called on this [[self_force]] object.
    real(wp), dimension(4), intent(out) :: force
    !! On return this 1d array contains the 4 components of the force,
    !! \((f_t,f_r,f_{\theta},f_{\phi})\).
    real(wp) :: time, chi, f_fac, df_fac_dt, d2f_fac_dt2

    f_fac = 0.0_wp

    if ( evolve_orbit ) then
      if ( use_chi .and. turn_on_force_smoothly ) then
        call orbit_info%get_chi ( chi )
        if ( chi >= evolve_after ) then
          call time_window ( chi-evolve_after, force_sigma, torder, &
                             f_fac, df_fac_dt, d2f_fac_dt2 )
        end if
      else if ( use_chi ) then
        f_fac = 1.0_wp
      end if

      if ( ( .not. use_chi ) .and. turn_on_force_smoothly ) then
        time = get_current_time ( )
        if ( time >= evolve_after ) then
          call time_window ( time-evolve_after, force_sigma, torder, &
                             f_fac, df_fac_dt, d2f_fac_dt2 )
        end if
      else if ( .not. use_chi ) then
        f_fac = 1.0_wp
      end if 
      force = f_fac*this%f
    else
      force = 0.0_wp
    end if


  end subroutine get_force


  subroutine set_accel ( this, at, ar, atheta, aphi, udota )
  !! Routine that sets the 4-acceleration and the inner product of the
  !! 4-velovity and the self-force.

    implicit none

    class(self_force), intent(inout) :: this
    !! The routine is called on this [[self_force]] object.
    real(wp), intent(in) :: at
    !! The \(t\)-component of the 4-acceleration, \(a^t\).
    real(wp), intent(in) :: ar
    !! The \(r\)-component of the 4-acceleration, \(a^r\).
    real(wp), intent(in) :: atheta
    !! The \(\theta\)-component of the 4-acceleration, \(a^{\theta}\).
    real(wp), intent(in) :: aphi
    !! The \(\phi\)-component of the 4-acceleration, \(a^{\phi}\).
    real(wp), intent(in) :: udota
    !! The inner product of the 4-velocity and the self-force,
    !! \(u^{\mu}f_{\mu}\).

    this%a = (/ at, ar, atheta, aphi /)
    this%udota = udota
  end subroutine set_accel


  subroutine set_daccel_dt ( this, dat_dt, dar_dt, datheta_dt, daphi_dt )
  !! Routine that sets the coordinate time derivative of the 4-acceleration.

    implicit none

    class(self_force), intent(inout) :: this
    !! The routine is called on this [[self_force]] object.
    real(wp), intent(in) :: dat_dt
    !! The time derivative of the \(t\)-component of the 4-acceleration,
    !! \(\dot{a}^t\).
    real(wp), intent(in) :: dar_dt
    !! The time derivative of the \(r\)-component of the 4-acceleration,
    !! \(\dot{a}^r\).
    real(wp), intent(in) :: datheta_dt
    !! The time derivative of the \(\theta\)-component of the 4-acceleration,
    !! \(\dot{a}^{\theta}\).
    real(wp), intent(in) :: daphi_dt
    !! The time derivative of the \(\phi\)-component of the 4-acceleration,
    !! \(\dot{a}^{\phi}\).

    this%da_dt = (/ dat_dt, dar_dt, datheta_dt, daphi_dt /)
  end subroutine set_daccel_dt


  subroutine set_d2accel_dt2 ( this, d2at_dt2, d2ar_dt2, &
                                     d2atheta_dt2, d2aphi_dt2 )
  !! Routine that sets the second coordinate time derivative of the
  !! 4-acceleration.

    implicit none

    class(self_force), intent(inout) :: this
    !! The routine is called on this [[self_force]] object.
    real(wp), intent(in) :: d2at_dt2
    !! The second time derivative of the \(t\)-component of the 4-acceleration,
    !! \(\ddot{a}^t\).
    real(wp), intent(in) :: d2ar_dt2
    !! The second time derivative of the \(r\)-component of the 4-acceleration,
    !! \(\ddot{a}^r\).
    real(wp), intent(in) :: d2atheta_dt2
    !! The second time derivative of the \(\theta\)-component of the
    !! 4-acceleration, \(\ddot{a}^{\theta}\).
    real(wp), intent(in) :: d2aphi_dt2
    !! The second time derivative of the \(\phi\)-component of the
    !! 4-acceleration, \(\ddot{a}^{\phi}\).

    this%d2a_dt2 = (/ d2at_dt2, d2ar_dt2, d2atheta_dt2, d2aphi_dt2 /)
  end subroutine set_d2accel_dt2


  subroutine get_accel ( this, accel )
  !! Routine that gets the 4-acceleration.

    implicit none

    class(self_force), intent(inout) :: this
    !! The routine is called on this [[self_force]] object.
    real(wp), dimension(4), intent(out) :: accel
    !! On return this 1d array contains the 4 components of the 4-acceleration,
    !! \((a^t,a^r,a^{\theta},a^{\phi})\).

    accel = this%a

  end subroutine get_accel


  subroutine get_daccel_dt ( this, daccel_dt )
  !! Routine that gets the coordinate time derivative of the 4-acceleration.

    implicit none

    class(self_force), intent(inout) :: this
    !! The routine is called on this [[self_force]] object.
    real(wp), dimension(4), intent(out) :: daccel_dt
    !! On return this 1d array contains the 4 components of the time derivative
    !! of the 4-acceleration,
    !! \((\dot{a}^t,\dot{a}^r,\dot{a}^{\theta},\dot{a}^{\phi})\).

    daccel_dt = this%da_dt

  end subroutine get_daccel_dt


  subroutine get_d2accel_dt2 ( this, d2accel_dt2 )
  !! Routine that gets the second coordinate time derivative of the
  !! 4-acceleration.

    implicit none

    class(self_force), intent(inout) :: this
    !! The routine is called on this [[self_force]] object.
    real(wp), dimension(4), intent(out) :: d2accel_dt2
    !! On return this 1d array contains the 4 components of the second time
    !! derivative of the 4-acceleration,
    !! \((\ddot{a}^t,\ddot{a}^r,\ddot{a}^{\theta},\ddot{a}^{\phi})\).

    d2accel_dt2 = this%d2a_dt2

  end subroutine get_d2accel_dt2


  subroutine output ( this )
  !! Routine that performs output of the self-force.

    use output_base
    use time_info, only : get_current_time

    implicit none

    class(self_force), intent(inout) :: this
    !! The routine is called on this [[self_force]] object.
    integer(ip) :: ioo_id, tmp_id
    character(len=14) :: filename = 'self_force.asc'

    ioo_id = this%ioo_id
    if (ioo_id<0) then
      ioo_id = next_available_io_id ()
      this%ioo_id = ioo_id
      print*,'Opening ', filename, ' with id ', ioo_id
      open(ioo_id, file=filename, status='replace', action='write')
    end if

    write(ioo_id, '(*(es23.15e3,1x))') get_current_time ( ), this%f, this%a, &
                                       this%da_dt, this%d2a_dt2
  end subroutine output

end module self_force_base
