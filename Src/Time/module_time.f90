module time_info
!! Module with variables and routines to keep track of time.

!! There are both working and quad precision copies of the time variable that
!! gets updated synchronously when the routines provided here are used.

  use kinds

  implicit none

  real(wp), private :: time
  !! Working precision copy of the time, \(t\).
  real(wp), private :: time_save
  !! Working precision backup copy of the time.
  real(wp), private :: dtime
  !! \(\Delta t\).
  real(qp), private :: qtime
  !! Quad precision copy of the time, \(t\).
  real(qp), private :: qtime_save
  !! Quad precision backup copy of the time.
  logical :: short_timesteps_active = .false.
  !! Variable to keep track of whether small \(\Delta t\) is used for quick
  !! but smooth turn on of the effective source.

contains

  subroutine init_time ( t0 )
  !! Routine to initialize the time variables.

    implicit none

    real(wp), intent(in) :: t0
    !! The initial time, \(t_0\).

    time = t0
    qtime = t0
 
  end subroutine init_time


  subroutine save_time ( )
  !! Routine to make a backup copy of the time.

    implicit none

    time_save = time
    qtime_save = qtime
  end subroutine save_time  


  subroutine restore_and_increment_time ( dt )
  !! Routine to restore and increment a backup copy of the time.
  !!
  !! Needed for some Runge-Kutta integrators.

    implicit none

    real(wp), intent(in) :: dt
    !! The \(\Delta t\) to use for the increment.

    time = time_save + dt
    qtime = qtime_save + dt

  end subroutine restore_and_increment_time


  subroutine set_dtime ( dt )
  !! Routine that sets \(\Delta t\).

    implicit none
 
    real(wp), intent(in) :: dt
    !! The value to use for \(\Delta t\).

    dtime = dt

  end subroutine set_dtime


  subroutine increment_time ( dt )
  !! Routine to increment \(t\) by \(\Delta t\).

    implicit none

    real(wp), intent(in) :: dt
    !! The value to use for \(\Delta t\).

    time = time + dt
    qtime = qtime + dt

  end subroutine increment_time


  function get_current_time ( ) result(t)
  !! Function to get the current working precision time, \(t\).

    implicit none

    real(wp) :: t
    !! Returns \(t\).

    t = time

  end function get_current_time


  function get_current_dtime ( ) result(dt)
  !! Function to get the current \(\Delta t\)

    implicit none

    real(wp) :: dt
    !! Returns \(\Delta t\).

    dt = dtime

  end function get_current_dtime


  function get_current_sub_dtime ( ) result(dt)
  !! Function to get the current \(\Delta t\) of an integrator substep.

    implicit none

    real(wp) :: dt
    !! Returns \(\Delta t\).

    dt = real(qtime-qtime_save,wp)

  end function get_current_sub_dtime


  function get_current_qtime ( ) result(t)
  !! Function to get the current quad precision time, \(t\).

    implicit none

    real(qp) :: t
    !! Returns \(t\).

    t = qtime

  end function get_current_qtime

end module time_info
