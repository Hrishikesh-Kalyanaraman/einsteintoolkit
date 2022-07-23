module rk4_integrator
!! Module that provides a 4th order low storage Runge-Kutta ODE integrator.
!!
!! The implementation is found in [[submodule_rk4_implementation.f90]].

  use kinds
  use method_of_lines

  implicit none

  type, extends(integrator) :: rk4
  !! A 4th order low storage Runge-Kutta ODE integrator.
    integer(ip) :: ntmp = 1
    !! 1 level of temporary storage are required.
  contains
    procedure :: ntemp => rk4_ntemp
    !! Routine to provide information about temporary storage levels is
    !! provided by [[rk4_ntemp]].
    procedure :: init => rk4_init
    !! Initialization routine is provided by [[rk4_init]].
    procedure :: step => rk4_step
    !! Stepping routine is provided by [[rk4_step]].
    procedure :: complete_step => rk4_complete_step
    !! Complete step routine is provided by [[rk4_complete_step]].
    procedure :: shutdown => rk4_shutdown
    !! Shut down routine is provided by [[rk4_shutdown]].
  end type rk4

  interface
    module function rk4_ntemp ( this ) result (ntemp)
    !! Routine that reports how many temporary storage levels are needed.
      class(rk4), intent(in) :: this
      !! The routine is called on this object.
      integer(ip) :: ntemp
      !! The return value is the number of required temorary storage levels.
    end function rk4_ntemp

    module subroutine rk4_init ( this, eqs )
    !! Routine that initializes the integrator.
      class(rk4), intent(inout) :: this
      !! The routine is called on this object.
      type(equation_pointer), dimension(:), intent(in) :: eqs
      !! A 1d-array of pointers to equations that will be integrated.
    end subroutine rk4_init

    module subroutine rk4_step ( this )
    !! Routine that takes a time step.
      class(rk4), intent(inout) :: this
      !! The routine is called on this object.
    end subroutine rk4_step

    module function rk4_complete_step ( this ) result (is_complete)
    !! Report if the integrator is in the process of taking an intermediate
    !! step or whether the step is completed.
      class(rk4), intent(inout) :: this
      !! The routine is called on this object.
      logical :: is_complete
      !! If .false. is returned the integrator is currently doing a substep.
    end function rk4_complete_step

    module subroutine rk4_shutdown ( this )
    !! Routine that shuts downs the integrator.
      class(rk4), intent(inout) :: this
      !! The routine is called on this object.
    end subroutine rk4_shutdown
  end interface

end module rk4_integrator
