module rk5_integrator
!! Module that provides a 5th order continuous explicit Runge-Kutta ODE
!! integrator (see Verner & Zennaro, 1995, Mathematics of Computation, 64, 211,
!! 1123-1146).
!!
!! The implementation is found in [[submodule_rk5_implementation.f90]].

  use kinds
  use method_of_lines

  implicit none

  type, extends(integrator) :: rk5
  !! A 5th order continuous explicit Runge-Kutta ODE integrator.
    integer(ip) :: ntmp = 7
    !! 7 levels of temporary storage are required.
  contains
    procedure :: ntemp => rk5_ntemp
    !! Routine to provide information about temporary storage levels is
    !! provided by [[rk5_ntemp]].
    procedure :: init => rk5_init
    !! Initialization routine is provided by [[rk5_init]].
    procedure :: step => rk5_step
    !! Stepping routine is provided by [[rk5_step]].
    procedure :: complete_step => rk5_complete_step
    !! Complete step routine is provided by [[rk5_complete_step]].
    procedure :: shutdown => rk5_shutdown
    !! Shut down routine is provided by [[rk5_shutdown]].
  end type rk5

  interface
    module function rk5_ntemp ( this ) result (ntemp)
    !! Routine that reports how many temporary storage levels are needed.
      class(rk5), intent(in) :: this
      !! The routine is called on this object.
      integer(ip) :: ntemp
     !! The return value is the number of required temorary storage levels.
    end function rk5_ntemp

    module subroutine rk5_init ( this, eqs )
    !! Routine that initializes the integrator.
      class(rk5), intent(inout) :: this
      !! The routine is called on this object.
      type(equation_pointer), dimension(:), intent(in) :: eqs
      !! A 1d-array of pointers to equations that will be integrated.
    end subroutine rk5_init

    module subroutine rk5_step ( this )
    !! Routine that takes a time step.
      class(rk5), intent(inout) :: this
      !! The routine is called on this object.
    end subroutine rk5_step

    module function rk5_complete_step ( this ) result (is_complete)
    !! Report if the integrator is in the process of taking an intermediate
    !! step or whether the step is completed.
      class(rk5), intent(inout) :: this
      !! The routine is called on this object.
      logical :: is_complete
      !! If .false. is returned the integrator is currently doing a substep.
    end function rk5_complete_step

    module subroutine rk5_shutdown ( this )
    !! Routine that shuts downs the integrator.
      class(rk5), intent(inout) :: this
      !! The routine is called on this object.
    end subroutine rk5_shutdown
  end interface

end module rk5_integrator
