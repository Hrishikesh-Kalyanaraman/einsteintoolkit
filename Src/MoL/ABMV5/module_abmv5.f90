module abmv5_integrator
!! Module that provides a 5th order Adams-Bashford-Moulton multi-value ODE
!! integrator.
!!
!! The implementation is found in [[submodule_abmv5_implementation.f90]].

  use kinds
  use method_of_lines

  implicit none

  type, extends(integrator) :: abmv5
  !! An Adams-Bashford-Moulton multi-value ODE integrator class.
    integer(ip) :: ntmp = 10
    !! 10 levels of temporary storage are required.
    integer(ip) :: order = 5
    !! The order of the method is 5.
    real(wp) :: last_dt
    !! Variable that keeps track of the last timestep.
  contains
    procedure :: ntemp => abmv5_ntemp
    !! Routine to provide information about temporary storage levels is
    !! provided by [[abmv5_ntemp]].
    procedure :: init => abmv5_init
    !! Initialization routine is provided by [[abmv5_init]].
    procedure :: step => abmv5_step
    !! Stepping routine is provided by [[abmv5_step]].
    procedure :: complete_step => abmv5_complete_step
    !! Complete step routine is provided by [[abmv5_complete_step]].
    procedure :: shutdown => abmv5_shutdown
    !! Shut down routine is provided by [[abmv5_shutdown]].
  end type abmv5

  interface
    module function abmv5_ntemp ( this ) result (ntemp)
    !! Routine that reports how many temporary storage levels are needed.
      class(abmv5), intent(in) :: this
      !! The routine is called on this object.
      integer(ip) :: ntemp
      !! The return value is the number of required temorary storage levels.
    end function abmv5_ntemp

    module subroutine abmv5_init ( this, eqs )
    !! Routine that initializes the integrator.
      class(abmv5), intent(inout) :: this
      !! The routine is called on this object.
      type(equation_pointer), dimension(:), intent(in) :: eqs
      !! A 1d-array of pointers to equations that will be integrated.
    end subroutine abmv5_init

    module subroutine abmv5_step ( this )
    !! Routine that takes a time step.
      class(abmv5), intent(inout) :: this
      !! The routine is called on this object.
    end subroutine abmv5_step

    module function abmv5_complete_step ( this ) result (is_complete)
    !! Report if the integrator is in the process of taking an intermediate
    !! step or whether the step is completed.
      class(abmv5), intent(inout) :: this
      !! The routine is called on this object.
      logical :: is_complete
      !! If .false. is returned the integrator is currently doing a substep.
    end function abmv5_complete_step

    module subroutine abmv5_shutdown ( this )
    !! Routine that shuts downs the integrator.
      class(abmv5), intent(inout) :: this
      !! The routine is called on this object.
    end subroutine abmv5_shutdown
  end interface

end module abmv5_integrator
