module method_of_lines
!! Module that defines an abstract class for the concept of an ODE integrator.
!! As this is just an abstract class, there is no implementation.

  use kinds
  use equations

  implicit none

  type, abstract :: integrator
  !! An abstract class that defines an equation time integrator.
    integer(ip) :: nequations
    !! The number of equations to integrate.
    class(equation_pointer), dimension(:), allocatable :: eqs
    !! A 1d array of equation pointers. Will be of length nequations.
    character(:), allocatable :: iname
    !! An allocatable string that will contain the name of the integrator.
    logical :: is_complete
    !! When this variable is .false. the integrator is in the process of
    !! doing a substep, while when .true. a step has been completed.
  contains
    procedure (integrator_ntemp_interface), deferred, pass :: ntemp
    !! A procedure that allows an integrator to tell the equations how much
    !! temporary storage is needed.
    procedure (integrator_init_interface), deferred, pass :: init
    !! A procedure that initializes an integrator.
    procedure (integrator_step_interface), deferred, pass :: step
    !! A procedure that takes one time step.
    procedure (integrator_complete_step), deferred, pass :: complete_step
    !! A procedure that reports if the integrator have completed a
    !! step or is in the process of taking an intermediate step.
    procedure (integrator_shutdown_interface), deferred, pass :: shutdown
    !! A procedure that shuts down an integrator.
  end type integrator

  abstract interface
    function integrator_ntemp_interface ( this ) result (ntemp)
    !! The return value is the number of temporary storage levels are needed.
      import :: integrator, ip
      class(integrator), intent(in) :: this
      !! The routine is called on this object.
      integer(ip) :: ntemp
      !! The number of temporary storage levels needed.
    end function integrator_ntemp_interface

    subroutine integrator_init_interface ( this, eqs )
    !! Initialize an integrator.
      import :: integrator, equation_pointer
      class(integrator), intent(inout) :: this
      !! The routine is called on this object.
      type(equation_pointer), dimension(:), intent(in) :: eqs
      !! A 1d array of pointers to the equations that will be integrated.
    end subroutine integrator_init_interface

    subroutine integrator_step_interface ( this )
    !! Take a time step.
      import :: integrator
      class(integrator), intent(inout) :: this
      !! The routine is called on this object.
    end subroutine integrator_step_interface

    function integrator_complete_step ( this ) result (is_complete)
    !! Report if the integrator is in the process of taking an intermediate
    !! step or whether the step is completed.
      import :: integrator
      class(integrator), intent(inout) :: this
      !! The routine is called on this object.
      logical :: is_complete
      !! If .false. is returned the integrator is currently doing a substep.
    end function integrator_complete_step

    subroutine integrator_shutdown_interface ( this )
    !! Shut down this integrator.
      import :: integrator
      class(integrator), intent(inout) :: this
      !! The routine is called on this object.
    end subroutine integrator_shutdown_interface
  end interface

end module method_of_lines
