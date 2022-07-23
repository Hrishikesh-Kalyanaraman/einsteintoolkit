module equations
!! Module that defines the abstract interface of an equation class.

  use kinds

  implicit none

  type, abstract :: equation
  !! An abstract equation interface.
    integer(ip) :: ntmp
    !! The number of temporary storage variables to allocate.
    character(:), allocatable :: ename
    !! The name of the system of equations
  contains
    procedure (eq_init_interface), deferred, pass :: init
    !! The initialization routine. Implementation is deferred to the
    !! derived class that actually implements an equation system.
    procedure (eq_rhs_interface), deferred, pass :: rhs
    !! The right hand side (RHS) routine. Implementation is deferred to the
    !! derived class that actually implements an equation system.
    procedure (eq_set_to_zero_interface), deferred, pass :: set_to_zero
    !! Routine to set the main, rhs or temporary variables to zero.
    !! Implementation is deferred to the derived class that actually
    !! implements an equation system.
    procedure (eq_update_vars_interface), deferred, pass :: update_vars
    !! Routine to update the main or temporary variables. Implementation is
    !! deferred to the derived class that actually implements an equation
    !! system.
    procedure (eq_save_globals_1), deferred, pass :: save_globals_1
    !! First of 2 routines to set global variables to interact with other
    !! equation systems. Implementation is deferred to the derived class that
    !! actually implements an equation system.
    procedure (eq_save_globals_2), deferred, pass :: save_globals_2
    !! Second of 2 routines to set global variables to interact with other
    !! equation systems. Implementation is deferred to the derived class that
    !! actually implements an equation system.
    procedure (eq_load_globals), deferred, pass :: load_globals
    !! Routine to get global variables from other equation systems.
    !! Implementation is deferred to the derived class that actually
    !! implements an equation system.
    procedure (eq_output), deferred, pass :: output
    !! Routine to perform output of the equation systems variables.
    !! Implementation is deferred to the derived class that actually
    !! implements an equation system.
    procedure (eq_print_data), deferred, pass :: print_data
    !! Routine to perform debug output of the equation systems variables.
    !! Implementation is deferred to the derived class that actually
    !! implements an equation system.
  end type equation

  abstract interface
    subroutine eq_init_interface ( this )
    !! The initialization routine. There are no additional input arguments.
    !! The routine gets it's input from the run-time parameters.
      import :: equation, ip
      class(equation), target, intent(inout) :: this
      !! The equation that is being initialized.
    end subroutine eq_init_interface

    subroutine eq_rhs_interface ( this )
    !! The RHS routine that sets the RHS variables from the current state of
    !! the system.
      import :: equation
      class(equation), intent(inout) :: this
      !! The equation for which the RHS is calculated.
    end subroutine eq_rhs_interface

    subroutine eq_set_to_zero_interface ( this, dest )
    !! Set either the main, RHS or temporary variables to zero.
      import :: equation, ip
      class(equation), intent(inout) :: this
      !! The routine is called on this equation.
      integer(ip), intent(in) :: dest
      !! Can be either -1 (RHS), 0 (main) or between 1 and
      !! [[equation:ntmp(variable)]] (temporary).
    end subroutine eq_set_to_zero_interface

    subroutine eq_update_vars_interface ( this, source, dest, source2, &
                                                 scalar, scalar2 )
    !! The routine that updates the variables defined in an equation. This
    !! is used to make the time integrator agnostic to how storage for the
    !! variables in the system of equations are set up.
    !!
    !! After completion:
    !! this%var(dest) = scalar\*this%var(source)+scalar2\*this%var(source2)
      import :: equation, ip, wp
      class(equation), target, intent(inout) :: this
      !! The routine is called on this equation.
      integer(ip), intent(in) :: source
      !! The first source for the update. Can be either -1 (RHS), 0 (main) or
      !! between 1 and [[equation:ntmp(variable)]] (temporary).
      integer(ip), intent(in) :: dest
      !! The destination for the update. Can be either -1 (RHS), 0 (main) or
      !! between 1 and [[equation:ntmp(variable)]] (temporary).
      integer(ip), optional, intent(in) :: source2
      !! The second source for the update. Can be either -1 (RHS), 0 (main) or
      !! between 1 and [[equation:ntmp(variable)]] (temporary).
      real(wp), optional, intent(in) :: scalar
      !! The scalar multiplying the first source.
      real(wp), optional, intent(in) :: scalar2
      !! The scalar multiplying the second source.
    end subroutine eq_update_vars_interface

    subroutine eq_save_globals_1 ( this )
    !! The routine where the equation can set global variables (for
    !! communication with other equations) for the first time.
      import :: equation
      class(equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine eq_save_globals_1

    subroutine eq_save_globals_2 ( this )
    !! The routine where the equation can set global variables (for
    !! communication with other equations) for the first time.
      import :: equation
      class(equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine eq_save_globals_2

    subroutine eq_load_globals ( this )
    !! The routine where the equation can get global variables (for
    !! communication with other equations).
      import :: equation
      class(equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine eq_load_globals

    subroutine eq_output ( this )
    !! The routine where output of the equation's variables is done.
      import :: equation
      class(equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine eq_output

    subroutine eq_print_data ( this )
    !! The routine where debug output (to stdout) of the equation's variables
    !! is done.
      import :: equation
      class(equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine eq_print_data
  end interface

  type :: equation_pointer
  !! A type with pointer that can point to the equation class and any derived
  !! class.
    class(equation), pointer :: p
    !! A pointer to the equation class.
  end type equation_pointer

end module equations
