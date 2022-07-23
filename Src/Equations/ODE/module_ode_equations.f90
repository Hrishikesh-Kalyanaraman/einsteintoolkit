module ode_equations
!! Module that defines an ODE equation. It is derived from the [[equation]]
!! class and provides the routines that are common for all ODE equation
!! systems while deferring implementation of the routines that are
!! specific to a given ODE equation system.
!!
!! The implementation is found in
!! [[submodule_ode_equations_implementation.f90]].

  use kinds
  use equations

  implicit none

  type, abstract, extends(equation) :: ode_equation
  !! A class derived from [[equation]] specific for ODE equations.
    integer(ip) :: nvars
    !! The number of variables in the ODE system.
    real(wp), dimension(:), allocatable :: var_data
    !! A 1d array of reals that contains the data variables.
    real(wp), dimension(:), allocatable :: rhs_data
    !! A 1d array of reals that contains the RHS variables.
    real(wp), dimension(:,:), allocatable :: tmp_data
    !! A 2d array of reals that contains the temporary variables needed by the
    !! time integrator.
    integer(ip) :: io_id
    !! The file unit id used for output for this system of ODE's.
  contains
    procedure (ode_eq_init_interface), deferred, pass :: init
    !! The [[equation:init]] routine is deferred.
    procedure (ode_eq_rhs_interface), deferred, pass :: rhs
    !! The [[equation:rhs]] routine is deferred.
    procedure :: set_to_zero => ode_set_to_zero
    !! The [[equation:set_to_zero]] routine is provided by
    !! [[ode_set_to_zero]].
    procedure :: update_vars => ode_update_vars
    !! The [[equation:update_vars]] routine is provided by
    !! [[ode_update_vars]].
    procedure (ode_eq_save_globals_1), deferred, pass :: save_globals_1
    !! The [[equation:save_globals_1]] routine is deferred.
    procedure (ode_eq_save_globals_2), deferred, pass :: save_globals_2
    !! The [[equation:save_globals_2]] routine is deferred.
    procedure (ode_eq_load_globals), deferred, pass :: load_globals
    !! The [[equation:load_globals]] routine is deferred.
    procedure (ode_eq_output), deferred, pass :: output
    !! The [[equation:output]] routine is deferred.
    procedure :: print_data => ode_print_data
    !! The [[equation:print_data]] routine is provided by [[ode_print_data]].
  end type ode_equation

  abstract interface
    subroutine ode_eq_init_interface ( this )
    !! The interface for the ODE version of [[equation:init]]. This interface
    !! is consistent with [[eq_init_interface]].
      import :: ode_equation, ip
      class(ode_equation), target, intent(inout) :: this
      !! The equation that is being initialized.
    end subroutine ode_eq_init_interface

    subroutine ode_eq_rhs_interface ( this )
    !! The interface for the ODE version of [[equation:rhs]]. This interface
    !! is consistent with [[eq_rhs_interface]].
      import :: ode_equation
      class(ode_equation), intent(inout) :: this
      !! The equation for which the RHS is calculated.
    end subroutine ode_eq_rhs_interface

    subroutine ode_eq_save_globals_1 ( this )
    !! The interface for the ODE version of [[equation:save_globals_1]].
    !! This interface is consistent with [[eq_save_globals_1]].
      import :: ode_equation
      class(ode_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine ode_eq_save_globals_1

    subroutine ode_eq_save_globals_2 ( this )
    !! The interface for the ODE version of [[equation:save_globals_2]].
    !! This interface is consistent with [[eq_save_globals_2]].
      import :: ode_equation
      class(ode_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine ode_eq_save_globals_2

    subroutine ode_eq_load_globals ( this )
    !! The interface for the ODE version of [[equation:load_globals]].
    !! This interface is consistent with [[eq_load_globals]].
      import :: ode_equation
      class(ode_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine ode_eq_load_globals

    subroutine ode_eq_output ( this )
    !! The interface for the ODE version of [[equation:output]].
    !! This interface is consistent with [[eq_output]].
      import :: ode_equation
      class(ode_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine ode_eq_output
  end interface

  interface
    module subroutine ode_set_to_zero ( this, dest )
    !! The interface for the ODE version of [[equation:set_to_zero]].
    !! This interface is consistent with [[eq_set_to_zero_interface]].
      class(ode_equation), intent(inout) :: this
      !! The routine is called on this equation.
      integer(ip), intent(in) :: dest
      !! Can be either -1 (RHS), 0 (main) or between 1 and
      !! [[equation:ntmp(variable)]] (temporary).
    end subroutine ode_set_to_zero

    module subroutine ode_update_vars ( this, source, dest, source2, &
                                              scalar, scalar2 )
    !! The interface for the ODE version of [[equation:update_vars]].
    !! This interface is consistent with [[eq_update_vars_interface]].
      class(ode_equation), target, intent(inout) :: this
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
    end subroutine ode_update_vars
    
    module subroutine ode_print_data ( this )
    !! The interface for the ODE version of [[equation:print_data]].
    !! This interface is consistent with [[eq_print_data]].
      class(ode_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine ode_print_data
  end interface
end module ode_equations
