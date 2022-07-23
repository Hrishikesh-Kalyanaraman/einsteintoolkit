module pde_equations
!! Module that defines a PDE equation. It is derived from the [[equation]]
!! class and provides the routines that are common for all PDE equation
!! systems while deferring implementation of the routines that are
!! specific to a given PDE equation system.
!!
!! The implementation is found in
!! [[submodule_pde_equations_implementation.f90]].

  use kinds
  use equations
  use grid_function

  implicit none

  type, abstract, extends(equation) :: cpde_equation
  !! A class derived from [[equation]] specific for PDE equations.
    integer(ip) :: nmodes
    !! The number of spherical harmonic modes. Should be moved to the specific
    !! implementation of the equation as this class should not be limited
    !! to systems using spherical harmonic decomposition.
    integer(ip) :: nvars
    !! The number of variables per mode. Should probably be changed to be the
    !! total number of variables in order to make this independent of
    !! spherical harmonic decomposition.
    type(cgf), dimension(:,:), allocatable :: eq_data
    !! A 2d array of complex grid functions that contains the data variables
    !! for the equation system. On allocation the size is 
    !! ([[cpde_equation:nvars]]:[[cpde_equation:nmodes]]). Should probably be
    !! changed to a 1d array to make this independent of spherical harmonic
    !! decomposition.
    type(cgf), dimension(:,:), allocatable :: eq_rhs_data
    !! A 2d array of complex grid functions that contains the rhs variables
    !! for the equation system. On allocation the size is 
    !! ([[cpde_equation:nvars]]:[[cpde_equation:nmodes]]). Should probably be
    !! changed to a 1d array to make this independent of spherical harmonic
    !! decomposition.
    type(cgf), dimension(:,:,:), allocatable :: eq_tmp_data
    !! A 3d array of complex grid functions that contains the temporary
    !! variables needed by the time integrator. On allocation the size is
    !! ([[cpde_equation:nvars]]:[[cpde_equation:nmodes]]:[[equation:ntmp]]).
    !! Should probably be changed to a 2d array to make this independent of
    !! spherical harmonic decomposition.
    type(cgf_pointer), dimension(:,:,:), allocatable :: data_pointer
    !! A 3d array of pointers to complex grid functions that points to
    !! the RHS variables (:,:,-1), the data variables (:,:,0) and the
    !! temporary variables (:,:,1:[[equation:ntmp]]) for the equation system.
    !! Should probably be changed to a 2d array to make this independent of
    !! spherical harmonic decomposition.
  contains 
    procedure (cpde_eq_init_interface), deferred, pass :: init
    !! The [[equation:init]] routine is deferred.
    procedure (cpde_eq_rhs_interface), deferred, pass :: rhs
    !! The [[equation:rhs]] routine is deferred.
    procedure :: set_to_zero => cpde_set_to_zero
    !! The [[equation:set_to_zero]] routine is provided by
    !! [[cpde_set_to_zero]].
    procedure :: update_vars => cpde_update_vars
    !! The [[equation:update_vars]] routine is provided by
    !! [[cpde_update_vars]].
    procedure (cpde_eq_save_globals_1), deferred, pass :: save_globals_1
    !! The [[equation:save_globals_1]] routine is deferred.
    procedure (cpde_eq_save_globals_2), deferred, pass :: save_globals_2
    !! The [[equation:save_globals_2]] routine is deferred.
    procedure (cpde_eq_load_globals), deferred, pass :: load_globals
    !! The [[equation:load_globals]] routine is deferred.
    procedure :: output => cpde_output
    !! The [[equation:output]] routine is provided by [[cpde_output]].
    procedure :: print_data => cpde_print_data
    !! The [[equation:print_data]] routine is provided by [[cpde_print_data]].
    procedure (cpde_apply_filter_interface), deferred, pass :: apply_filter
    !! A routine unique to PDE equations that allows a filter to be applied
    !! to the evolved variables.
  end type cpde_equation

  abstract interface
    subroutine cpde_eq_init_interface ( this )
    !! The interface for the PDE version of [[equation:init]]. This interface
    !! is consistent with [[eq_init_interface]].
      import :: cpde_equation, ip
      class(cpde_equation), target, intent(inout) :: this
      !! The equation that is being initialized.
    end subroutine cpde_eq_init_interface

    subroutine cpde_eq_rhs_interface ( this )
    !! The interface for the PDE version of [[equation:rhs]]. This interface
    !! is consistent with [[eq_rhs_interface]].
      import :: cpde_equation
      class(cpde_equation), intent(inout) :: this
      !! The equation for which the RHS is calculated.
    end subroutine cpde_eq_rhs_interface

    subroutine cpde_eq_save_globals_1 ( this )
    !! The interface for the PDE version of [[equation:save_globals_1]].
    !! This interface is consistent with [[eq_save_globals_1]].
      import :: cpde_equation
      class(cpde_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine cpde_eq_save_globals_1

    subroutine cpde_eq_save_globals_2 ( this )
    !! The interface for the PDE version of [[equation:save_globals_2]].
    !! This interface is consistent with [[eq_save_globals_2]].
      import :: cpde_equation
      class(cpde_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine cpde_eq_save_globals_2

    subroutine cpde_eq_load_globals ( this )
    !! The interface for the PDE version of [[equation:load_globals]].
    !! This interface is consistent with [[eq_load_globals]].
      import :: cpde_equation
      class(cpde_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine cpde_eq_load_globals

    subroutine cpde_apply_filter_interface ( this )
    !! The interface for the [[cpde_equation:apply_filter]] routine.
      import :: cpde_equation
      class(cpde_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine cpde_apply_filter_interface
  end interface

  interface
    module subroutine cpde_set_to_zero ( this, dest )
    !! The interface for the PDE version of [[equation:set_to_zero]].
    !! This interface is consistent with [[eq_set_to_zero_interface]].
      class(cpde_equation), intent(inout) :: this
      !! The routine is called on this equation.
      integer(ip), intent(in) :: dest
      !! Can be either -1 (RHS), 0 (main) or between 1 and
      !! [[equation:ntmp(variable)]] (temporary).
    end subroutine cpde_set_to_zero

    module subroutine cpde_update_vars ( this, source, dest, source2, &
                                              scalar, scalar2 )
    !! The interface for the PDE version of [[equation:update_vars]].
    !! This interface is consistent with [[eq_update_vars_interface]].
      class(cpde_equation), target, intent(inout) :: this
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
    end subroutine cpde_update_vars

    module subroutine cpde_output ( this )
    !! The interface for the PDE version of [[equation:output]].
    !! This interface is consistent with [[eq_output]].
      class(cpde_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine cpde_output

    module subroutine cpde_print_data ( this )
    !! The interface for the PDE version of [[equation:print_data]].
    !! This interface is consistent with [[eq_print_data]].
      class(cpde_equation), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine cpde_print_data
  end interface

end module pde_equations
