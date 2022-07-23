module analytic_circular_orbit

  use kinds
  use equations

  implicit none

  type, extends(equation) :: circular_orbit
    real(wp) :: r
    real(wp) :: omega
    real(wp) :: amp
    real(wp) :: sigma
    real(wp) :: t0
    real(wp) :: En
    real(wp) :: Lz
    real(wp) :: chiomega
    real(wp) :: phi_initial
    integer(ip) :: io_id
  contains
    procedure :: init => co_init
    procedure :: rhs => co_rhs
    procedure :: set_to_zero => co_set_to_zero
    procedure :: update_vars => co_update_vars
    procedure :: save_globals_1 => co_save_globals_1
    procedure :: save_globals_2 => co_save_globals_2
    procedure :: load_globals => co_load_globals
    procedure :: output => co_output
    procedure :: print_data => co_print_data
  end type circular_orbit

  interface
    module subroutine co_init ( this )
      class(circular_orbit), target, intent(inout) :: this
    end subroutine co_init

    module subroutine co_rhs ( this )
      class(circular_orbit), intent(inout) :: this
    end subroutine co_rhs

    module subroutine co_set_to_zero ( this, dest )
      class(circular_orbit), intent(inout) :: this
      integer(ip), intent(in) :: dest
    end subroutine co_set_to_zero

    module subroutine co_update_vars ( this, source, dest, source2, &
                                              scalar, scalar2 )
      class(circular_orbit), target, intent(inout) :: this
      integer(ip), intent(in) :: source, dest
      integer(ip), optional, intent(in) :: source2
      real(wp), optional, intent(in) :: scalar, scalar2
    end subroutine co_update_vars

    module subroutine co_save_globals_1 ( this )
      class(circular_orbit), intent(inout) :: this
    end subroutine co_save_globals_1

    module subroutine co_save_globals_2 ( this )
      class(circular_orbit), intent(inout) :: this
    end subroutine co_save_globals_2

    module subroutine co_load_globals ( this )
      class(circular_orbit), intent(inout) :: this
    end subroutine co_load_globals 

    module subroutine co_output ( this )
      class(circular_orbit), intent(inout) :: this
    end subroutine co_output 

    module subroutine co_print_data ( this )
      class(circular_orbit), intent(inout) :: this
    end subroutine co_print_data

  end interface

end module analytic_circular_orbit
