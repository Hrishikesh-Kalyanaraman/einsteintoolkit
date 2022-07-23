module world_tube
!! Module that defines a world-tube class.
!!
!! The implementation is found in [[submodule_world_tube_implementation.f90]].

  use kinds
  use grid_function

  implicit none

  type ::  wtube
  !! A world tube class.
    integer(ip) :: n
    !! The number of elements on the grid.
    integer(ip) :: wsize
    !! The number of elements that is inside the world-tube.
    type(rgf) :: win
    !! A real grid function that contains a window function used to define the
    !! world-tube.
    type(rgf) :: dwin
    !! A real grid function that contains the radial derivative of the window
    !! function.
    type(rgf) :: d2win
    !! A real grid function that contains the second radial derivative of the
    !! window function.
    type(igfb) :: boundary_info
    !! An integer boundary grid function that contains information about
    !! whether an element boundary coincides with the world-tube.
    integer(ip) :: windex1
    !! The element index of the left boundary of the world-tube.
    integer(ip) :: windex2
    !! The element index of the right boundary of the world-tube.
    contains
      procedure, non_overridable :: is_boundary
      !! Routine to decide whether a given element boundary is a world-tube
      !! boundary.
  end type wtube

  type(wtube) :: wt
  !! A world-tube object that can be made available by use association.

  interface wtube
    module procedure init_wtube
    !! The world-tube constructor.
  end interface wtube

  interface
    module function init_wtube ( )
    !! The interface for the world-tube constructor.
      type(wtube) :: init_wtube
      !! The world tube being constructed.
    end function init_wtube

    module function is_boundary ( this, n, dir )
    !! A routine to determine whether a given element boundary is a world-tube
    !! boundary.
      class(wtube), intent(in) :: this
      !! The routine is called on this world-tube object.
      integer(ip), intent(in) :: n
      !! The index of the element to check.
      integer(ip), intent(in) :: dir
      !! The direction of the boundary within the element: -1 for left.
      !! boundary, +1 for the right boundary.
      logical :: is_boundary
      !! On return, .true. if this is a world-tube boundary and .false. if it
      !! is not.
    end function is_boundary

  end interface

end module world_tube
