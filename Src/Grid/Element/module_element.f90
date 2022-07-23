module element
!! Module that defines the concept of a DG element and various associated
!! routines.
!!
!! The implementation is found in [[submodule_element_implementation.f90]])
  use kinds

  type, abstract :: element_data
  !! An abstract class for the data in an element.
    integer(ip) :: order
    !! The order of the element.
  end type element_data

  type, extends(element_data) :: element_rdata
  !! A real data type instance of the abstract element data type.
    real(wp), dimension(:), allocatable :: var
    !! A real 1d array that will contain the data of the element.
    contains
      final :: deallocate_rdata
      !! The finalizer that will deallocate the 1d data array.
  end type element_rdata

  type, extends(element_data) :: element_cdata
  !! A complex data type instance of the abstract element data type.
    complex(wp), dimension(:), allocatable :: var
    !! A complex 1d array that will contain the data of the element.
    contains
      final :: deallocate_cdata
      !! The finalizer that will deallocate the 1d data array.
  end type element_cdata

  type, abstract :: element_boundary_data
  !! An abstract class for the boundary data in an element.
  end type element_boundary_data

  type, extends(element_boundary_data) :: element_boundary_idata
  !! An integer data type instance of the abstract element boundary data type.
    integer(ip), dimension(2) :: bvar
    !! An integer 1d array of size 2 that will contain the boundary data of
    !! the element.
  end type element_boundary_idata

  type, extends(element_boundary_data) :: element_boundary_rdata
  !! A real data type instance of the abstract element boundary data type.
    real(wp), dimension(2) :: bvar
    !! A real 1d array of size 2 that will contain the boundary data of
    !! the element.
  end type element_boundary_rdata

  type, extends(element_boundary_data) :: element_boundary_cdata
  !! A complex data type instance of the abstract element boundary data type.
    complex(wp), dimension(2) :: bvar
    !! A complex 1d array of size 2 that will contain the boundary data of
    !! the element.
  end type element_boundary_cdata

  interface element_rdata 
  !! The constructor for the real element data class.
    module procedure allocate_rdata
  end interface

  interface element_cdata 
  !! The constructor for the complex element data class.
    module procedure allocate_cdata
  end interface

  interface element_boundary_idata 
  !! The constructor for the integer element boundary data class.
    module procedure init_boundary_idata
  end interface

  interface element_boundary_rdata 
  !! The constructor for the real element boundary data class.
    module procedure init_boundary_rdata
  end interface

  interface element_boundary_cdata 
  !! The constructor for the complex element boundary data class.
    module procedure init_boundary_cdata
  end interface

  interface
    module function allocate_rdata ( order )
    !! The interface to the implementation of the constructor for the
    !! real element data class.
      type(element_rdata) :: allocate_rdata
      !! The return object has to be of the element_rdata class.
      integer(ip), intent(in) :: order
      !! The order of the element. The data array will be of size order+1.
    end function allocate_rdata

    module subroutine deallocate_rdata ( this )
    !! The interface to the implementation of the finalizer for the
    !! real element data class.
      type(element_rdata), intent(inout) :: this
      !! The argument has to be of the element_rdata class.
    end subroutine deallocate_rdata

    module function allocate_cdata ( order )
    !! The interface to the implementation of the constructor for the
    !! complex element data class.
      type(element_cdata) :: allocate_cdata
      !! The return object has to be of the element_cdata class.
      integer(ip), intent(in) :: order
      !! The order of the element. The data array will be of size order+1.
    end function allocate_cdata

    module subroutine deallocate_cdata ( this )
    !! The interface to the implementation of the finalizer for the
    !! complex element data class.
      type(element_cdata), intent(inout) :: this
      !! The argument has to be of the element_cdata class.
    end subroutine deallocate_cdata

    module function init_boundary_idata ( idata )
    !! The interface to the implementation of the constructor for the
    !! integer boundary element data class.
      type(element_boundary_idata) :: init_boundary_idata
      !! The return object has to be of the element_boundary_idata class.
      integer(ip), dimension(2), intent(in) :: idata
      !! The object gets initialized with the data in this array.
    end function init_boundary_idata

    module function init_boundary_rdata ( rdata )
    !! The interface to the implementation of the constructor for the
    !! real boundary element data class.
      type(element_boundary_rdata) :: init_boundary_rdata
      !! The return object has to be of the element_boundary_rdata class.
      real(wp), dimension(2), intent(in) :: rdata
      !! The object gets initialized with the data in this array.
    end function init_boundary_rdata

    module function init_boundary_cdata ( cdata )
    !! The interface to the implementation of the constructor for the
    !! complex boundary element data class.
      type(element_boundary_cdata) :: init_boundary_cdata
      !! The return object has to be of the element_boundary_cdata class.
      complex(wp), dimension(2), intent(in) :: cdata
      !! The object gets initialized with the data in this array.
    end function init_boundary_cdata
  end interface

end module element
