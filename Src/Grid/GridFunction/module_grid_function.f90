module grid_function
!! Module that defines the concept of a grid function and the interface of the
!! associated routines.
!!
!! The implementation is found in
!! [[submodule_grid_function_implementation.f90]].
  use kinds
  use element

  implicit none

  type, abstract :: gf
  !! An abstract grid function class.
    integer(ip) :: n
    !! The number of elements.
    character(:), allocatable :: vname
    !! The name of the grid function.
    integer(ip) :: io_id
    !! The file unit used for output of this grid function.
  end type gf
  
  type, extends(gf) :: rgf
  !! A real data instance of the abstract grid function class. Note this is
  !! not complete, as it has not been needed for evolution yet.
    type(element_rdata), dimension(:), allocatable :: elems
    !! A 1d array of the real element data class.
    contains
      procedure, non_overridable :: output => output_rgf
      !! Generic type bound procedure for output.
      final :: deallocate_rgf
      !! The finalizer.
  end type rgf

  type, extends(gf) :: cgf
  !! A complex data instance of the abstract grid function class.
    type(element_cdata), dimension(:), allocatable :: elems
    !! A 1d array of the complex element data class.
    contains
      procedure, non_overridable :: output => output_cgf
      !! Generic type bound procedure for output.
      procedure, non_overridable :: copy => copy_cgf
      !! Generic type bound procedure for copying the data from one complex
      !! grid function to another.
      procedure, non_overridable :: zero => zero_cgf
      !! Generic type bound procedure for setting a complex grid function
      !! to zero.
      procedure, non_overridable :: mult_sc => mult_sc_cgf
      !! Generic type bound procedure for multplying a complex grid function
      !! with a real scalar.
      procedure, non_overridable :: sc_mult_gf => sc_mult_gf_cgf
      !! Generic type bound procedure for multplying a complex grid function
      !! with a real scalar and storing the result in another complex
      !! grid function.
      procedure, non_overridable :: add_gf => add_gf_cgf
      !! Generic type bound procedure for adding 2 complex grid functions
      !! together and storing the result in the first one.
      procedure, non_overridable :: add_sc_mult_gf => add_sc_mult_gf_cgf
      !! Generic type bound procedure adding together a grid function and
      !! a scalar multiplying another grid function and storing the result
      !! in the first one.
      procedure, non_overridable :: mult_sc_add_sc_mult_gf => mult_sc_add_sc_mult_gf_cgf
      !! Generic type bound procedure for adding together a scalar multiplying
      !! a grid function with another scalar multiplying another grid function
      !! and storing the result in the first one.
      procedure, non_overridable :: gf1_plus_sc_mult_gf2 => gf1_plus_sc_mult_gf2_cgf
      !! Generic type bound procedure for storing the result of adding a
      !! second grid function with a scalar multiplying a third grid function
      !! in a grid function.
      procedure, non_overridable :: sc_mult_gf1_plus_sc_mult_gf2 => sc_mult_gf1_plus_sc_mult_gf2_cgf
      !! Generic type bound procedure for storing the result of adding a
      !! scalar multiplying a second grid function with another scalar
      !! multiplying a third grid function in a grid function.
      final :: deallocate_cgf
      !! The finalizer.
  end type cgf

  type, abstract :: gf_pointer
  !! An abstract class of a pointer to a grid function.
  end type gf_pointer

  type, extends(gf_pointer) :: rgf_pointer
  !! A real data type instance of a pointer to a grid function class.
    class(rgf), pointer :: p
    !! The pointer to a real grid function.
  end type rgf_pointer

  type, extends(gf_pointer) :: cgf_pointer
  !! A complex data type instance of a pointer to a grid function class.
    class(cgf), pointer :: p
    !! The pointer to a complex grid function.
  end type cgf_pointer

  type, abstract :: gfb
  !! An abstract class of a grid function of element boundary data.
    integer(ip) :: n
    !! The number of elements in the grid function.
    character(:), allocatable :: vname
    !! The name of the boundary data.
    integer(ip) :: iob_id
    !! The file unit used for output of this object.
  end type gfb

  type, extends(gfb) :: igfb
  !! An integer data type instance of a boundary grid funcion class.
    type(element_boundary_idata), dimension(:), allocatable :: elems
    !! A 1d array of element_boundary_idata type.
    contains
      procedure, non_overridable :: output => output_igfb
      !! Generic type bound procedure for output.
      final :: deallocate_igfb
      !! The finalizer.
  end type igfb

  type, extends(gfb) :: rgfb
  !! A real data type instance of a boundary grid funcion class.
    type(element_boundary_rdata), dimension(:), allocatable :: elems
    !! A 1d array of element_boundary_rdata type.
    contains
      procedure, non_overridable :: output => output_rgfb
      !! Generic type bound procedure for output.
      final :: deallocate_rgfb
      !! The finalizer.
  end type rgfb

  type, extends(gfb) :: cgfb
  !! A complex data type instance of a boundary grid funcion class.
    type(element_boundary_cdata), dimension(:), allocatable :: elems
    !! A 1d array of element_boundary_cdata type.
    contains
      procedure, non_overridable :: output => output_cgfb
      !! Generic type bound procedure for output.
      final :: deallocate_cgfb
      !! The finalizer.
  end type cgfb

  interface rgf
  !! The constructor for a real data type grid function.
    module procedure init_rgf
  end interface rgf

  interface cgf
  !! The constructor for a complex data type grid function.
    module procedure init_cgf
  end interface cgf

  interface igfb
  !! The constructor for an integer data type boundary grid function.
    module procedure init_igfb
  end interface igfb

  interface rgfb
  !! The constructor for a real data type boundary grid function.
    module procedure init_rgfb
  end interface rgfb

  interface cgfb
  !! The constructor for a complex data type boundary grid function.
    module procedure init_cgfb
  end interface cgfb

  interface
    module function init_rgf ( n, order, var_name )
    !! The interface for a constructor for a real data type grid function.
      type(rgf) :: init_rgf
      !! The return type has to be of type rgf.
      character(len=*), intent(in) :: var_name
      !! The name to be assigned to this grid function.
      integer(ip), intent(in) :: n
      !! The number of elements in this grid function.
      integer(ip), intent(in) :: order
      !! The order of the elements in this grid function.
    end function init_rgf

    module subroutine output_rgf ( this, coord )
    !! The interface for an output routine for a real data type grid function.
      class(rgf), intent(inout) :: this
      !! Has to be called with an rgf class.
      type(rgf), intent(in) :: coord
      !! A real grid function that contains the coordinates.
    end subroutine output_rgf

    module subroutine deallocate_rgf ( this )
    !! The interface for a finalize for a real data type grid function.
      type(rgf) :: this
      !! Has to be called with an rgf class.
    end subroutine deallocate_rgf

    module function init_cgf ( n, order, var_name )
    !! The interface for a constructor for a complex data type grid function.
      type(cgf) :: init_cgf
      !! The return type has to be of type cgf.
      character(len=*), intent(in) :: var_name
      !! The name to be assigned to this grid function.
      integer(ip), intent(in) :: n
      !! The number of elements in this grid function.
      integer(ip), intent(in) :: order
      !! The order of the elements in this grid function.
    end function init_cgf

    module subroutine output_cgf ( this, coord )
    !! The interface for an output routine for a complex data type grid
    !! function.
      class(cgf), intent(inout) :: this
      !! Has to be called with an cgf class.
      type(rgf), intent(in) :: coord
      !! A real grid function that contains the coordinates.
    end subroutine output_cgf

    module subroutine copy_cgf ( this, gf )
    !! The interface for a procedure for copying the data from one complex
    !! grid function to another.
    !!
    !! this%elems(:)%var(:) = gf%elems(:)%var(:)
      class(cgf), intent(inout) :: this
      !! The routine is called on this object.
      class(cgf), intent(in) :: gf 
      !! The grid function to copy.
    end subroutine copy_cgf

    module subroutine zero_cgf ( this )
      !! The interface for a procedure for setting a complex grid function
      !! to zero.
      !!
      !! this%elems(:)%var(:) = 0
      class(cgf), intent(inout) :: this
      !! The routine is called on this object.
    end subroutine zero_cgf

    module subroutine mult_sc_cgf ( this, scalar )
    !! The interface for a procedure for multplying a complex grid function
    !! with a real scalar.
    !!
    !! this%elems(:)%var(:) = scalar*this%elems(:)%var(:)
      class(cgf), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: scalar 
      !! The scalar to multiply with.
    end subroutine mult_sc_cgf

    module subroutine sc_mult_gf_cgf ( this, scalar, gf )
    !! The interface for a procedure for multplying a complex grid function
    !! with a real scalar and storing the result in another complex
    !! grid function.
    !!
    !! this%elems(:)%var(:) = scalar*gf%elems(:)%var(:)
      class(cgf), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: scalar 
      !! The scalar to multiply with.
      class(cgf), intent(in) :: gf
      !! The grid function to multiply with.
    end subroutine sc_mult_gf_cgf

    module subroutine add_gf_cgf ( this, gf )
    !! The interface for a procedure for adding 2 complex grid functions
    !! together and storing the result in the first one.
    !!
    !! this%elems(:)%var(:) = this%elems(:)%var(:)+gf%elems(:)%var(:)

      class(cgf), intent(inout) :: this
      !! The routine is called on this object.
      class(cgf), intent(in) :: gf
      !! The grid function to add.
    end subroutine add_gf_cgf

    module subroutine add_sc_mult_gf_cgf ( this, scalar, gf )
    !! The interface for a procedure adding together a grid function and
    !! a scalar multiplying another grid function and storing the result
    !! in the first one.
    !!
    !! this%elems(:)%var(:) = this%elems(:)%var(:)+scalar*gf%elems(:)%var(:)
      class(cgf), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: scalar
      !! The scalar to multiply with.
      class(cgf), intent(in) :: gf
      !! The second grid function.
    end subroutine add_sc_mult_gf_cgf

    module subroutine mult_sc_add_sc_mult_gf_cgf ( this, scalar1, scalar2, gf )
    !! The interface for a procedure for adding together a scalar multiplying
    !! a grid function with another scalar multiplying another grid function
    !! and storing the result in the first one.
    !!
    !! this%elems(:)%var(:) = scalar1\*this%elems(:)%var(:)+scalar2\*gf%elems(:)%var(:)
      class(cgf), intent(inout) :: this
      !! The routine is called on this object.
      real(wp), intent(in) :: scalar1
      !! The first scalar.
      real(wp), intent(in) :: scalar2
      !! The second scalar.
      class(cgf), intent(in) :: gf
      !! The second grid function.
    end subroutine mult_sc_add_sc_mult_gf_cgf

    module subroutine gf1_plus_sc_mult_gf2_cgf ( this, gf1, scalar, gf2 )
    !! The interface for a procedure for storing the result of adding a
    !! second grid function with a scalar multiplying a third grid function
    !! in a grid function.
    !!
    !! this%elems(:)%var(:) = gf1%elems(:)%var(:)+scalar*gf2%elems(:)%var(:)
      class(cgf), intent(inout) :: this
      !! The routine is called on this object.
      class(cgf), intent(in) :: gf1
      !! The second grid function.
      class(cgf), intent(in) :: gf2
      !! The third grid function.
      real(wp), intent(in) :: scalar
      !! The scalar.
    end subroutine gf1_plus_sc_mult_gf2_cgf

    module subroutine sc_mult_gf1_plus_sc_mult_gf2_cgf ( this, gf1, scalar1, &
                                                         gf2, scalar2 )
    !! The interface for a procedure for storing the result of adding a
    !! scalar multiplying a second grid function with another scalar
    !! multiplying a third grid function in a grid function.
    !!
    !! this%elems(:)%var(:) = scalar1\*gf1%elems(:)%var(:)+scalar2\*gf2%elems(:)%var(:)
      class(cgf), intent(inout) :: this
      !! The routine is called on this object.
      class(cgf), intent(in) :: gf1
      !! The second grid function.
      class(cgf), intent(in) :: gf2
      !! The third grid function.
      real(wp), intent(in) :: scalar1
      !! The first scalar.
      real(wp), intent(in) :: scalar2
      !! The second scalar.
    end subroutine sc_mult_gf1_plus_sc_mult_gf2_cgf

    module subroutine deallocate_cgf ( this )
    !! The interface for a finalizer for a complex data type grid function.
      type(cgf) :: this
      !! The object to finalize
    end subroutine deallocate_cgf

    module function init_igfb ( n, var_name )
    !! The interface for a constructor for an integer data type boundary grid
    !! function.
      type(igfb) :: init_igfb
      !! The object to construct.
      character(len=*), intent(in) :: var_name
      !! The name of the object.
      integer(ip), intent(in) :: n
      !! The number of elements.
    end function init_igfb

    module subroutine output_igfb ( this, coord )
    !! The interface for an output routine for an integer data type boundary
    !!grid function.
      class(igfb), intent(inout) :: this
      !! The routine is called on this object.
      type(rgf), intent(in) :: coord
      !! A real grid function that contains the coordinates.
    end subroutine output_igfb

    module subroutine deallocate_igfb ( this )
    !! The interface for a finalizer for an integer data type boundary
    !! grid function.
      type(igfb) :: this
      !! The object to finalize.
    end subroutine deallocate_igfb

    module function init_rgfb ( n, var_name )
    !! The interface for a constructor for a real data type boundary grid
    !! function.
      type(rgfb) :: init_rgfb
      !! The object to construct.
      character(len=*), intent(in) :: var_name
      !! The name of the object.
      integer(ip), intent(in) :: n
      !! The number of elements.
    end function init_rgfb

    module subroutine output_rgfb ( this, coord )
    !! The interface for an output routine for a real data type boundary
    !!grid function.
      class(rgfb), intent(inout) :: this
      !! The routine is called on this object.
      type(rgf), intent(in) :: coord
      !! A real grid function that contains the coordinates.
    end subroutine output_rgfb

    module subroutine deallocate_rgfb ( this )
    !! The interface for a finalizer for a real data type boundary
    !! grid function.
      type(rgfb) :: this
      !! The object to finalize.
    end subroutine deallocate_rgfb

    module function init_cgfb ( n, var_name )
    !! The interface for a constructor for a complex data type boundary grid
    !! function.
      type(cgfb) :: init_cgfb
      !! The object to construct.
      character(len=*), intent(in) :: var_name
      !! The name of the object.
      integer(ip), intent(in) :: n
      !! The number of elements.
    end function init_cgfb

    module subroutine output_cgfb ( this, coord )
    !! The interface for an output routine for a complex data type boundary
    !!grid function.
      class(cgfb), intent(inout) :: this
      !! The routine is called on this object.
      type(rgf), intent(in) :: coord
      !! A real grid function that contains the coordinates.
    end subroutine output_cgfb

    module subroutine deallocate_cgfb ( this )
    !! The interface for a finalizer for a complex data type boundary
    !! grid function.
      type(cgfb) :: this
      !! The object to finalize.
    end subroutine deallocate_cgfb

  end interface

end module grid_function
