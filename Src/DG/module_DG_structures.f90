module DG_structures
!! Module that defines the basic DG reference element type and the interface
!! to it's constructor and finalizer as well as routines for calculating
!! characteristic fluxes.
!!
!! The implementation (found in [[submodule_DG_implementation.f90]]) is based
!! on 'Nodal Discontinuous Galerkin Methods' by Hesthaven and Warburton.

  use kinds

  implicit none

  type ref_element
  !! The reference element.
    integer(ip) :: n 
    !! Order of this element.
    real(wp), dimension(:), allocatable :: r
    !! Node location within this element.
    real(wp), dimension(:), allocatable :: w
    !! Weights for integration.
    real(wp), dimension(:,:), allocatable :: v
    !! Vandermonde matrix for this element.
    real(wp), dimension(:,:), allocatable :: dr
    !! Derivative matrix for this element.
    real(wp), dimension(:,:), allocatable :: lift
    !! Lift matrix for this element.
    real(wp), dimension(:,:), allocatable :: filter
    !! An optional smoothing filter.
    logical :: have_filter
    !! Set to .true. if a filter has been initialized.

    contains
      generic :: characteristic_flux => char_flux_real, &
                                        char_flux_complex
      !! Use characteristic_flux for either char_flux_real or
      !! char_flux_complex.
      procedure :: char_flux_real
      !! The real version of the characteristic flux routine.
      procedure :: char_flux_complex
      !! The complex version of the characteristic flux routine.
      final :: deallocate_ref_element
      !! The finalizer.
  end type ref_element

  interface ref_element
    module procedure init_ref_element
    !! The reference element constructor.
  end interface ref_element

  interface

    module function init_ref_element ( order, sorder, nc )
      !! Interface for the reference element constructor.
      type(ref_element) :: init_ref_element
      !! The return type has to be a reference element.
      integer(ip), intent(in) :: order
      !! The order of the reference element.
      integer(ip), optional :: sorder
      !! The order of a smoothing filter. Both this and nc have to be
      !! present before a filter is constructed.
      integer(ip), optional :: nc
      !! No filtering of orders below this value. Both this and sorder have
      !! to be present before a filter is constructed.
    end function init_ref_element

    module subroutine deallocate_ref_element ( refel )
      !! Interface for the finalizer.
      type(ref_element) :: refel
      !! The argument has to be a reference element type.
    end subroutine deallocate_ref_element
    
    module function char_flux_real ( this, nvar, order, uint, uext, flux, &
                                     lambda, s, sinv, debug_output )
      !! Interface for the characteristic flux routine for real data types.
      class(ref_element), intent(in) :: this
      !! Has to be a reference element type.
      integer(ip), intent(in) :: nvar
      !! The number of variables the characteristic flux has to be computed for.
      integer(ip), intent(in) :: order
      !! The order of the reference element. Used to declare the size of the
      !! return array.
      real(wp), dimension(2,nvar), intent(in) :: uint
      !! The boundary data internal to this element.
      real(wp), dimension(2,nvar), intent(in) :: uext
      !! The boundary data external to this element.
      real(wp), dimension(2,nvar), intent(in) :: flux
      !! The boundary fluxes.
      real(wp), dimension(2,nvar), intent(in) :: lambda
      !! The boundary characteristic speeds.
      real(wp), dimension(2,nvar,nvar), intent(in) :: s
      !! The matrix that converts from characteristic to evolved variables.
      real(wp), dimension(2,nvar,nvar), intent(in) :: sinv
      !! The matrix that converts from evolved to characteristic variables.
      logical, intent(in) :: debug_output
      !! Currently ignored for the real version of this routine.
      real(wp), dimension(order+1,nvar) :: char_flux_real
      !! At return the numerical flux has been lifted to the whole element.
    end function char_flux_real

    module function char_flux_complex ( this, nvar, order, uint, uext, flux, &
                                        lambda, s, sinv, debug_output )
      !! Interface for the characteristic flux routine for complex data types.
      class(ref_element), intent(in) :: this
      !! Has to be a reference element type.
      integer(ip), intent(in) :: nvar
      !! The number of variables the characteristic flux has to be computed for.
      integer(ip), intent(in) :: order
      !! The order of the reference element. Used to declare the size of the
      !! return array.
      complex(wp), dimension(2,nvar), intent(in) :: uint
      !! The boundary data internal to this element.
      complex(wp), dimension(2,nvar), intent(in) :: uext
      !! The boundary data external to this element.
      complex(wp), dimension(2,nvar), intent(in) :: flux
      !! The boundary fluxes.
      real(wp), dimension(2,nvar), intent(in) :: lambda
      !! The boundary characteristic speeds.
      real(wp), dimension(2,nvar,nvar), intent(in) :: s
      !! The matrix that converts from characteristic to evolved variables.
      real(wp), dimension(2,nvar,nvar), intent(in) :: sinv
      !! The matrix that converts from evolved to characteristic variables.
      logical, intent(in) :: debug_output
      !! If .true., produce debug output
      complex(wp), dimension(order+1,nvar) :: char_flux_complex
      !! At return the numerical flux has been lifted to the whole element.
    end function char_flux_complex

  end interface
  
end module DG_structures
