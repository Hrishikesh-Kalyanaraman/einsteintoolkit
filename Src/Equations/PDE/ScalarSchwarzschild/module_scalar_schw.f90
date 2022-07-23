module scalar_schw
!! Module that defines an equation class for evolving the spherical harmonic
!! decomposed field of a scalar point charge moving in  a Schwarzschild
!! spacetime. From the starting point of the wave equation in Tortoise
!! coordinates, the class supports the use of hyperboloidal coordinates in
!! regions near the horizon and near future null infinity as well as time
!! dependent coordinates for the case of a particle on a non-circular orbit.
!! The class interfaces with an effective source provided as C++ routines.
!!
!! The wave equation is written in first order in time and space form with
!! the evolved variables being \(\Psi,
!! \Upsilon=\frac{\partial\Psi}{\partial t},
!! \Pi=\frac{\partial\Psi}{\partial\rho}\) where \(\rho\) is the radial
!! coordinate and \(t\) is the time.
!!
!! The implementation is found in [[submodule_scalar_schw_implementation.f90]].
  use kinds
  use DG_structures
  use grid_function
  use pde_equations
  use scalar_schw_eff
  use orbit_base
  use time_dependent_coordinate
  use iso_c_binding

  implicit none

  type, extends(cpde_equation) :: scal_schw
  !! A class that defines the evolution equations for a scalar field produced
  !! by a point charge in Schwarzschild coordinates decomposed into
  !! spherical harmonics.
    integer(c_int), dimension(:), allocatable :: ll
    !! A 1d array of c_int containing the l values of all the evolved modes.
    integer(c_int), dimension(:), allocatable :: mm
    !! A 1d array of c_int containing the m values of all the evolved modes.
    type(ref_element) :: refelem
    !! The [[ref_element]] contains the derivative matrix needed for
    !! approximating partial derivatives as well as the routines needed
    !! for calculating characteristic fluxes.
    type(cgf), dimension(:,:), allocatable :: eq_flux_data
    !! A 2d array of complex grid functions that stores flux data for the
    !! system.
    type(rgf), dimension(:), allocatable :: eq_coeffs
    !! A 1d array of real grid functions that stores the equation coefficients.
    type(rgf), dimension(:), allocatable :: eq_lcoeffs
    !! A 1d array of real grid functions that stores the l-dependent
    !! potential.
    type(rgfb), dimension(:), allocatable :: eq_lambda
    !! A 1d array of real boundary grid functions that stores the
    !! characteristic speeds at the element boundaries.
    type(rgfb), dimension(:,:), allocatable :: eq_s
    !! A 2d array of real boundary grid functions that stores the
    !! matrix used to convert from characteristic to evolved variables.
    type(rgfb), dimension(:,:), allocatable :: eq_sinv
    !! A 2d array of real boundary grid functions that stores the
    !! matrix used to convert from evolve to characteristic variables.
    type(rgf) :: r_schw
    !! A real grid function that stores the Schwarzschild radial coordinate.
    type(rgf) :: r_star
    !! A real grid function that stores the tortoise radial coordinate.
    type(scal_schw_eff) :: effs
    !! The effective source for a scalar point charge on generic orbits in a
    !! Schwarzschild spacetime.
    type(tdc) :: time_dep_coord
    !! A time dependent coordinate class for the scalar wave equation in
    !! Tortoise coordinates.
  contains
    procedure :: init => scal_schw_init
    !! The [[equation:init]] routine is provided by [[scal_schw_init]].
    procedure :: rhs => scal_schw_rhs
    !! The [[equation:rhs]] routine is provided by [[scal_schw_rhs]].
    procedure :: save_globals_1 => scal_schw_save_globals_1
    !! The [[equation:save_globals_1]] routine is provided by
    !! [[scal_schw_save_globals_1]].
    procedure :: save_globals_2 => scal_schw_save_globals_2
    !! The [[equation:save_globals_2]] routine is provided by
    !! [[scal_schw_save_globals_2]].
    procedure :: load_globals => scal_schw_load_globals
    !! The [[equation:load_globals]] routine is provided by
    !! [[scal_schw_load_globals]].
    procedure :: apply_filter => scal_schw_apply_filter
    !! The [[cpde_equation:apply_filter]] routine is provided by
    !! [[scal_schw_apply_filter]].
    generic :: flux => scal_schw_flux
    !! A routine to calculate fluxes is provided by [[scal_schw_flux]] under
    !! the generic name [[scal_schw:flux]].
    procedure :: scal_schw_flux
    !! Provides a routine to calculate fluxes.
    procedure :: read_all_modes
    !! Provides a routine to read in external initial data.
    procedure :: output_coords
    !! Provides a routine to output the coordinates in a format needed for
    !! a frequency domain python code that can provide initial data.
    procedure :: tortoise_to_hyperboloidal
    !! Provides a routine to transform from tortoise to hyperboloid that
    !! is used when reading in external data.
  end type scal_schw

  complex(wp), dimension(:,:), allocatable :: drho
  !! A 2d array of reals used to store the radial derivative of \(\Upsilon\)
  !! for all modes in a single element. The size when allocated is
  !! ([[parameters:order]]+1:[[cpde_equation:nmodes]]).
  complex(wp), dimension(:,:), allocatable :: dpi
  !! A 2d array of reals used to store the radial derivative of \(\Pi\)
  !! for all modes in a single element. The size when allocated is
  !! ([[parameters:order]]+1:[[cpde_equation:nmodes]]).
  complex(wp), dimension(:,:,:), allocatable :: flux_result
  !! A 3d array of complex used to store fluxes for all modes at the boundary
  !! of a single element. The size when allocated is
  !! ([[parameters:order]]+1:2:[[cpde_equation:nmodes]]).

  interface
    module subroutine scal_schw_init ( this )
    !! The interface for the [[scal_schw]] version of [[equation:init]].
    !! This interface is consistent with [[eq_init_interface]].
      class(scal_schw), target, intent(inout) :: this
      !! The equation that is being initialized.
    end subroutine scal_schw_init

    module subroutine scal_schw_rhs ( this )
    !! The interface for the [[scal_schw]] version of [[equation:rhs]]. This
    !! interface is consistent with [[eq_rhs_interface]].
      class(scal_schw), intent(inout) :: this
      !! The equation for which the RHS is calculated.
    end subroutine scal_schw_rhs

    module subroutine scal_schw_save_globals_1 ( this )
    !! The interface for the [[scal_schw]] version of
    !! [[equation:save_globals_1]]. This interface is consistent with
    !! [[eq_save_globals_1]].
      class(scal_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine scal_schw_save_globals_1

    module subroutine scal_schw_save_globals_2 ( this )
    !! The interface for the [[scal_schw]] version of
    !! [[equation:save_globals_2]]. This interface is consistent with
    !! [[eq_save_globals_2]].
      class(scal_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine scal_schw_save_globals_2

    module subroutine scal_schw_load_globals ( this )
    !! The interface for the [[scal_schw]] version of
    !! [[equation:load_globals]]. This interface is consistent with
    !! [[eq_load_globals]].
      class(scal_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine scal_schw_load_globals

    module subroutine scal_schw_apply_filter ( this )
    !! The interface for the [[scal_schw]] version of
    !! [[cpde_equation:apply_filter]]. 
      class(scal_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine scal_schw_apply_filter

    module subroutine scal_schw_flux ( this )
    !! The interface for the [[scal_schw_flux]] routine for calculating
    !! fluxes. Can be called under the generic name [[scal_schw:flux]].
      class(scal_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine scal_schw_flux

    module subroutine read_all_modes ( this )
    !! The interface for the [[read_all_modes]] routine that reads
    !! external initial data.
      class(scal_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine read_all_modes

    module subroutine output_coords ( this )
    !! The interface for the [[output_coords]] routine that writes
    !! out the coordinates for an external python frequency domain code
    !! that calculates initial data.
      class(scal_schw), intent(inout) :: this
      !! The routine is called on this equation.
    end subroutine output_coords

    module subroutine tortoise_to_hyperboloidal ( this, elem, node, inner, &
                                                  dpsidt, dpsidr )
    !! The interface for the [[tortoise_to_hyperboloidal]] routine that
    !! converts from Tortoise to hyperboloidal coordinates when reading
    !! in external initial data.
      class(scal_schw), intent(inout) :: this
      !! The routine is called on this equation.
      integer(ip), intent(in) :: elem
      !! The index of the element where the transformation is calculated.
      integer(ip), intent(in) :: node
      !! The index of the node within that element where the transformation is
      !! calculated.
      logical, intent(in) :: inner
      !! If .true. use the transformation for the inner region. If .false. use
      !! the transformation for the outer region.
      complex(wp), intent(inout) :: dpsidt
      !! On input holds the time derivative in Tortoise coordinates. On output
      !! holds the time derivative in hyperboloidal coordinates.
      complex(wp), intent(inout) :: dpsidr
      !! On input holds the radial derivative in Tortoise coordinates. On
      !! output holds the radial derivative in hyperboloidal coordinates.
    end subroutine tortoise_to_hyperboloidal

    module function convert_var_name ( var_name, mode )
    !! A helper function that produces unique file names by adding the
    !! mode number to a variable name.
      character(len=*), intent(in) :: var_name
      !! The input variable name.
      integer(ip), intent(in) :: mode
      !! the mode number.
      character(len=len(var_name)+5) :: convert_var_name
      !! The return value is the converted variable name.
    end function convert_var_name
  end interface

end module scalar_schw
