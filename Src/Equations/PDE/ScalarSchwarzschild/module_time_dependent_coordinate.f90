module time_dependent_coordinate
!! Module that defines the interface for a class that provides the routines
!! necessary for the use of time dependent coordinates. This is currently
!! tied to supporting the scalar wave equation on a Schwarzschild background
!! provided in [[scalar_schw]].
!!
!! Starting from Tortoise coordinates \((t,r_*)\) the equations are
!! transformed to time dependent coordinate \((\lambda,\xi)\) as descibed in
!! Field, Hesthaven & Lau, Class. Quant. Grav. 26 (2009) 165010.

!! The implementation is provided in
!! [[submodule_time_dependent_coordinate_implementation.f90]].

  use kinds
  use grid_function

  implicit none

  type :: tdc
  !! A class that defines a time dependent coordinate transformation for the
  !! scalar wave equation on s Schwarzschild spacetime.
    type(rgfb) :: dxdlambda_b
    !! A real boundary grid function that stores
    !! \(\frac{\partial r_*}{\partial \lambda}\) on the element boundaries.
    type(rgfb) :: dxdxi_b
    !! A real boundary grid function that stores
    !! \(\frac{\partial r_*}{\partial \xi}\) on the element boundaries.
    real(wp) :: maxspeed
    !! The maximum coordinate speed (needed for CFL timestep condition).
    type(rgf) :: dxdlambda
    !! A real grid function that stores
    !! \(\frac{\partial r_*}{\partial \lambda}\) everywhere.
    type(rgf) :: dxdxi
    !! A real grid function that stores
    !! \(\frac{\partial r_*}{\partial \xi}\) everywhere.
    type(rgf) :: d2xdlambda2
    !! A real grid function that stores
    !! \(\frac{\partial^2 r_*}{\partial \lambda^2}\) everywhere.
    type(rgf) :: d2xdxi2
    !! A real grid function that stores
    !! \(\frac{\partial^2 r_*}{\partial \xi^2}\) everywhere.
    type(rgf) :: d2xdlambdadxi
    !! A real grid function that stores
    !! \(\frac{\partial^2 r_*}{\partial \lambda \partial \xi}\) everywhere.
    type(rgf) :: rm2m
    !! A real grid function that stores \(r-2M\) everywhere.
  contains
    procedure :: init => tdc_init
    !! The initialization routine.
    procedure :: set_coefficients => tdc_set_coefficients
    !! The routine that sets the coefficients for the wave equation.
    procedure :: tdc_to_tortoise_cvec => tdc_tdc_to_tortoise_cvec
    !! Routine to transform from \(\lambda,\xi)\) to \((t,r_*)\) for complex
    !! vector arguments at the boundary of an element.
    procedure :: tdc_to_tortoise_cscal => tdc_tdc_to_tortoise_cscal
    !! Routine to transform from \(\lambda,\xi)\) to \((t,r_*)\) for complex
    !! scalar arguments at the boundary of an element.
    procedure :: tdc_to_tortoise_rscal => tdc_tdc_to_tortoise_rscal
    !! Routine to transform from \(\lambda,\xi)\) to \((t,r_*)\) for real
    !! scalar arguments at the boundary of an element.
    generic :: tdc_to_tortoise => tdc_to_tortoise_cvec, tdc_to_tortoise_cscal, &
                                  tdc_to_tortoise_rscal
    !! Generic name for routines that transform from \(\lambda,\xi)\) to
    !! \((t,r_*)\) at the boundary of an element.
    procedure :: tortoise_to_tdc_cvecb => tdc_tortoise_to_tdc_cvecb
    !! Routine to transform from \((t,r_*)\) to \(\lambda,\xi)\) for complex
    !! vector arguments at the boundary of an element.
    procedure :: tortoise_to_tdc_cscalb => tdc_tortoise_to_tdc_cscalb
    !! Routine to transform from \((t,r_*)\) to \(\lambda,\xi)\) for complex
    !! scalar arguments at the boundary of an element.
    procedure :: tortoise_to_tdc_cscal => tdc_tortoise_to_tdc_cscal
    !! Routine to transform from \((t,r_*)\) to \(\lambda,\xi)\) for complex
    !! scalar arguments at any node in any element.
    generic :: tortoise_to_tdc => tortoise_to_tdc_cvecb, &
                                  tortoise_to_tdc_cscalb, &
                                  tortoise_to_tdc_cscal
    !! Generic name for routines that transform from \((t,r_*)\) to
    !! \(\lambda,\xi)\).
  end type tdc

  interface
    module subroutine tdc_init ( this )
    !! The interface for the initialization routine.
      class(tdc), intent(inout) :: this
      !! This time dependent coordinate object is being initialized.
    end subroutine tdc_init

    module subroutine tdc_set_coefficients ( this, coeffs, lcoeffs, lambda, &
                                             s, sinv, rho, rstar, rschw, ll )
    !! The interface for the routine that sets the coefficient for the
    !! scalar wave equation in Schwarzschild spacetime.
      use iso_c_binding
      class(tdc), intent(inout) :: this
      !! The routine is called on this time dependent coordinate object.
      type(rgf), dimension(:), intent(inout) :: coeffs
      !! A 1d array of real grid functions containing the wave equation
      !! coefficients (should be [[scal_schw:eq_coeffs]]).
      !! On output it will contain \((c_{\xi\xi}, c_{\lambda\xi},
      !! c_{\xi}, c_{\lambda})\).
      type(rgf), dimension(:), intent(inout) :: lcoeffs
      !! A 1d array of real grid functions containing the \(\ell\) dependent
      !! coefficients (should be [[scal_schw:eq_lcoeffs]]).
      !! On output it will contain the potential for all modes.
      type(rgfb), dimension(:), intent(inout) :: lambda
      !! A 1d array of real boundary grid functions containing the
      !! characteristic speeds at the boundary of the elements (should be
      !! [[scal_schw:eq_lambda]]). On output it will be updated.
      type(rgfb), dimension(:,:), intent(inout) :: s
      !! A 2d array of real boundary grid functions containing the matrix used
      !! to convert from characteristic to evolved variables (should be
      !! [[scal_schw:eq_s]]). On output it will be updated.
      type(rgfb), dimension(:,:), intent(inout) :: sinv
      !! A 2d array of real boundary grid functions containing the matrix used
      !! to convert from evolved to characteristic variables (should be
      !! [[scal_schw:eq_sinv]]). On output it will be updated.
      type(rgf), intent(in) :: rho
      !! A real grid function containing the computational radial coordinate,
      !! \(\rho\).
      type(rgf), intent(inout) :: rstar
      !! A real grid function containing the tortoise radial coordinate,
      !! \(r_*\). On output this will be updated.
      type(rgf), intent(inout) :: rschw
      !! A real grid function containing the Schwarzschild radial coodrinate,
      !! \(r_{\mathrm{schw}}\). On output this will be updated.
      integer(c_int), dimension(:), intent(in) :: ll
      !! A 1d array of c_int containing the \(\ell\)-values of all the modes.
    end subroutine tdc_set_coefficients

    module subroutine tdc_tdc_to_tortoise_cvec ( this, ielem, dir, &
                                                 dudlambda, dudxi, &
                                                 dudt, dudrstar )
    !! Routine to transform from time dependent \((\lambda,\xi)\) to tortoise
    !! coordinates \((t,r_*)\) for complex vector input at element boundaries.
      class(tdc), intent(inout) :: this
      !! The routine is called on this time dependent coordinate object.
      integer(ip), intent(in) :: ielem
      !! The index of the element that contains the coordinate transformation
      !! information.
      integer(ip), intent(in) :: dir
      !! The index of the boundary within the element that contains the
      !! coordinate transformation. Left boundary: 1, right boundary: 2.
      complex(wp), dimension(:), intent(in) :: dudlambda
      !! The vector of \(\frac{\partial\Psi}{\partial\lambda}\) values to
      !! transform.
      complex(wp), dimension(:), intent(in) :: dudxi
      !! The vector of \(\frac{\partial\Psi}{\partial\xi}\) values to
      !! transform.
      complex(wp), dimension(:), intent(out) :: dudt
      !! On output contains the \(\frac{\partial\Psi}{\partial t}\) vector.
      complex(wp), dimension(:), intent(out) :: dudrstar
      !! On output contains the \(\frac{\partial\Psi}{\partial r_*}\) vector.
    end subroutine tdc_tdc_to_tortoise_cvec

    module subroutine tdc_tdc_to_tortoise_cscal ( this, ielem, dir, &
                                                  dudlambda, dudxi, &
                                                  dudt, dudrstar )
    !! Routine to transform from time dependent \((\lambda,\xi)\) to tortoise
    !! coordinates \((t,r_*)\) for complex scalar input at element boundaries.
      class(tdc), intent(inout) :: this
      !! The routine is called on this time dependent coordinate object.
      integer(ip), intent(in) :: ielem
      !! The index of the element that contains the coordinate transformation
      !! information.
      integer(ip), intent(in) :: dir
      !! The index of the boundary within the element that contains the
      !! coordinate transformation. Left boundary: 1, right boundary: 2.
      complex(wp), intent(in) :: dudlambda
      !! The value of \(\frac{\partial\Psi}{\partial\lambda}\) to transform.
      complex(wp), intent(in) :: dudxi
      !! The value of \(\frac{\partial\Psi}{\partial\xi}\) to transform.
      complex(wp), intent(out) :: dudt
      !! On output contains the \(\frac{\partial\Psi}{\partial t}\) value.
      complex(wp), intent(out) :: dudrstar
      !! On output contains the \(\frac{\partial\Psi}{\partial r_*}\) value.
    end subroutine tdc_tdc_to_tortoise_cscal

    module subroutine tdc_tdc_to_tortoise_rscal ( this, ielem, dir, &
                                                  dudlambda, dudxi, &
                                                  dudt, dudrstar )
    !! Routine to transform from time dependent \((\lambda,\xi)\) to tortoise
    !! coordinates \((t,r_*)\) for real scalar input at element boundaries.
      class(tdc), intent(inout) :: this
      !! The routine is called on this time dependent coordinate object.
      integer(ip), intent(in) :: ielem
      !! The index of the element that contains the coordinate transformation
      !! information.
      integer(ip), intent(in) :: dir
      !! The index of the boundary within the element that contains the
      !! coordinate transformation. Left boundary: 1, right boundary: 2.
      real(wp), intent(in) :: dudlambda
      !! The value of \(\frac{\partial\Psi}{\partial\lambda}\) to transform.
      real(wp), intent(in) :: dudxi
      !! The value of \(\frac{\partial\Psi}{\partial\xi}\) to transform.
      real(wp), intent(out) :: dudt
      !! On output contains the \(\frac{\partial\Psi}{\partial t}\) value.
      real(wp), intent(out) :: dudrstar
      !! On output contains the \(\frac{\partial\Psi}{\partial r_*}\) value.
    end subroutine tdc_tdc_to_tortoise_rscal

    module subroutine tdc_tortoise_to_tdc_cvecb ( this, ielem, dir, &
                                                 dudt, dudrstar, &
                                                 dudlambda, dudxi )
    !! Routine to transform from tortoise coordinates \((t,r_*)\) to time
    !! dependent \((\lambda,\xi)\) for complex vector input at element
    !! boundaries.
      class(tdc), intent(inout) :: this
      !! The routine is called on this time dependent coordinate object.
      integer(ip), intent(in) :: ielem
      !! The index of the element that contains the coordinate transformation
      !! information.
      integer(ip), intent(in) :: dir
      !! The index of the boundary within the element that contains the
      !! coordinate transformation. Left boundary: 1, right boundary: 2.
      complex(wp), dimension(:), intent(in) :: dudt
      !! The vector of \(\frac{\partial\Psi}{\partial t}\) to transform.
      complex(wp), dimension(:), intent(in) :: dudrstar
      !! The vector of \(\frac{\partial\Psi}{\partial r_*}\) to transform.
      complex(wp), dimension(:), intent(out) :: dudlambda
      !! On output contains the \(\frac{\partial\Psi}{\partial\lambda}\) vector.
      complex(wp), dimension(:), intent(out) :: dudxi
      !! On output contains the \(\frac{\partial\Psi}{\partial\xi}\) vector.
    end subroutine tdc_tortoise_to_tdc_cvecb

    module subroutine tdc_tortoise_to_tdc_cscalb ( this, ielem, dir, &
                                                  dudt, dudrstar, &
                                                  dudlambda, dudxi )
    !! Routine to transform from tortoise coordinates \((t,r_*)\) to time
    !! dependent \((\lambda,\xi)\) for complex scalar input at element
    !! boundaries.
      class(tdc), intent(inout) :: this
      !! The routine is called on this time dependent coordinate object.
      integer(ip), intent(in) :: ielem
      !! The index of the element that contains the coordinate transformation
      !! information.
      integer(ip), intent(in) :: dir
      !! The index of the boundary within the element that contains the
      !! coordinate transformation. Left boundary: 1, right boundary: 2.
      complex(wp), intent(in) :: dudt
      !! The value of \(\frac{\partial\Psi}{\partial t}\) to transform.
      complex(wp), intent(in) :: dudrstar
      !! The value of \(\frac{\partial\Psi}{\partial r_*}\) to transform.
      complex(wp), intent(out) :: dudlambda
      !! On output contains the \(\frac{\partial\Psi}{\partial\lambda}\) value.
      complex(wp), intent(out) :: dudxi
      !! On output contains the \(\frac{\partial\Psi}{\partial\xi}\) value.
    end subroutine tdc_tortoise_to_tdc_cscalb

    module subroutine tdc_tortoise_to_tdc_cscal ( this, elem, node, &
                                                  dpsidt, dpsidr )
    !! Routine to transform from tortoise coordinates \((t,r_*)\) to time
    !! dependent \((\lambda,\xi)\) for complex scalar input.
      class(tdc), intent(inout) :: this
      !! The routine is called on this time dependent coordinate object.
      integer(ip), intent(in) :: elem
      !! The element index of the point to transform.
      integer(ip), intent(in) :: node
      !! The node index of the point to transform.
      complex(wp), intent(inout) :: dpsidt
      !! On input contains: \(\frac{\partial\Psi}{\partial t}\). On output
      !! contains \(\frac{\partial\Psi}{\partial\lambda}\).
      complex(wp), intent(inout) :: dpsidr
      !! On input contains: \(\frac{\partial\Psi}{\partial r_*}\). On output
      !! contains \(\frac{\partial\Psi}{\partial\xi}\).
    end subroutine tdc_tortoise_to_tdc_cscal
  end interface

end module time_dependent_coordinate
