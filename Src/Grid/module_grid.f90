module grid
!! Module that defines the variables used to set up a grid for self-force
!! calculations. This also includes an initialization routine.

  use kinds
  use element
  use grid_function

  implicit none

  real(wp) :: rho_min
  !! The minimal coordinate point on the grid. With hyperboloidal coordinates
  !! this is the location of the horizon. Gets initialized to the value of the
  !! parameter [[parameters:Sminus]].
  real(wp) :: rho_max
  !! The maximal coordinate point on the grid. With hyperboloidal coordinates
  !! this is the location of scri+. Gets initialized so that the particle
  !! location is the center of the grid.
  real(wp) :: rho_particle
  !! The coordinate location of the center of the grid. When a particle is
  !! used it gets initialized to the tortoise coordinate location of the
  !! midpoint of the orbit. When a particle is not used it gets initialized
  !! to the value of the parameter [[parameters:r_center]] converted to
  !! tortoise coordinates.
  real(wp) :: Tminus
  !! The tortoise coordinate location of the inner boundary of the region
  !! where time dependent coordinates are used.
  real(wp) :: Tplus
  !! The tortoise coordinate location of the outer boundary of the region
  !! where time dependent coordinates are used.
  real(wp) :: delta_rho_min
  !! The smalles coordinate distance between any two DG nodes on the grid.
  type(rgf) :: rho
  !! A real grid function containing the computational coordinate,
  !!\(\rho\).
  type(rgf) :: drhodr
  !! A real grid function containing \(\frac{d\rho}{dr}\) where
  !! \(r\) is the node coordinates in the reference element (\(r\in[-1,1]\)).
  type(rgf) :: drdrho
  !! A real grid function containing \(\frac{dr}{d\rho}\).
  integer(ip), dimension(2) :: Tminus_ind
  !! A size 2 integer array containing the element indices of the two
  !! elements that border the lower boundery of the time dependent coordinate
  !! region.
  integer(ip), dimension(2) :: Tplus_ind
  !! A size 2 integer array containing the element indices of the two
  !! elements that border the upper boundery of the time dependent coordinate
  !! region.
  integer(ip), dimension(2) :: particle_element
  !! A size 2 integer array containing the element indices of the two
  !! elements that has the particle on their boundary.
  integer(ip), dimension(2) :: particle_node
  !! A size 2 integer array contining the node indeces of the particle
  !! location in the two elements that contain the particle.

  interface
    module subroutine init_grid_coordinates ( relem  )
    !! The routine that initializes the variables describing the grid.
    !!
    !! All values used to set up the grid are taken from the input
    !! parameters defined in [[parameters]].
      use DG_structures

      type(ref_element), intent(in) :: relem
      !! The reference element of the DG elements that is the building block
      !! for grid functions.
    end subroutine 
  end interface

end module grid
