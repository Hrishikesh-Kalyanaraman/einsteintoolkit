module self_force_observer
!! Module that defines the interface for a self-force observer class for a
!! scalar charge in orbit around a Schwarzschild black hole, i.e. it provides
!! extraction of the self-force from the variables in [[scalar_schw]]).
!!
!! The implementation is found in
!! [[submodule_self_force_observer_implementation.f90]]

  use scalar_schw
  use observers

  implicit none

  type, extends(observer) :: sf_observer
  !! A class that defines an observer of the self-force for a scalar charge
  !! in orbit around a Schwarzschild black hole.
    type(scal_schw), pointer :: p
    !! A pointer to the [[scal_schw]] equation that provides the data.
    real(wp), dimension(:), allocatable :: fl
    !! The \(\ell\)-modes of the regular field at the particle location.
    real(wp), dimension(:), allocatable :: ftl
    !! The \(\ell\)-modes of the time derivative of the regular field at the
    !! particle location.
    real(wp), dimension(:), allocatable :: fphil
    !! The \(\ell\)-modes of the \(\phi\) derivative of the regular field at
    !! the particle location.
    real(wp), dimension(:), allocatable :: frl
    !! The \(\ell\)-modes of the radial derivative of the regular field at the
    !! particle location.
  contains
    procedure :: init => sf_init
    !! The [[observer:init]] routine is provided by [[sf_init]]
    procedure :: extract => sf_extract
    !! The [[observer:extract]] routine is provided by [[sf_extract]]
    procedure :: output => sf_output 
    !! The [[observer:output]] routine is provided by [[sf_output]]
    final :: close_sf_observer
    !! The finalizer.
  end type sf_observer

  type(sf_observer) :: sfobs

  interface
    module subroutine sf_init ( this, rad, coord, object )
    !! The interface for the [[sf_observer]] version of [[observer:init]].
    !! This interface is consistent with [[obs_init_interface]].
      class(sf_observer), intent(inout) :: this
      !! The self-force observer that is being initialized.
      real(wp), dimension(:), intent(in) :: rad
      !! A 1d-array containing the radii where the observations have to be
      !! performed. Obviously it only makes sense to pass in the particle
      !! location.
      type(rgf), intent(in) :: coord
      !! A real grid function containing the coordinates of the grid.
      class(*), target, intent(in) :: object
      !! The object on which observations have to be done. This has to be of
      !! type [[scal_schw]].
    end subroutine sf_init

    module subroutine sf_extract ( this, effs )
    !! The interface for the [[sf_observer]] version of [[observer:extract]].
    !! This interface is consistent with [[obs_extract_interface]].
      class(sf_observer), intent(inout) :: this
      !! The routine is called on this self-force observer.
      class(eff_source), optional, pointer, intent(in) :: effs
      !! Optionally pass in an effective source object when regularization has
      !! to be done.
    end subroutine sf_extract

    module subroutine sf_output ( this )
    !! The interface for the [[sf_output]] version of [[observer:output]].
    !! This interface is consistent with [[obs_output_interface]].
      class(sf_observer), intent(inout) :: this
      !! The routine is called on this self-force observer.
    end subroutine sf_output

    module subroutine close_sf_observer ( this )
    !! The interface for the finalizer.
      type(sf_observer) :: this
      !! The self-force observer to be finalized.
    end subroutine close_sf_observer
  end interface
end module self_force_observer
