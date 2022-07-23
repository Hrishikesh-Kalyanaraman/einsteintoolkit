module observers
!! Module that defines the abstract interface of an observer class as well
!! as some simple observers for extracting real and complex data values from
!! grid functions.
!!
!! The implementation is found in [[submodule_observers_implementation.f90]].
  use kinds
  use grid_function
  use effective_source

  implicit none

  type, abstract :: observer
  !! An abstract observer interface.
    integer(ip) :: nradii
    !! The number of radii the observer should observe at.
    integer(ip) :: ioo_id
    !! The file unit number this observer should use for output.
    character(:), allocatable :: vname
    !! The name of the observer.
    real(wp), dimension(:), allocatable :: radii
    !! A 1d-array containing the radii that the observer should observe at.
    !! On allocation the size is [[observer:nradii]].
    integer(ip), dimension(:), allocatable :: elem_index
    !! A 1d-array containing the index of all the elements that contains
    !! [[observer:radii]].
    integer(ip), dimension(:), allocatable :: node_index
    !! A 1d-array containing the node index within all the elements that
    !! contains [[observer:radii]].
  contains
    procedure(obs_init_interface), deferred, pass :: init
    !! The initialization routine. Implementation is deferred to the
    !! derived class that actually implements an observer.
    procedure(obs_extract_interface), deferred, pass :: extract
    !! The extraction routine that performs the observation. Implementation is
    !! deferred to the derived class that actually implements an observer.
    procedure(obs_output_interface), deferred, pass :: output
    !! The output routine. Implementation is deferred to the derived class
    !! that actually implements an observer.
  end type observer

  abstract interface
    subroutine obs_init_interface ( this, rad, coord, object )
    !! The initialization routine.
      import :: observer, wp, rgf
      class(observer), intent(inout) :: this
      !! The observer that is being initialized.
      real(wp), dimension(:), intent(in) :: rad
      !! A 1d-array containing the radii where the observations have to be
      !! performed.
      type(rgf), intent(in) :: coord
      !! A real grid function containing the coordinates of the grid.
      class(*), target, intent(in) :: object
      !! The object on which observations have to be done.
    end subroutine obs_init_interface

    subroutine obs_extract_interface ( this, effs )
    !! The extraction routine.
      import :: observer, eff_source
      class(observer), intent(inout) :: this
      !! The routine is called on this observer.
      class(eff_source), optional, pointer, intent(in) :: effs
      !! Optionally pass in an effective source object when regularization has
      !! to be done.
    end subroutine obs_extract_interface

    subroutine obs_output_interface ( this )
    !! The output routine.
      import :: observer, wp
      class(observer), intent(inout) :: this
      !! The routine is called on this observer.
    end subroutine obs_output_interface
  end interface
 
  type, extends(observer) :: robserver
  !! An observer class that extracts information from real grid functions.
    type(rgf), pointer :: p
    !! A pointer to the real grid function that observations will be performed
    !! on.
    real(wp), dimension(:), allocatable :: extract_data
    !! A 1d-array that will hold the extracted data. On allocation the size is
    !! [[robserver:nradii]].
  contains
    procedure :: init => robs_init
    !! The initialization routine is provided by [[robs_init]].
    procedure :: extract => robs_extract
    !! The extraction routine is provided by [[robs_extract]].
    procedure :: output => robs_output
    !! The output routine is provided by [[robs_output]].
    final :: close_robserver
    !! The finalizer.
  end type robserver

  type, extends(observer) :: cobserver
  !! An observer class that extracts information from complex grid functions.
    type(cgf), pointer :: p
    !! A pointer to the complex grid function that observations will be
    !! performed on.
    complex(wp), dimension(:), allocatable :: extract_data
    !! A 1d-array that will hold the extracted data. On allocation the size is
    !! [[cobserver:nradii]].
  contains
    procedure :: init => cobs_init
    !! The initialization routine is provided by [[cobs_init]].
    procedure :: extract => cobs_extract
    !! The extraction routine is provided by [[cobs_extract]].
    procedure :: output => cobs_output
    !! The output routine is provided by [[cobs_output]].
    final :: close_cobserver
    !! The finalizer.
  end type cobserver

  interface
    module subroutine robs_init ( this, rad, coord, object )
    !! The initialization routine for the real grid function observer.
      class(robserver), intent(inout) :: this
      !! The observer that is being initialized.
      real(wp), dimension(:), intent(in) :: rad
      !! A 1d array of real values containing the radii for which observations
      !! will be performed.
      type(rgf), intent(in) :: coord
      !! A real grid function containing the coordinates.
      class(*), target, intent(in) :: object
      !! The real grid function for which observations will be performed.
    end subroutine robs_init

    module subroutine robs_extract ( this, effs )
    !! The extraction routine for the real grid function observer.
      class(robserver), intent(inout) :: this
      !! The routine is called on this observer.
      class(eff_source), optional, pointer, intent(in) :: effs
      !! Optionally pass in an effective source object when regularization has
      !! to be done.
    end subroutine robs_extract

    module subroutine robs_output ( this )
    !! The output routine for the real grid function observer.
      class(robserver), intent(inout) :: this
      !! The routine is called on this observer.
    end subroutine robs_output

    module subroutine close_robserver ( this )
    !! The finalizer for the real grid function observer.
      type(robserver) :: this
      !! The real grid function observer to be finalized..
    end subroutine close_robserver

    module subroutine cobs_init ( this, rad, coord, object )
    !! The initialization routine for the complex grid function observer.
      class(cobserver), intent(inout) :: this
      !! The observer that is being initialized.
      real(wp), dimension(:), intent(in) :: rad
      !! A 1d array of real values containing the radii for which observations
      !! will be performed.
      type(rgf), intent(in) :: coord
      !! A real grid function containing the coordinates.
      class(*), target, intent(in) :: object
      !! The complex grid function for which observations will be performed.
    end subroutine cobs_init

    module subroutine cobs_extract ( this, effs )
    !! The extraction routine for the complex grid function observer.
      class(cobserver), intent(inout) :: this
      !! The routine is called on this observer.
      class(eff_source), optional, pointer, intent(in) :: effs
      !! Optionally pass in an effective source object when regularization has
      !! to be done.
    end subroutine cobs_extract

    module subroutine cobs_output ( this )
    !! The output routine for the complex grid function observer.
      class(cobserver), intent(inout) :: this
      !! The routine is called on this observer.
    end subroutine cobs_output

    module subroutine close_cobserver ( this )
    !! The finalizer for the complex grid function observer.
      type(cobserver) :: this
      !! The routine is called on this observer.
    end subroutine close_cobserver

    module subroutine find_indices ( rad, coord, elem_index, node_index )
    !! A helper routine used by initialization routines to locate the
    !! element and node indices corresponding to the provided radii
    !!
    !! Currently the provided radii has to be within a distance or \(10^{-12}\)
    !! of a node location of the grid. If not the routine will abort the run.
      real(wp), dimension(:), intent(in) :: rad
      !! 1d array containing the radii to locate.
      type(rgf), intent(in) :: coord
      !! A real grid function containing the coordinates.
      integer(ip), dimension(:), allocatable, intent(out) :: elem_index
      !! On output contains the indices of the elements containing the
      !! provided radii.
      integer(ip), dimension(:), allocatable, intent(out) :: node_index
      !! On output contains the indices of the nodes inside the elements
      !! containing the provided radii.
    end subroutine find_indices

  end interface
end module observers
