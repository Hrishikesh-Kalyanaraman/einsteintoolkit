module singular_observer

  use observers
  use effective_source, only: eff_source

  implicit none

  type, extends(observer) :: sing_observer
    integer(ip) :: nl
    real(wp), dimension(:,:), allocatable :: psi
    real(wp), dimension(:,:,:), allocatable :: dpsidt, dpsidr
    class(eff_source), pointer :: p
  contains
    procedure :: init => sobs_init
    procedure :: extract => sobs_extract
    procedure :: output => sobs_output
    final :: close_sing_observer
  end type sing_observer

  interface
    module subroutine sobs_init ( this, rad, coord, object )
      class(sing_observer), intent(inout) :: this
      real(wp), dimension(:), intent(in) :: rad
      type(rgf), intent(in) :: coord
      class(*), target, intent(in) :: object
    end subroutine sobs_init

    module subroutine sobs_extract ( this, effs )
      class(sing_observer), intent(inout) :: this
      class(eff_source), optional, pointer, intent(in) :: effs
      !! Optionally pass in an effective source object when regularization has
      !! to be done.
    end subroutine sobs_extract

    module subroutine sobs_output ( this )
      class(sing_observer), intent(inout) :: this
    end subroutine sobs_output

    module subroutine close_sing_observer ( this )
      type(sing_observer) :: this 
    end subroutine close_sing_observer

  end interface
end module singular_observer

