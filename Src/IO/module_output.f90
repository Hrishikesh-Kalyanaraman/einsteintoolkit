module output_base
!! Module that provides basic IO functionality for keeping track of the
!! next available file unit.

  use kinds

  integer(ip), private :: next_id = 20
  !! This is kept private in order to ensure that everything is consistent.
  !! The first available file unit is 20.

contains

  function next_available_io_id () result(io_id)
  !! Function that provides the next available file unit.

    integer(ip) :: io_id
    !! The returned value is not used for output for any other quantity.

!$OMP CRITICAL
    io_id = next_id
    next_id = next_id + 1
!$OMP END CRITICAL

  end function next_available_io_id


  subroutine release_io_id ( io_id )
  !! Function that releases the last assigned file unit.

    integer(ip), intent(inout) :: io_id
    !! The file unit to attempt to release. If this is not the last assigned
    !! file unit, the routine aborts the execution.

!$OMP CRITICAL
    if ( io_id == next_id-1 ) then
      next_id = next_id - 1
      close(io_id)
      print*,'io_id: ', io_id, ' has been released'
      io_id = -1
    else
      print*,'release_io_id called with invalid io_id'
      stop
    end if
!$OMP END CRITICAL
  end subroutine release_io_id

end module output_base
