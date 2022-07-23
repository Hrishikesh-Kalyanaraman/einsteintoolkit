submodule(element) element_implementation
!! The implementation of the interfaces defined in [[element]].

contains

  module procedure allocate_rdata  

    implicit none

    integer :: allocation_status

!    print*,'Allocating real element'
    allocate_rdata%order = order
    if ( .not. allocated(allocate_rdata%var) ) then
      allocate ( allocate_rdata%var(order+1), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Allocation error for real element data'
        stop
      else
        allocate_rdata%var = 0.0_wp
      end if
    end if
  end procedure allocate_rdata


  module procedure deallocate_rdata

    implicit none

    integer :: allocation_status

!    print*,'Deallocating real element'
    if ( allocated ( this%var ) ) then
      deallocate( this%var, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for real element data'
        stop
      end if
    end if
  end procedure deallocate_rdata


  module procedure allocate_cdata  

    implicit none

    integer :: allocation_status

!    print*,'Allocating complex element'
    allocate_cdata%order = order
    if ( .not. allocated(allocate_cdata%var) ) then
      allocate ( allocate_cdata%var(order+1), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Allocation error for element complex data'
        stop
      else
        allocate_cdata%var = cmplx(0.0_wp,0.0_wp)
      end if
    end if
  end procedure allocate_cdata


  module procedure deallocate_cdata

    implicit none

    integer :: allocation_status

!    print*,'Deallocating vector element'
    if ( allocated ( this%var ) ) then
      deallocate( this%var, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for element complex data'
        stop
      end if
    end if
  end procedure deallocate_cdata


  module procedure init_boundary_idata

    implicit none

    init_boundary_idata%bvar = idata

  end procedure init_boundary_idata


  module procedure init_boundary_rdata

    implicit none

    init_boundary_rdata%bvar = rdata

  end procedure init_boundary_rdata


  module procedure init_boundary_cdata

    implicit none

    init_boundary_cdata%bvar = cdata

  end procedure init_boundary_cdata

end submodule element_implementation
