submodule(world_tube) world_tube_implementation
!! Implementation of the interfaces in [[world_tube]].

contains

  module procedure init_wtube

    use parameters, only : n_elems, order, world_tube_width
    implicit none

    integer(ip) :: i, j

    integer :: allocation_status

    if (world_tube_width<0) then
      print*,'Initialization of world tube called with world_tube_width < 0'
      stop
    end if

    init_wtube%wsize = 2*world_tube_width

    init_wtube%n = n_elems

    init_wtube%win = rgf ( n_elems, order, 'win' )
    init_wtube%dwin = rgf ( n_elems, order, 'dwin' )
    init_wtube%d2win = rgf ( n_elems, order, 'd2win' )

    init_wtube%boundary_info = igfb ( n_elems, 'boundary_info' )

    init_wtube%windex1 = n_elems/2+1-world_tube_width
    init_wtube%windex2 = n_elems/2+world_tube_width
    if ( init_wtube%wsize > 0 ) then
      ! If we are using a real world tube...
      do i = 1, n_elems
        if ( i>=init_wtube%windex1 .and. i<=init_wtube%windex2) then
          init_wtube%win%elems(i)%var(:) = 1.0_wp
          init_wtube%dwin%elems(i)%var(:) = 0.0_wp
          init_wtube%d2win%elems(i)%var(:) = 0.0_wp
        else
          init_wtube%win%elems(i)%var(:) = 0.0_wp
          init_wtube%dwin%elems(i)%var(:) = 0.0_wp
          init_wtube%d2win%elems(i)%var(:) = 0.0_wp
        end if
        if ( i==init_wtube%windex1-1 ) then
          init_wtube%boundary_info%elems(i)%bvar(1) = izero
          init_wtube%boundary_info%elems(i)%bvar(2) = +1
        else if ( i==init_wtube%windex1 ) then
          init_wtube%boundary_info%elems(i)%bvar(1) = -1
          init_wtube%boundary_info%elems(i)%bvar(2) = izero
        else if ( i==init_wtube%windex2 ) then
          init_wtube%boundary_info%elems(i)%bvar(1) = izero
          init_wtube%boundary_info%elems(i)%bvar(2) = -1
        else if ( i==init_wtube%windex2+1) then
          init_wtube%boundary_info%elems(i)%bvar(1) = +1
          init_wtube%boundary_info%elems(i)%bvar(2) = izero
        else
          init_wtube%boundary_info%elems(i)%bvar(:) = izero
        end if
      end do
    else
      ! If the worldtube has zero width the boundary info contains information
      ! about whether we need the left (-1) or right (+1) derivative of the
      ! singular field.
      do i = 1, n_elems
        init_wtube%win%elems(i)%var(:) = 0.0_wp
        init_wtube%dwin%elems(i)%var(:) = 0.0_wp
        init_wtube%d2win%elems(i)%var(:) = 0.0_wp
      end do
      init_wtube%boundary_info%elems(n_elems/2)%bvar(1) = izero
      init_wtube%boundary_info%elems(n_elems/2)%bvar(2) = -1
      init_wtube%boundary_info%elems(n_elems/2+1)%bvar(1) = 1
      init_wtube%boundary_info%elems(n_elems/2+1)%bvar(2) = izero
    end if

  end procedure init_wtube


  module procedure is_boundary

    implicit none

    integer(ip), dimension(-1:1) :: bindex = (/ 1, 0, 2 /)

    select case ( dir )
    case (-1, 1)
      if ( this%boundary_info%elems(n)%bvar(bindex(dir)) /= izero ) then
        is_boundary = .true.
      else
        is_boundary = .false.
      end if
    case default
      print*,'world_tube%is_boundary called with invalid dir'
      stop
    end select

  end procedure is_boundary


end submodule world_tube_implementation
