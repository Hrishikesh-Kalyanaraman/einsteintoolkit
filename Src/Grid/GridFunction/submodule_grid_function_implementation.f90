submodule(grid_function) grid_function_implementation
!! The implementation of the interfaces defined in [[grid_function]].

contains

  module procedure init_rgf

    implicit none

    integer(ip) :: i

    integer :: allocation_status

    init_rgf%n = n
    init_rgf%io_id = -1
    if ( .not. allocated(init_rgf%vname) ) then
      allocate ( character(len(trim(var_name))) :: init_rgf%vname, &
                   stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating grid function ', var_name
        stop
      end if
    end if
    init_rgf%vname = var_name
    if ( .not. allocated(init_rgf%elems) ) then
      allocate ( init_rgf%elems(n), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Allocation error for real element data'
        stop
      end if 
    end if
    do i = 1, n
      init_rgf%elems(i) = element_rdata(order)
    end do
  end procedure init_rgf


  module procedure output_rgf

    use output_base
    use time_info

    implicit none

    character(len=len_trim(this%vname)+4) :: filename
    integer(ip) :: i, j, io_id

    io_id = this%io_id
    if (io_id < 0 ) then
      filename = trim(this%vname) // '.asc'
      io_id = next_available_io_id ()
      this%io_id = io_id
      print*,'Opening ', filename, ' with id ', io_id
      open(io_id, file=filename, status='replace', action='write')
    end if
    write(io_id,*)
    write(io_id,*)
    write(io_id,*) '#time = ', get_current_time ( )

    do i = 1, this%n
      do j = 1, this%elems(i)%order+1
        write(io_id,'(*(es23.15e3,1x))') coord%elems(i)%var(j), &
                                        this%elems(i)%var(j)
      end do
    end do

  end procedure output_rgf


  module procedure deallocate_rgf

    implicit none

    integer :: allocation_status

    if ( allocated(this%elems) ) then
      deallocate ( this%elems, this%vname, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for real element data'
        stop
      end if 
    end if

    if (this%io_id >=0) then
      close(this%io_id)
    end if  

  end procedure deallocate_rgf


  module procedure init_cgf

    implicit none

    integer(ip) :: i

    integer :: allocation_status

!    print*,'Setting number of elements'
    init_cgf%n = n
    init_cgf%io_id = -1
!    print*,'allocating array of elements'
    if ( .not. allocated(init_cgf%vname) ) then
      allocate ( character(len(trim(var_name))) :: init_cgf%vname, &
                   stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating grid function ', var_name
        stop
      end if
    end if
    init_cgf%vname = var_name
    if ( .not. allocated(init_cgf%elems) ) then
      allocate ( init_cgf%elems(n), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Allocation error for real element data'
        stop
      end if 
    end if
    do i = 1, n
      init_cgf%elems(i) = element_cdata(order)
    end do
!    print*,'Done initializing elements'
  end procedure init_cgf


  module procedure output_cgf

    use output_base
    use time_info

    implicit none

    character(len=len_trim(this%vname)+4) :: filename
    integer(ip) :: i, j, io_id

    io_id = this%io_id
    if (io_id<0) then 
      filename = trim(this%vname) // '.asc'
      io_id = next_available_io_id()
      this%io_id = io_id
      print*,'Opening ', filename, ' with id ', io_id
      open(io_id, file=filename, status='replace', action='write')
    end if
    write(io_id,*)
    write(io_id,*)
    write(io_id,*) '#time = ', get_current_time ( )

    do i = 1, this%n
      do j = 1, this%elems(i)%order+1
        write(io_id,'(*(es23.15e3,1x))') coord%elems(i)%var(j), &
                                     real(this%elems(i)%var(j),wp), &
                                     aimag(this%elems(i)%var(j))
      end do
    end do

  end procedure output_cgf


  module procedure deallocate_cgf

    implicit none

    integer :: allocation_status

    if ( allocated(this%elems) ) then
      deallocate ( this%elems, this%vname, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for real element data'
        stop
      end if 
    end if

    if (this%io_id>=0) then
      close(this%io_id)
    end if

  end procedure deallocate_cgf


  module procedure init_igfb

    implicit none

    integer(ip) :: i

    integer :: allocation_status

    init_igfb%n = n
    init_igfb%iob_id = -1
    if ( .not. allocated(init_igfb%vname) ) then
      allocate ( character(len(trim(var_name))) :: init_igfb%vname, &
                   stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating boundary grid function ', var_name
        stop
      end if
    end if
    init_igfb%vname = var_name
    if ( .not. allocated(init_igfb%elems) ) then
      allocate ( init_igfb%elems(n), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Allocation error for integer element boundary data'
        stop
      end if 
    end if
    do i = 1, n
      init_igfb%elems(i) = element_boundary_idata( (/ izero, izero /) )
    end do
  end procedure init_igfb


  module procedure output_igfb

    use output_base
    use time_info

    implicit none

    character(len=len_trim(this%vname)+4) :: filename
    integer(ip) :: i, j, iob_id
    integer(ip), dimension(2) :: ind

     iob_id = this%iob_id
     if (iob_id<0) then
       filename = trim(this%vname) // '.asc'
       iob_id = next_available_io_id ()
       this%iob_id = iob_id
       print*,'Opening ', filename, ' with id ', iob_id
       open(iob_id, file=filename, status='replace', action='write')
     end if
     write(iob_id,*)
     write(iob_id,*)
     write(iob_id,*) '#time = ', get_current_time ( )

     ind(1) = 1
     do i = 1, this%n
       ind(2) = coord%elems(i)%order + 1
       do j = 1, 2
         write(iob_id,'(es23.15e3,1x,i0,1x)') coord%elems(i)%var(ind(j)), &
                                          this%elems(i)%bvar(j)
       end do
     end do

  end procedure output_igfb


  module procedure deallocate_igfb

    implicit none

    integer :: allocation_status

    if ( allocated(this%elems) ) then
      deallocate ( this%elems, this%vname, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for real element data'
        stop
      end if 
    end if

    if (this%iob_id>=0) then
      close(this%iob_id)
    end if
  end procedure deallocate_igfb


  module procedure init_rgfb

    implicit none

    integer(ip) :: i

    integer :: allocation_status

    init_rgfb%n = n
    init_rgfb%iob_id = -1
    if ( .not. allocated(init_rgfb%vname) ) then
      allocate ( character(len(trim(var_name))) :: init_rgfb%vname, &
                   stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating boundary grid function ', var_name
        stop
      end if
    end if
    init_rgfb%vname = var_name
    if ( .not. allocated(init_rgfb%elems) ) then
      allocate ( init_rgfb%elems(n), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Allocation error for real element data'
        stop
      end if 
    end if
    do i = 1, n
      init_rgfb%elems(i) = element_boundary_rdata( (/ rzero, rzero /) )
    end do
  end procedure init_rgfb


  module procedure output_rgfb

    use output_base
    use time_info

    implicit none

    character(len=len_trim(this%vname)+4) :: filename
    integer(ip) :: i, j, iob_id
    integer(ip), dimension(2) :: ind

     iob_id = this%iob_id
     if (iob_id<0) then
       filename = trim(this%vname) // '.asc'
       iob_id = next_available_io_id ()
       this%iob_id = iob_id
       print*,'Opening ', filename, ' with id ', iob_id
       open(iob_id, file=filename, status='replace', action='write')
     end if
     write(iob_id,*)
     write(iob_id,*)
     write(iob_id,*) '#time = ', get_current_time ( )
     print*,'Writing to file with id: ', iob_id

     ind(1) = 1
     do i = 1, this%n
       ind(2) = coord%elems(i)%order + 1
       do j = 1, 2
         write(iob_id,'(*(es23.15e3,1x))') coord%elems(i)%var(ind(j)), &
                                          this%elems(i)%bvar(j)
       end do
     end do

  end procedure output_rgfb


  module procedure deallocate_rgfb

    implicit none

    integer :: allocation_status

    if ( allocated(this%elems) ) then
      deallocate ( this%elems, this%vname, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for real element data'
        stop
      end if 
    end if

    if (this%iob_id>=0) then
      close(this%iob_id)
    end if
  end procedure deallocate_rgfb


  module procedure init_cgfb

    implicit none

    integer(ip) :: i

    integer :: allocation_status

    init_cgfb%n = n
    init_cgfb%iob_id = -1
    if ( .not. allocated(init_cgfb%vname) ) then
      allocate ( character(len(trim(var_name))) :: init_cgfb%vname, &
                   stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating boundary grid function ', var_name
        stop
      end if
    end if
    init_cgfb%vname = var_name
    if ( .not. allocated(init_cgfb%elems) ) then
      allocate ( init_cgfb%elems(n), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Allocation error for real element data'
        stop
      end if 
    end if
    do i = 1, n
      init_cgfb%elems(i) = element_boundary_cdata ( (/ czero, czero /) )
    end do
  end procedure init_cgfb


  module procedure output_cgfb

    use output_base
    use time_info

    implicit none

    character(len=len_trim(this%vname)+4) :: filename
    integer(ip) :: i, j, iob_id
    integer(ip), dimension(2) :: ind

    iob_id = this%iob_id
    if (iob_id <0) then 
      filename = trim(this%vname) // '.asc'
      iob_id = next_available_io_id ()
      this%iob_id = iob_id
      print*,'Opening ', filename, ' with id ', iob_id
      open(iob_id, file=filename, status='replace', action='write')
    end if
    write(iob_id,*)
    write(iob_id,*)
    write(iob_id,*) '#time = ', get_current_time ( )

    ind(1) = 1
    do i = 1, this%n
      ind(2) = coord%elems(i)%order + 1
      do j = 1, 2
        write(20,'(*(es23.15e3,1x))') coord%elems(i)%var(ind(j)), &
                                     real(this%elems(i)%bvar(j),wp), &
                                     aimag(this%elems(i)%bvar(j))
      end do
    end do

    close(iob_id)

  end procedure output_cgfb


  module procedure mult_sc_cgf

    implicit none

    integer(ip) :: i, j, n
    
    n = this%n
    do i = 1, n
      this%elems(i)%var  = scalar*this%elems(i)%var
    end do
  end procedure mult_sc_cgf


  module procedure sc_mult_gf_cgf

    implicit none

    integer(ip) :: i, n1, n2, order1, order2

    n1 = this%n
    n2 = gf%n
    if ( n1 /= n2 ) then
      print*,'Error (sc_mult_gf_cgf): Attempting to scale a grid function and store it in a grid function with different number of elements'
      stop
    end if
    do i = 1, n1
      this%elems(i)%var = scalar*gf%elems(i)%var
    end do
  end procedure sc_mult_gf_cgf


  module procedure copy_cgf

    implicit none

    integer(ip) :: i, j, n1, n2, order1, order2

    n1 = this%n
    n2 = gf%n
    if ( n1 /= n2 ) then
      print*,'Error (copy_cgf): Attempting to copy between 2 grid functions with different number of elements'
      stop
    end if
    do i = 1, n1
      order1 = this%elems(i)%order
      order2 = gf%elems(i)%order
      if ( order1 /= order2 ) then
        print*,'Error (add_gf_cgf): Attempting to copy between 2 grid functions with different order elements'
        stop
      end if
      this%elems(i)%var = gf%elems(i)%var
    end do
  end procedure copy_cgf


  module procedure zero_cgf

    implicit none

    integer(ip) :: i

    do i = 1, this%n
      this%elems(i)%var = czero
    end do

  end procedure zero_cgf


  module procedure add_gf_cgf

    implicit none

    integer(ip) :: i, j, n1, n2, order1, order2

    n1 = this%n
    n2 = gf%n
    if ( n1 /= n2 ) then
      print*,'Error (add_gf_cgf): Attempting to add 2 grid functions with different number of elements'
      stop
    end if
    do i = 1, n1
      order1 = this%elems(i)%order
      order2 = gf%elems(i)%order
      if ( order1 /= order2 ) then
        print*,'Error (add_gf_cgf): Attempting to add 2 grid functions with different order elements'
        stop
      end if
      this%elems(i)%var = this%elems(i)%var + gf%elems(i)%var
    end do
  end procedure add_gf_cgf


  module procedure add_sc_mult_gf_cgf

    implicit none

    integer(ip) :: i, j, n1, n2, order1, order2

    n1 = this%n
    n2 = gf%n
    if ( n1 /= n2 ) then
      print*,'Error (add_sc_mult_gf): Attempting to add grid functions with different number of elements'
      stop
    end if
    do i = 1, n1
      order1 = this%elems(i)%order
      order2 = gf%elems(i)%order
      if ( order1 /= order2 ) then
        print*,'Error (add_sc_mult_gf): Attempting to add 2 grid functions with different order elements'
        stop
      end if
      this%elems(i)%var = this%elems(i)%var + scalar*gf%elems(i)%var
    end do

  end procedure add_sc_mult_gf_cgf


  module procedure mult_sc_add_sc_mult_gf_cgf

    implicit none

    integer(ip) :: i, j, n1, n2, order1, order2

    n1 = this%n
    n2 = gf%n
    if ( n1 /= n2 ) then
      print*,'Error (add_sc_mult_gf_cgf): Attempting to add grid functions with different number of elements'
      stop
    end if
    do i = 1, n1
      order1 = this%elems(i)%order
      order2 = gf%elems(i)%order
      if ( order1 /= order2 ) then
        print*,'Error (add_sc_mult_gf_cgf): Attempting to add 2 grid functions with different order elements'
        stop
      end if
      this%elems(i)%var = scalar1*this%elems(i)%var + scalar2*gf%elems(i)%var
    end do

  end procedure mult_sc_add_sc_mult_gf_cgf


  module procedure gf1_plus_sc_mult_gf2_cgf

    implicit none

    integer(ip) :: i, j, n1, n2, n3, order1, order2, order3

    n1 = this%n
    n2 = gf1%n
    n3 = gf2%n
    if ( n1 /= n2 .or. n1 /= n3 ) then
      print*,'Error (gf1_plus_sc_mult_gf2_cgf): Attempting to operate on grid functions with different number of elements'
      stop
    end if
    do i = 1, n1
      order1 = this%elems(i)%order
      order2 = gf1%elems(i)%order
      order3 = gf2%elems(i)%order
      if ( order1 /= order2 .or. order1 /= order3 ) then
        print*,'Error (gf1_plus_sc_mult_gf2_cgf): Attempting to operate on grid functions with different order elements'
        stop
      end if
      this%elems(i)%var = gf1%elems(i)%var + scalar*gf2%elems(i)%var
    end do
  end procedure gf1_plus_sc_mult_gf2_cgf


  module procedure sc_mult_gf1_plus_sc_mult_gf2_cgf

    implicit none

    integer(ip) :: i, j, n1, n2, n3, order1, order2, order3

    n1 = this%n
    n2 = gf1%n
    n3 = gf2%n
    if ( n1 /= n2 .or. n1 /= n3 ) then
      print*,'Error (sc_mult_gf1_plus_sc_mult_gf2_cgf): Attempting to operate on grid functions with different number of elements'
      stop
    end if
    do i = 1, n1
      order1 = this%elems(i)%order
      order2 = gf1%elems(i)%order
      order3 = gf2%elems(i)%order
      if ( order1 /= order2 .or. order1 /= order3 ) then
        print*,'Error (sc_mult_gf1_plus_sc_mult_gf2_cgf): Attempting to operate on grid functions with different order elements'
        stop
      end if
      this%elems(i)%var = scalar1*gf1%elems(i)%var + scalar2*gf2%elems(i)%var
    end do
  end procedure sc_mult_gf1_plus_sc_mult_gf2_cgf


  module procedure deallocate_cgfb

    implicit none

    integer :: allocation_status

    if ( allocated(this%elems) ) then
      deallocate ( this%elems, this%vname, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for real element data'
        stop
      end if 
    end if

    if (this%iob_id>=0) then
      close(this%iob_id)
    end if
  end procedure deallocate_cgfb


end submodule grid_function_implementation
