submodule(observers) observers_implementation
!! The implementation of the interfaces defined in [[observers]].

  implicit none

contains

  module procedure robs_init

    implicit none

    select type(object)
    type is (rgf)
      this%p => object
    class default
      print*,'Real grid function observer initialized with invalid type'
      stop
    end select
    allocate ( this%vname, source = this%p%vname )
    this%nradii = size(rad)
    this%radii = rad
    this%ioo_id = -1

    call find_indices ( rad, coord, this%elem_index, this%node_index )
    allocate ( this%extract_data(size(rad)) )

  end procedure robs_init
  

  module procedure robs_extract

    implicit none

    integer(ip) :: i, j, k

    do k = 1, this%nradii
      i = this%elem_index(k)
      j = this%node_index(k)
      this%extract_data(k) = this%p%elems(i)%var(j)
    end do
  end procedure robs_extract


  module procedure robs_output
  end procedure robs_output


  module procedure close_robserver
  end procedure close_robserver


  module procedure cobs_init

    implicit none

    select type(object)
    type is (cgf)
      this%p => object
    class default
      print*,'Compleex grid function observer initialized with invalid type'
      stop
    end select
    allocate ( this%vname, source = this%p%vname )
    this%nradii = size(rad)
    this%radii = rad
    this%ioo_id = -1

    print*,'rad = ', rad
    call find_indices ( rad, coord, this%elem_index, this%node_index )
    allocate ( this%extract_data(size(rad)) )

  end procedure cobs_init
  

  module procedure cobs_extract

    implicit none

    integer(ip) :: i, j, k

    do k = 1, this%nradii
      i = this%elem_index(k)
      j = this%node_index(k)
      this%extract_data(k) = this%p%elems(i)%var(j)
    end do
  end procedure cobs_extract


  module procedure cobs_output

    use output_base
    use time_info

    implicit none

    character(len=len_trim(this%vname)+12) :: filename
    character(:), allocatable :: key
    integer(ip) :: ioo_id

    ioo_id = this%ioo_id
    if (ioo_id<0) then
      filename = trim(this%vname) // '.extract.asc' 
      ioo_id = next_available_io_id ()
      this%ioo_id = ioo_id
      print*,'Opening ', filename, ' with id ', ioo_id
      open(ioo_id, file=filename, status='replace', action='write')
      key = generate_key ( this%radii )
      write(ioo_id,*) key
    end if
    
    write(ioo_id,'(*(es23.15e3,1x))') get_current_time ( ), this%extract_data 

  end procedure cobs_output


  module procedure close_cobserver

    deallocate (this%vname, this%radii, this%elem_index, this%node_index)
    if (this%ioo_id>=0) then
      close(this%ioo_id)
    end if
  end procedure close_cobserver

 
  module procedure find_indices

    implicit none

    integer(ip) :: i, j, k, n
    logical, dimension(size(rad)) :: found 
    real(wp) :: dist

    found = .false.
    n = size(rad)
 
    allocate ( elem_index(n), node_index(n) )

    do i = 1, coord%n
      do j = 1, coord%elems(i)%order+1
        do k = 1, n
          if (.not. found(k) ) then
            dist = abs( coord%elems(i)%var(j) - rad(k) )
            if ( rad(k)/=0.0_wp ) then
              dist = dist / abs( rad(k) )                            
            end if
            if ( dist < 1e-12 ) then
              found(k) = .true.
              elem_index(k) = i
              node_index(k) = j
            end if
          end if
        end do
      end do
    end do      
 
    if ( any (.not. found) ) then
      print*,'Error: No valid extraction point found for a cobserver'
      stop
    end if
  end procedure find_indices


  function generate_key ( rad ) result(key)
  !! Function that generates a head key or heading to add to the top of
  !! output files for an observer.

    real(wp), dimension(:), intent(in) :: rad
    !! The radii for this observer.
    character(:), allocatable :: key
    !! The return value is the key.
    character(9+size(rad)*30+1) :: tmp
    character(30) :: next
    integer(ip) :: i

    tmp = trim('#1: time') 
!    print*,'tmp = ', tmp
    do i = 1, size(rad)
      write(next,'(i1,a1,i1,a4,es23.15e3)') 2*i, ',', 2*i+1, ': r=', rad(i)
!      print*,'next = ', next
      tmp = trim(tmp)//' '//trim(next)
!      print*,'tmp = ', tmp
    end do
    key = tmp 

  end function generate_key

end submodule observers_implementation 
