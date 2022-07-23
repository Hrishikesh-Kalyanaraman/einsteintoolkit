submodule(pde_equations) pde_equations_implementation
!! The implementation of the interfaces defined in [[pde_equations]].

  implicit none

contains

  module procedure cpde_set_to_zero

    implicit none

    integer(ip) :: i, j, nmodes, nvars

    if (dest<-1 .or. dest>this%ntmp) then
      print*,'Error: set_to_zero called with invalid destination'
      stop
    end if

    nmodes = this%nmodes
    nvars = this%nvars
    do j = 1, nmodes
      do i = 1, nvars
        call this%data_pointer(i,j,dest)%p%zero()
      end do
    end do

  end procedure cpde_set_to_zero


  module procedure cpde_update_vars

    implicit none

! For source, source2 and dest:
! 0 => data
! -1 => rhs_data
! 1-... => tmp_data()

    integer(ip) :: i, j, k, nmodes, nvars, ntmp

    nmodes = this%nmodes
    nvars = this%nvars
    ntmp = this%ntmp

!    print*,'PDE equation update vars called'
    if (source<=-2 .or. source>ntmp) then
      print*,'Error: update_vars called with incorrect source argument'
      stop
    end if
    if (dest<=-1 .or. dest>ntmp) then
      print*,'Error: update_vars called with incorrect dest argument'
      stop
    end if
    if ( present(source2) .and. ( source2<=-2 .or. source2>ntmp) ) then
      print*,'Error: update_vars called with incorrect source2 argument'
      stop
    end if
    if ( present(source2) .and. (source2==dest) ) then
      print*,'Error: when source2 is present in update_vars it has to be different than dest'
      stop
    end if

! First take care of cases with only source and destination.
    if ( .not. present(source2) ) then
      if ( source == dest ) then
        if ( .not. present(scalar) ) then
          print*,'Error: if source and dest is the same in update_vars, scalar has to be present'
          stop
        else if ( scalar==1.0_wp ) then
          print*,'Error: calling update_vars with same source and destination and unit scalar is a noop and should not occur'
          stop
        else
!$OMP PARALLEL DO private(j,i) shared(nmodes,nvars,this)
          do j = 1, nmodes
            do i = 1, nvars
              call this%data_pointer(i,j,dest)%p%mult_sc( scalar )
            end do
          end do
!$OMP END PARALLEL DO
        end if
      else
        if ( .not. present(scalar) ) then
!$OMP PARALLEL DO private(j,i) shared(nmodes,nvars,this)
          do j = 1, nmodes
            do i = 1, nvars
              call this%data_pointer(i,j,dest)%p%copy( &
                               this%data_pointer(i,j,source)%p )
            end do
          end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO private(j,i) shared(nmodes,nvars,this)
          do j = 1, nmodes
            do i = 1, nvars
              call this%data_pointer(i,j,dest)%p%sc_mult_gf( scalar, &
                              this%data_pointer(i,j,source)%p )
            end do
          end do
!$OMP END PARALLEL DO
        end if
      end if
    else
      if (source==dest ) then
        if ( ( .not. present(scalar) ) .and. ( .not. present(scalar2) ) ) then
!$OMP PARALLEL DO private(j,i) shared(nmodes,nvars,this)
          do j = 1, nmodes
            do i = 1, nvars
              call this%data_pointer(i,j,dest)%p%add_gf( &
                              this%data_pointer(i,j,source2)%p )
            end do
          end do
!$OMP END PARALLEL DO
        else if ( ( present(scalar) ) .and. ( .not. present(scalar2) ) ) then
          print*,'Error: calling update_vars with scalar present and scalar2 not present is currently not supported'
          stop
        else if ( ( .not. present(scalar) ) .and. ( present(scalar2) ) ) then
!$OMP PARALLEL DO private(j,i) shared(nmodes,nvars,this)
          do j = 1, nmodes
            do i = 1, nvars
              call this%data_pointer(i,j,dest)%p%add_sc_mult_gf ( &
                              scalar2, this%data_pointer(i,j,source2)%p )
            end do
          end do
!$OMP END PARALLEL DO
        else if ( ( present(scalar) ) .and. ( present(scalar2) ) ) then
!$OMP PARALLEL DO private(j,i) shared(nmodes,nvars,this)
          do j = 1, nmodes
            do i = 1, nvars
              call this%data_pointer(i,j,dest)%p%mult_sc_add_sc_mult_gf ( &
                              scalar, scalar2, &
                              this%data_pointer(i,j,source2)%p )
            end do
          end do
!$OMP END PARALLEL DO
        end if
      else
        if ( present(scalar) .and. present(scalar2) ) then
!$OMP PARALLEL DO private(j,i) shared(nmodes,nvars,this)
          do j = 1, nmodes
            do i = 1, nvars
              call &
              this%data_pointer(i,j,dest)%p%sc_mult_gf1_plus_sc_mult_gf2 ( &
                              this%data_pointer(i,j,source)%p, scalar, &
                              this%data_pointer(i,j,source2)%p, scalar2 )
            end do
          end do
!$OMP END PARALLEL DO
        end if
        if ( ( present(scalar) ) .and. ( .not. present(scalar2) ) ) then
          print*,'Error: calling update_vars with source 2 present, source and dest different, scalar present and scalar 2 not present is currently not supported'
          stop
        end if
        if ( ( .not. present(scalar) ) .and. ( present(scalar2) ) ) then
!$OMP PARALLEL DO private(j,i) shared(nmodes,nvars,this)
          do j = 1, nmodes
            do i = 1, nvars
              call this%data_pointer(i,j,dest)%p%gf1_plus_sc_mult_gf2 ( &
                              this%data_pointer(i,j,source)%p, scalar2, &
                              this%data_pointer(i,j,source2)%p )
            end do
          end do
!$OMP END PARALLEL DO
        end if
        if ( ( .not. present(scalar) ) .and. ( .not. present(scalar2) ) ) then
          print*,'Error: calling update_vars with source 2 present, source and dest different and both scalar and scalar2 not present is currently not supported'
          stop
        end if
      end if
    end if

  end procedure cpde_update_vars


  module procedure cpde_output

    implicit none

  end procedure cpde_output


  module procedure cpde_print_data

    implicit none

    integer(ip) :: j

    print*,'pde data (1) = ', real(this%eq_data(1,1)%elems(16)%var,wp)
    print*,'pde data (2) = ', real(this%eq_data(2,1)%elems(16)%var,wp)
    print*,'pde data (3) = ', real(this%eq_data(3,1)%elems(16)%var,wp)
    print*,'rhs data (1) = ', real(this%eq_rhs_data(1,1)%elems(16)%var,wp)
    print*,'rhs data (2) = ', real(this%eq_rhs_data(2,1)%elems(16)%var,wp)
    print*,'rhs data (3) = ', real(this%eq_rhs_data(3,1)%elems(16)%var,wp)
    do j = 1, this%ntmp
      print*,'tmp data(', j, ',1) = ', real(this%eq_tmp_data(1,1,j)%elems(16)%var,wp)
      print*,'tmp data(', j, ',2) = ', real(this%eq_tmp_data(2,1,j)%elems(16)%var,wp)
      print*,'tmp data(', j, ',3) = ', real(this%eq_tmp_data(3,1,j)%elems(16)%var,wp)
    end do
  end procedure cpde_print_data

end submodule pde_equations_implementation
