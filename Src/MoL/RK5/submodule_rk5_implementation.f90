! Need to check this for consistency with regards to save_globals and time
! updates.
submodule(rk5_integrator) rk5_implementation
!! The implementation of the interfaces defined in [[rk5_integrator]].

contains

  module procedure rk5_ntemp

    implicit none

    ntemp = this%ntmp

  end procedure rk5_ntemp


  module procedure rk5_init

    implicit none

    integer(ip) :: i

    this%is_complete = .true.
    this%nequations = size(eqs)
! The following line currently needs to be replaced by the next 6 in order
! to avoid a segmentation fault. I suspect this is a problem with the gnu
! compiler.
!    allocate ( this%eqs, source=eqs )
    allocate ( this%eqs(this%nequations) )
    do i = 1, this%nequations
      print*,'Setting up RK5 pointer to ', eqs(i)%p%ename
      this%eqs(i)%p => eqs(i)%p
    end do
    print*,'RK5 initialized with ', this%nequations, ' equations'

    do i = 1, this%nequations
      call this%eqs(i)%p%save_globals_1 ( )
    end do
    do i = 1, this%nequations
      call this%eqs(i)%p%save_globals_2 ( )
    end do
        
    allocate ( character(3) :: this%iname )
    this%iname = 'rk5'

  end procedure rk5_init


  module procedure rk5_shutdown

    implicit none

    integer(ip) :: i

    do i = 1, this%nequations
      nullify ( this%eqs(i)%p )
    end do

    deallocate ( this%eqs, this%iname )

    print*,'RK5 integrator shutdown'

  end procedure rk5_shutdown


  module procedure rk5_step

    use time_info, only : save_time, get_current_dtime, increment_time, &
                          restore_and_increment_time

    implicit none

    integer(ip), parameter :: nsteps = 7
    real(wp), dimension(nsteps-1,nsteps-1), parameter :: &
      rk5a = reshape ( (/ 0.2_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, &                       ! First row
                          3.0_wp/32.0_wp, 5.0_wp/32.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, &               ! Second row
                          1.0_wp/8.0_wp, -5.0_wp/8.0_wp, 1.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, &               ! Third row
                          -1.0_wp/8.0_wp, 25.0_wp/8.0_wp, -3.0_wp, &
                          0.5_wp, 0.0_wp, 0.0_wp, &               ! Fourth row
                          207.0_wp/1372.0_wp, -405.0_wp/686.0_wp, &
                          297.0_wp/343.0_wp, 1485.0_wp/9604.0_wp, &
                          297.0_wp/4802.0_wp, 0.0_wp, &           ! Fifth row
                          49.0_wp/36096.0_wp, 1855.0_wp/12032.0_wp, &
                          7945.0_wp/16544.0_wp, -12845.0_wp/24064.0_wp, &
                          -315.0_wp/24064.0_wp, 156065.0_wp/198528.0_wp &
                           /), &                                  ! Sixth row
                            (/ nsteps-1, nsteps-1 /) )
    real(wp), dimension(nsteps), parameter :: &
      rk5b =  (/ 83.0_wp/945.0_wp, 0.0_wp, 248.0_wp/825.0_wp, &
                 41.0_wp/180.0_wp, 1.0_wp/36.0_wp, 2401.0_wp/38610_wp, &
                 6016.0_wp/20475.0_wp /)
    real(wp), dimension(nsteps), parameter :: &
      rk5c = (/ 0.2_wp, 0.25_wp, 0.5_wp, 0.5_wp, 9.0_wp/14.0_wp, &
                7.0_wp/8.0_wp, 1.0_wp /)
    integer(ip) :: i, j, k, l, nequations
    real(wp) :: localtime
    real(wp) :: dt, dtime

    nequations = this%nequations

    do j = 1, nequations
! Load the globals before calculating the first RHS.
      call this%eqs(j)%p%load_globals ( )
    end do

    do j = 1, nequations
! Set the temporary storage space to zero.
      do k = 1, this%ntmp
        call this%eqs(j)%p%set_to_zero ( k )
      end do
! Copy the state vector to the last temporage storage space
      call this%eqs(j)%p%update_vars ( 0, nsteps )

! Calculate the RHS for the first sub-step
      call this%eqs(j)%p%rhs ( )
    end do

    call save_time ()
    dtime = get_current_dtime ( )
! Loop over the first nsteps-1 sub-steps.
    do i = 1, nsteps-1
      dt = rk5c(i)*dtime
      call restore_and_increment_time ( dt )

      do j = 1, nequations
! Store the RHS in the i'th temporary storage space.
        call this%eqs(j)%p%update_vars ( -1, i )
! Add the first RHS to the stored initial state and store it in the variable
! space.
        call this%eqs(j)%p%update_vars ( nsteps, 0, source2=1, &
                                         scalar2=dtime*rk5a(1,i) )
        do k = 2, i
! Add the rest of the stored RHS to the new intermediate state.
          call this%eqs(j)%p%update_vars ( 0, 0, source2=k, &
                                           scalar2=dtime*rk5a(k,i) )
        end do
      end do

! Call save globals after the new intermediate state has been computed.
      do j = 1, nequations
        call this%eqs(j)%p%save_globals_1 ( )
      end do

      do j = 1, nequations
        call this%eqs(j)%p%save_globals_2 ( )
      end do
 
      do j = 1, nequations
! Call load globals before the next RHS is computed.
        call this%eqs(j)%p%load_globals ( )
      end do
 
      do j = 1, nequations
! Call the RHS routine with this intermediate state.
        call this%eqs(j)%p%rhs()
      end do
    end do

! Take the last sub-step.
    dt = rk5c(nsteps)*dtime
    call restore_and_increment_time ( dt )

    do j = 1, nequations
! Add the first RHS to the stored initial state and store it in the variable
! space.
      call this%eqs(j)%p%update_vars ( nsteps, 0, source2=1, &
                                       scalar2=dtime*rk5b(1) )

      do k = 2, nsteps-1
! Add the remaining RHS that was stored in temporary space to the state vector.
        call this%eqs(j)%p%update_vars ( 0, 0, source2=k, &
                                         scalar2=dtime*rk5b(k) )
      end do

! Add the last RHS (not stored in temporary space) to form the final state
! vector at the full timestep.
      call this%eqs(j)%p%update_vars ( 0, 0, source2=1, &
                                       scalar2=dtime*rk5b(nsteps) )
    end do

! Call save globals after the full time step has been computed.
    do j = 1, nequations
      call this%eqs(j)%p%save_globals_1 ( )
    end do

    do j = 1, nequations
      call this%eqs(j)%p%save_globals_2 ( )
    end do
 
    call restore_and_increment_time ( dtime )

  end procedure rk5_step


  module procedure rk5_complete_step

  end procedure rk5_complete_step
end submodule rk5_implementation
