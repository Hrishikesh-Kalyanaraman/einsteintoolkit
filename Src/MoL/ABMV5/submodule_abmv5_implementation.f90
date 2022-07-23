submodule(abmv5_integrator) abmv5_implementation
!! The implementation of the interfaces defined in [[abmv5_integrator]].

contains

  module procedure abmv5_ntemp

    implicit none

    ntemp = this%ntmp

  end procedure abmv5_ntemp


  module procedure abmv5_init

    use rk4_integrator
    use time_info, only : get_current_dtime, get_current_time, &
                          set_dtime, init_time
    use numerics, only : factorial

    implicit none

    real(wp), dimension(5,4) :: deriv_coeff = reshape ( (/ &
      -25.0_wp/12.0_wp, 4.0_wp, -3.0_wp, 4.0_wp/3.0_wp, -1.0_wp/4.0_wp, &
      35.0_wp/12.0_wp, -26.0_wp/3.0_wp, 19.0_wp/2.0_wp, -14.0_wp/3.0_wp, &
      11.0_wp/12.0_wp, &
      -5.0_wp/2.0_wp, 9.0_wp, -12.0_wp, 7.0_wp, -3.0_wp/2.0_wp, &
      1.0_wp, -4.0_wp, 6.0_wp, -4.0_wp, 1.0_wp /), &
      (/ 5,4 /) )
    type(rk4) :: start_rk
    real(wp) :: local_dt, dt_save, local_t, t_save, fac
    integer(ip) :: nequations, i, j, k, l

    this%is_complete = .true.
    this%nequations = size(eqs)
    nequations = this%nequations

    allocate ( this%eqs(this%nequations) )
    do j = 1, this%nequations
      print*, 'Setting up ABMV5 pointer to ', eqs(j)%p%ename
      this%eqs(j)%p => eqs(j)%p
    end do
    print*, 'ABMV5 initialized with ', this%nequations, ' equations.'

! Use rk4 to start the Adams-Bashford-Moulton scheme.
    call start_rk%init ( eqs )

    dt_save = get_current_dtime ( )
    t_save = get_current_time ( )

! As rk4 is less accurate than 5th order ABMV, we use a smaller delta t.
    local_dt = 0.5_wp * dt_save
    call set_dtime ( local_dt )

    do j = 1, nequations
! Calculate the RHS at the initial time.
      call this%eqs(j)%p%rhs ( )
! Store the state in the last temporary variable. RK4 will use the first
! temporary variable, so we have to copy this back once the RK4 steps are done.
      call this%eqs(j)%p%update_vars ( 0, this%ntmp )
! Store the scaled RHS in the second temporary variable.
      call this%eqs(j)%p%update_vars ( -1, 2, scalar=dt_save )
    end do
 
    do i = 3, this%order+1
! Take an RK4 step.
!      call set_dtime ( local_dt )
      call start_rk%step ( )
! For all the equations...
      do j = 1, nequations
! Calculate the RHS at the current time.
        call this%eqs(j)%p%rhs ( )
! Store the scaled RHS in the appropriate temporary variable.
        call this%eqs(j)%p%update_vars ( -1, i, scalar=dt_save )
      end do
! As delta t might have changed, reset it to local_dt. We don't need to adjust
! delta t, since we are only evolving for a short time without having to worry
! about stability.
      call set_dtime ( local_dt )
    end do

    call start_rk%shutdown ( )

! Now restore the state from the last temporary storage.
    do j = 1, nequations
      call this%eqs(j)%p%update_vars ( this%ntmp, 1 )
    end do

! Calculate time derivatives and store them in the final slots of the temporary
! storage.
    do k = 1, this%order-1
! First store the first coefficient times the RHS at the initial time.
      do j = 1, nequations
        call this%eqs(j)%p%update_vars ( 2, k+6, scalar=deriv_coeff(1,k) )
      end do
! Then add the derivative coefficients times the RHS at the following times.
      do l = 2, this%order
        do j = 1, nequations
          call this%eqs(j)%p%update_vars ( k+6, k+6, source2=l+1, &
                                           scalar2=deriv_coeff(l,k) )
        end do
      end do
    end do

! Scale the higher time derivatives and store them in the temporary variable.
    do k = 3, this%order+1
      do j = 1, nequations
! The time derivatives are calculated based on time steps of local_dt, while
! we need to store derivatives scaled with the correct power of dt_save.
        fac = (dt_save/local_dt)**(k-2)/factorial(k-1)
        call this%eqs(j)%p%update_vars ( k+4, k, scalar=fac )
      end do
    end do

! Restore the initial state.
    do j = 1, nequations
      call this%eqs(j)%p%update_vars ( 1, 0 )
    end do

! Restore t and delta t.
    call init_time ( t_save )
    call set_dtime ( dt_save )
    this%last_dt = dt_save

    allocate ( character(5) :: this%iname )
    this%iname = 'abmv5'

  end procedure abmv5_init


  module procedure abmv5_shutdown

    implicit none

    integer(ip) :: i

    do i = 1, this%nequations
      nullify ( this%eqs(i)%p )
    end do

    deallocate ( this%eqs, this%iname )

    print*,'ABMV5 integrator shutdown'

  end procedure ABMV5_shutdown


  module procedure abmv5_step

    use time_info, only : get_current_dtime, increment_time

    implicit none

    integer(ip), parameter :: order = 5
    integer(ip), parameter :: ns = order+1
    real(wp), dimension(ns), parameter :: &
      rcoeff = (/ 251.0_wp/720.0_wp, 1.0_wp, 25.0_wp/24.0_wp, &
                  35.0_wp/72.0_wp, 5.0_wp/48.0_wp, 1.0_wp/120.0_wp /)
    real(wp), dimension(3:ns), parameter :: &
      fac1 = (/ 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp /)
    real(wp), dimension(4:ns), parameter :: &
      fac2 = (/ 3.0_wp, 6.0_wp, 10.0_wp /)
    real(wp), dimension(5:ns), parameter :: fac3 = (/ 4.0_wp, 10.0_wp /)
    real(wp), dimension(ns:ns), parameter :: fac4 = (/ 5.0_wp  /)
    real(wp), dimension(ns) :: dtscale
    real(wp) :: dtime, dtfac
    integer(ip) :: nequations, i, j

    nequations = this%nequations

    dtime = get_current_dtime ( )
! Calculate the ratio between the current timestep and the last timestep
    dtfac = dtime/this%last_dt

! Calculate scale factors for rescaling of time derivatives to match the new
! timestep.
    do i = 1, ns
      dtscale(i) = dtfac**(i-1)
    end do

! Store the current timestep as the last timestep for the next step.
    this%last_dt = dtime

! Rescale the time derivatives.
    if ( dtfac/=1.0_wp ) then
      do j = 1, nequations
        do i = 2, ns
          call this%eqs(j)%p%update_vars ( i, i, scalar=dtscale(i) )
        end do
      end do
    end if

! The Adams-Bashford predictor step...
    do j = 1, nequations
! Copy the statevector in the temporary variable to the data variable.
! var_data = tmp_data(1)
      call this%eqs(j)%p%update_vars ( 1, 0 )
      do i = 2, ns
! var_data = var_data + tmp_data(i)
        call this%eqs(j)%p%update_vars ( 0, 0, source2=i )
! The higher derivatives does not depend on the RHS calculated for the
! corrector step. Calculating these quantities here means that we can
! reuse scalartmp(j,2) and tempk (:,:,j,2,:) to store alpha.
        if ( i>2 ) then
! tmp_data(2) = tmp_data(2)+fac1(i)*tmp_data(i)
          call this%eqs(j)%p%update_vars ( 2, 2, source2=i, scalar2=fac1(i) )
        end if
        if ( i>3 ) then
! tmp_data(3) = tmp_data(3)+fac2(i)*tmp_data(i)
          call this%eqs(j)%p%update_vars ( 3, 3, source2=i, scalar2=fac2(i) )
        end if
        if ( i>4 ) then
! tmp_data(4) = tmp_data(4)+fac3(i)*tmp_data(i)
          call this%eqs(j)%p%update_vars ( 4, 4, source2=i, scalar2=fac3(i) )
        end if
        if ( i>5 ) then
! tmp_data(5) = tmp_data(5)+fac4(i)*tmp_data(i)
          call this%eqs(j)%p%update_vars ( 5, 5, source2=i, scalar2=fac4(i) )
        end if
      end do
! Store the predicted state in the first temporary storage variable.
! tmp_data(1) = var_data
      call this%eqs(j)%p%update_vars ( 0, 1 )
    end do

! Increment time before the RHS call
    call increment_time ( dtime )

! Call save globals after the new intermediate state has been computed.
    do j = 1, nequations
      call this%eqs(j)%p%save_globals_1 ( )
    end do

    do j = 1, nequations
      call this%eqs(j)%p%save_globals_2 ( )
    end do

    do j = 1, nequations
! Call load globals before the RHS is computed for the corrector step.
      call this%eqs(j)%p%load_globals ( )
    end do

    do j = 1, nequations
! Call the RHS routine with the predicted state in preparation for the
! corrector state.
      call this%eqs(j)%p%rhs()
    end do

! The Adams-Moulton corrector step.

    do j = 1, nequations
! Calculate alpha and store it in the temporary variable ns+1.
! tmp_data(ns+1) = dtime*rhs_data - 1.0*tmp_data(2)
      call this%eqs(j)%p%update_vars ( -1, ns+1, source2=2, scalar=dtime, &
                                       scalar2=-1.0_wp)
    end do

    do j = 1, nequations
! Use alpha to correct the statevector.
! tmp_data(1) = tmp_data(1) + rcoeff(1)*tmp_data(ns+1)
      call this%eqs(j)%p%update_vars ( 1, 1, source2=ns+1, scalar2=rcoeff(1) )
      do i = 3, ns
! tmp_data(i) = tmp_data(i) + rcoeff(i)*tmp_data(ns+1)
        call this%eqs(j)%p%update_vars ( i, i, source2=ns+1, scalar2=rcoeff(i) )
      end do
! Copy the statevector in the temporary variable to the data variable.
! var_data = tmp_data(1)
      call this%eqs(j)%p%update_vars ( 1, 0 )
    end do

! Call save globals after the new intermediate state has been computed.
    do j = 1, nequations
      call this%eqs(j)%p%save_globals_1 ( )
    end do

    do j = 1, nequations
      call this%eqs(j)%p%save_globals_2 ( )
    end do

    do j = 1, nequations
! Call load globals before the RHS is computed for the corrector step.
      call this%eqs(j)%p%load_globals ( )
    end do

    do j = 1, nequations
! Call the RHS routine with the corrected state in preparation for the
! next predictor step. This RHS evaluation is performed at the same time
! as the last one, so no need to increment time.
      call this%eqs(j)%p%rhs()
    end do

    do j = 1, nequations
! Store the rescaled rhs in temporary variable 2.
! tmp_data(2) = dtime*rhs_data
      call this%eqs(j)%p%update_vars ( -1, 2, scalar=dtime )
    end do

  end procedure abmv5_step


  module procedure abmv5_complete_step

  end procedure abmv5_complete_step

end submodule abmv5_implementation
