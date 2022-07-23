submodule(rk4_integrator) rk4_implementation
!! The implementation of the interfaces defined in [[rk4_integrator]].

contains

  module procedure rk4_ntemp

    implicit none

    ntemp = this%ntmp

  end procedure rk4_ntemp


  module procedure rk4_init

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
      print*,'Setting up RK4 pointer to ', eqs(i)%p%ename
      this%eqs(i)%p => eqs(i)%p
    end do
    print*,'RK4 initialized with ', this%nequations, ' equations'

    do i = 1, this%nequations
      call this%eqs(i)%p%save_globals_1 ( )
    end do
        
    do i = 1, this%nequations
      call this%eqs(i)%p%save_globals_2 ( )
    end do

    allocate ( character(3) :: this%iname )
    this%iname = 'rk4' 
  end procedure rk4_init


  module procedure rk4_shutdown

    implicit none

    integer(ip) :: i

    do i = 1, this%nequations
      nullify ( this%eqs(i)%p )
    end do

    deallocate ( this%eqs, this%iname )

    print*,'RK4 integrator shutdown'

  end procedure rk4_shutdown


  module procedure rk4_step

    use time_info, only : save_time, get_current_dtime, increment_time, &
                          restore_and_increment_time

    implicit none

    integer(ip), parameter :: nsteps = 5
    real(wp), dimension(nsteps), parameter :: &
      rk4a = (/ 0.0_wp, &
                -567301805773.0_wp/1357537059087.0_wp, &
                -2404267990393.0_wp/2016746695238.0_wp, &
                -3550918686646.0_wp/2091501179385.0_wp, &
                -1275806237668.0_wp/842570457699.0_wp /)
    real(wp), dimension(nsteps), parameter :: &
      rk4b = (/ 1432997174477.0_wp/9575080441755.0_wp, &
                5161836677717.0_wp/13612068292357.0_wp, &
                1720146321549.0_wp/2090206949498.0_wp, &
                3134564353537.0_wp/4481467310338.0_wp, &
                2277821191437.0_wp/14882151754819.0_wp /)
    real(wp), dimension(nsteps), parameter :: &
      rk4c = (/ 1432997174477.0_wp/9575080441755.0_wp, &
                2526269341429.0_wp/6820363962896.0_wp, &
                2006345519317.0_wp/3224310063776.0_wp, &
                2802321613138.0_wp/2924317926251.0_wp, &
                1.0_wp /)
    integer(ip) :: i, j, k, l, nequations
    real(wp) :: localtime
    real(wp) :: dt, dtime

    nequations = this%nequations
    do j = 1, nequations
      call this%eqs(j)%p%set_to_zero(1)
    end do

    call save_time ()
    dtime = get_current_dtime ( )
    do i = 1, nsteps
!      print*,'step ', i

      do j = 1, nequations
        call this%eqs(j)%p%load_globals ( )
      end do
 
      do j = 1, nequations
!        print*,'j = ', j
        call this%eqs(j)%p%rhs()
        call this%eqs(j)%p%update_vars ( 1, 1, source2=-1, &
                                         scalar=rk4a(i), scalar2=dtime )
        call this%eqs(j)%p%update_vars ( 0, 0, source2=1, scalar2=rk4b(i) )
      end do
      dt = rk4c(i)*dtime
      call restore_and_increment_time ( dt )

      if (i<nsteps) then
        this%is_complete = .false.
      else
        this%is_complete = .true.
      end if

      do j = 1, nequations
        call this%eqs(j)%p%save_globals_1 ( )
      end do

      do j = 1, nequations
        call this%eqs(j)%p%save_globals_2 ( )
      end do
    end do
    
!    call restore_and_increment_time ( dtime )

  end procedure rk4_step


  module procedure rk4_complete_step

    implicit none

    is_complete = this%is_complete
    
  end procedure rk4_complete_step

end submodule rk4_implementation
