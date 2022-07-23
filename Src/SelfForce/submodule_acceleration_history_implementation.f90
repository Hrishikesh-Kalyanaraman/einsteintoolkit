submodule(acceleration_history) acceleration_history_implementation
!! The implementation of the interfaces defined in [[acceleration_history]].

contains

  module procedure init_timelevels

    use parameters, only : use_hermite_extrapolation, hermite_order, &
                           use_smooth_derivs_for_hermite

    implicit none

    integer :: allocation_status

    if ( use_hermite_extrapolation ) then
      if ( hermite_order /= 3 .and. hermite_order /= 5 ) then
        print*, &
          'With Hermite extraplation the hermite_order has to be 3 or 5'
        stop
      end if
      if ( use_smooth_derivs_for_hermite ) then
        if ( extrap_order+1>nlevels ) then
          print*,'Error: accel_history initialized with extrap_order+1>nlevels'
          stop
        end if
      else
        if ( extrap_order+2>nlevels ) then
          print*,'With Hermite extrapolation nlevels>=extrap_order+2'
          stop
        end if
      end if
    else
      if ( extrap_order+1>nlevels ) then
        print*,'Error: accel_history initialized with extrap_order+1>nlevels'
        stop
      end if
    end if

    init_timelevels%nlevels = nlevels
    init_timelevels%extrap_order = extrap_order
    init_timelevels%nfilled = 0

    if ( .not. allocated(init_timelevels%at) ) then
      allocate ( init_timelevels%at(nlevels), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocation at'
        stop
      end if
      init_timelevels%at = 0.0_wp
    end if
    if ( .not. allocated(init_timelevels%ar) ) then
      allocate ( init_timelevels%ar(nlevels), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocation ar'
        stop
      end if
      init_timelevels%ar = 0.0_wp
    end if
    if ( .not. allocated(init_timelevels%aphi) ) then
      allocate ( init_timelevels%aphi(nlevels), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocation aphi'
        stop
      end if
      init_timelevels%aphi = 0.0_wp
    end if
    if ( .not. allocated(init_timelevels%dt) ) then
      allocate ( init_timelevels%dt(nlevels), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocation dt'
        stop
      end if
      init_timelevels%dt = 0.0_wp
    end if
    if ( .not. allocated(init_timelevels%t_rel) ) then
      allocate ( init_timelevels%t_rel(nlevels), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocation t_rel'
        stop
      end if
      init_timelevels%t_rel = 0.0_wp
    end if
    if ( .not. allocated(init_timelevels%xm) ) then
      allocate ( init_timelevels%xm(extrap_order+1,nlevels), &
                 stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocation xm'
        stop
      end if
      init_timelevels%xm = 0.0_wp
    end if
    if ( .not. allocated(init_timelevels%yv) ) then
      allocate ( init_timelevels%yv(nlevels), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocation yv'
        stop
      end if
      init_timelevels%yv = 0.0_wp
    end if
    if ( .not. allocated(init_timelevels%cv) ) then
      allocate ( init_timelevels%cv(extrap_order+1), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocation cv'
        stop
      end if
      init_timelevels%cv = 0.0_wp
    end if
    if ( .not. allocated(init_timelevels%covm) ) then
      allocate ( init_timelevels%covm(extrap_order+1,extrap_order+1), &
                 stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocation covm'
        stop
      end if
      init_timelevels%covm = 0.0_wp
    end if
  end procedure init_timelevels


  module procedure deallocate_timelevels

    implicit none

    integer :: allocation_status

    if ( allocated(this%at) ) then
      deallocate ( this%at, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for at'
        stop
      end if
    end if
    if ( allocated(this%ar) ) then
      deallocate ( this%ar, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for ar'
        stop
      end if
    end if
    if ( allocated(this%aphi) ) then
      deallocate ( this%aphi, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for aphi'
        stop
      end if
    end if
    if ( allocated(this%dt) ) then
      deallocate ( this%dt, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for dt'
        stop
      end if
    end if
    if ( allocated(this%t_rel) ) then
      deallocate ( this%t_rel, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for t_rel'
        stop
      end if
    end if
    if ( allocated(this%xm) ) then
      deallocate ( this%xm, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for xm'
        stop
      end if
    end if
    if ( allocated(this%yv) ) then
      deallocate ( this%yv, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for yv'
        stop
      end if
    end if
    if ( allocated(this%cv) ) then
      deallocate ( this%cv, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for cv'
        stop
      end if
    end if
    if ( allocated(this%covm) ) then
      deallocate ( this%covm, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Deallocation error for covm'
        stop
      end if
    end if

  end procedure deallocate_timelevels


  module procedure cycle_timelevels

    implicit none

    integer(ip) :: nf, nl

    nl = this%nlevels
    if ( this%nfilled > 0 ) then
      nf = min(this%nfilled,nl-1)
      this%at(nl-nf:nl-1) = this%at(nl-nf+1:nl)
      this%ar(nl-nf:nl-1) = this%ar(nl-nf+1:nl)
      this%aphi(nl-nf:nl-1) = this%aphi(nl-nf+1:nl)
      this%dt(nl-nf:nl-1) = this%dt(nl-nf+1:nl)
      this%t_rel(nl-nf:nl-1) = this%t_rel(nl-nf+1:nl)-dtime
    end if
    this%at(nl) = accel(1)
    this%ar(nl) = accel(2)
    this%aphi(nl) = accel(4)
    this%dt(nl) = dtime
    this%t_rel(nl) = 0.0_wp
    this%nfilled = min ( this%nfilled+1, nl )

  end procedure cycle_timelevels


  module procedure extrapolate

    use gsl_interface, only : multifit_linear
    use time_info, only : get_current_time
    use parameters, only : use_hermite_extrapolation, hermite_order, &
                           use_smooth_derivs_for_hermite
    use hermite, only : H3_interpolate, H5_interpolate

    implicit none

    integer(ip) :: nl, nf, order, ndat
    integer(ip) :: i, j, ierr
    real(wp) :: current_t, chisqr
    logical :: first = .true.
    real(wp), dimension(3,2) :: derivs

!    if (first) then
!      open(84,file='extrap.dat',status='replace',action='write')
!      first = .false.
!    end if
    nl = this%nlevels
    nf = this%nfilled
    order = max(min(this%extrap_order,nf-1),0)
    ndat = min(nl,nf)

    extrapolate = 0.0_wp

    if ( use_hermite_extrapolation .and. &
         ( .not. use_smooth_derivs_for_hermite) ) then
      if ( hermite_order==3 .and. nf<5 ) then
        return
      end if
      if ( hermite_order==5 .and. nf<7 ) then
        return
      end if
    end if

    if ( use_hermite_extrapolation .and. use_smooth_derivs_for_hermite ) then
!      if ( nf < order+2 ) then
      if ( nf < nl ) then
        return
      end if
    end if

    if ( .not. use_hermite_extrapolation ) then
      if ( nf < order+1 ) then
        return
      end if
    end if
!    if ( .not. use_hermite_extrapolation ) then
!      do i = 1, ndat
!        current_t = this%t_rel(nl-ndat+i)
!        do j = 1, order+1
!          this%xm(j,i) = current_t**(j-1)
!        end do
!      end do
!    end if

    ! First a^t
    if ( use_hermite_extrapolation ) then
      if ( use_smooth_derivs_for_hermite ) then
!        print*,'nl = ', nl
!        print*,'nf = ', nf
!        print*,'ndat = ', ndat
!        print*,'t = ', this%t_rel(nl-ndat+1:nl)
!        print*,'y = ', this%at(nl-ndat+1:nl)
        ! Evaluate the derivative at the next to last timelevel using the
        ! first ndat-1 tabulated values..
        derivs(:,1) = this%calc_smooth_derivs ( this%t_rel(nl-ndat+1:nl-1), &
                                                this%at(nl-ndat+1:nl-1), &
                                                this%t_rel(nl-1), order )
        ! Evaluate the derivative at the timelevel using the last ndat-1 
        !tabulated values..
        derivs(:,2) = this%calc_smooth_derivs ( this%t_rel(nl-ndat+2:nl), &
                                                this%at(nl-ndat+2:nl), &
                                                this%t_rel(nl), order )
        select case ( hermite_order )
        case (3)
          extrapolate(1,:) = H3_interpolate ( this%t_rel(nl-4:nl), &
                                              this%at(nl-4:nl), dtime, &
                                              dydx=derivs(2,:) )
        case (5)
          extrapolate(1,:) = H5_interpolate ( this%t_rel(nl-6:nl), &
                                              this%at(nl-6:nl), dtime, &
                                              dydx=derivs(2,:), &
                                              d2ydx2=derivs(3,:) )
        case default
          print*,'Error: Hermite interpolation with the wrong order'
          stop
        end select
      else
        select case ( hermite_order )
        case (3)
          extrapolate(1,:) = H3_interpolate ( this%t_rel(nl-4:nl), &
                                              this%at(nl-4:nl), dtime )
        case (5)
          extrapolate(1,:) = H5_interpolate ( this%t_rel(nl-6:nl), &
                                              this%at(nl-6:nl), dtime )
        case default
          print*,'Error: Hermite interpolation with the wrong order'
          stop
        end select
      end if
      print*,'a^t extrapolated = ', extrapolate(1,:)
    else
      extrapolate(1,:) = this%calc_smooth_derivs ( this%t_rel(nl-ndat+1:nl), &
                                                   this%at(nl-ndat+1:nl), &
                                                   dtime, order )
!      do i = 1, ndat
!        this%yv(i) = this%at(nl-ndat+i)
!      end do
!
!      call multifit_linear ( ndat, order+1, this%xm(:,1:ndat), &
!                             this%yv(1:ndat), this%cv, &
!                             this%covm, chisqr, ierr )
!
!      extrapolate(1,:) = (/ this%cv(order+1), order*this%cv(order+1), &
!                          order*(order-1)*this%cv(order+1) /)
!
!      do i = order, 1, -1
!        extrapolate(1,1) = extrapolate(1,1)*dtime + this%cv(i)
!        if (i>1) then
!          extrapolate(1,2) = extrapolate(1,2)*dtime + (i-1)*this%cv(i)
!        end if
!        if (i>2) then
!          extrapolate(1,3) = extrapolate(1,3)*dtime + (i-1)*(i-2)*this%cv(i)
!        end if
!      end do
    end if

    ! Second a^r
    if ( use_hermite_extrapolation ) then
      if ( use_smooth_derivs_for_hermite ) then
        ! Evaluate the derivative at the next to last timelevel using the
        ! first ndat-1 tabulated values..
        derivs(:,1) = this%calc_smooth_derivs ( this%t_rel(nl-ndat+1:nl-1), &
                                                this%ar(nl-ndat+1:nl-1), &
                                                this%t_rel(nl-1), order )
        ! Evaluate the derivative at the timelevel using the last ndat-1 
        !tabulated values..
        derivs(:,2) = this%calc_smooth_derivs ( this%t_rel(nl-ndat+2:nl), &
                                                this%ar(nl-ndat+2:nl), &
                                                this%t_rel(nl), order )
        select case ( hermite_order )
        case (3)
          extrapolate(2,:) = H3_interpolate ( this%t_rel(nl-4:nl), &
                                              this%ar(nl-4:nl), dtime, &
                                              dydx=derivs(2,:) )
        case (5)
          extrapolate(2,:) = H5_interpolate ( this%t_rel(nl-6:nl), &
                                              this%ar(nl-6:nl), dtime, &
                                              dydx=derivs(2,:), &
                                              d2ydx2=derivs(3,:) )
        case default
          print*,'Error: Hermite interpolation with the wrong order'
          stop
        end select
      else
        select case ( hermite_order )
        case (3)
          extrapolate(2,:) = H3_interpolate ( this%t_rel(nl-4:nl), &
                                              this%ar(nl-4:nl), dtime )
        case (5)
          extrapolate(2,:) = H5_interpolate ( this%t_rel(nl-6:nl), &
                                              this%ar(nl-6:nl), dtime )
        case default
          print*,'Error: Hermite interpolation with the wrong order'
          stop
        end select
      end if
!      write(84,*) get_current_time(), this%dt, this%t_rel, dtime, this%ar, extrapolate(2,:)
    else
      extrapolate(2,:) = this%calc_smooth_derivs ( this%t_rel(nl-ndat+1:nl), &
                                                   this%ar(nl-ndat+1:nl), &
                                                   dtime, order )
!      do i = 1, ndat
!        this%yv(i) = this%ar(nl-ndat+i)
!      end do
!
!      call multifit_linear ( ndat, order+1, this%xm(:,1:ndat), &
!                             this%yv(1:ndat), this%cv, &
!                             this%covm, chisqr, ierr )
!      write(84,*) get_current_time(), this%dt, this%t_rel, dtime, this%ar, this%cv
!
!      extrapolate(2,:) = (/ this%cv(order+1), order*this%cv(order+1), &
!                            order*(order-1)*this%cv(order+1) /)
!
!      do i = order, 1, -1
!        extrapolate(2,1) = extrapolate(2,1)*dtime + this%cv(i)
!        if (i>1) then
!          extrapolate(2,2) = extrapolate(2,2)*dtime + (i-1)*this%cv(i)
!        end if
!        if (i>2) then
!          extrapolate(2,3) = extrapolate(2,3)*dtime + (i-1)*(i-2)*this%cv(i)
!        end if
!      end do
    end if
    
    ! Third a^phi
    if ( use_hermite_extrapolation ) then
      if ( use_smooth_derivs_for_hermite ) then
        ! Evaluate the derivative at the next to last timelevel using the
        ! first ndat-1 tabulated values..
        derivs(:,1) = this%calc_smooth_derivs ( this%t_rel(nl-ndat+1:nl-1), &
                                                this%aphi(nl-ndat+1:nl-1), &
                                                this%t_rel(nl-1), order )
        ! Evaluate the derivative at the timelevel using the last ndat-1 
        !tabulated values..
        derivs(:,2) = this%calc_smooth_derivs ( this%t_rel(nl-ndat+2:nl), &
                                                this%aphi(nl-ndat+2:nl), &
                                                this%t_rel(nl), order )
        select case ( hermite_order )
        case (3)
          extrapolate(4,:) = H3_interpolate ( this%t_rel(nl-4:nl), &
                                              this%aphi(nl-4:nl), dtime, &
                                              dydx=derivs(2,:) )
        case (5)
          extrapolate(4,:) = H5_interpolate ( this%t_rel(nl-6:nl), &
                                              this%aphi(nl-6:nl), dtime, &
                                              dydx=derivs(2,:), &
                                              d2ydx2=derivs(3,:) )
        case default
          print*,'Error: Hermite interpolation with the wrong order'
          stop
        end select
      else
        select case ( hermite_order )
        case (3)
          extrapolate(4,:) = H3_interpolate ( this%t_rel(nl-4:nl), &
                                              this%aphi(nl-4:nl), dtime )
        case (5)
          extrapolate(4,:) = H5_interpolate ( this%t_rel(nl-6:nl), &
                                              this%aphi(nl-6:nl), dtime )
        case default
          print*,'Error: Hermite interpolation with the wrong order'
          stop
        end select
      end if
    else
      extrapolate(4,:) = this%calc_smooth_derivs ( this%t_rel(nl-ndat+1:nl), &
                                                   this%aphi(nl-ndat+1:nl), &
                                                   dtime, order )
!      do i = 1, ndat
!        this%yv(i) = this%aphi(nl-ndat+i)
!      end do
!
!      call multifit_linear ( ndat, order+1, this%xm(:,1:ndat), &
!                             this%yv(1:ndat), this%cv, &
!                             this%covm, chisqr, ierr )
!
!      extrapolate(4,:) = (/ this%cv(order+1), order*this%cv(order+1), &
!                            order*(order-1)*this%cv(order+1) /)
!
!      do i = order, 1, -1
!        extrapolate(4,1) = extrapolate(4,1)*dtime + this%cv(i)
!        if (i>1) then
!          extrapolate(4,2) = extrapolate(4,2)*dtime + (i-1)*this%cv(i)
!        end if
!        if (i>2) then
!          extrapolate(4,3) = extrapolate(4,3)*dtime + (i-1)*(i-2)*this%cv(i)
!        end if
!      end do
    end if

  end procedure extrapolate


  module procedure calc_smooth_derivs

    use gsl_interface, only : multifit_linear

    implicit none

    integer(ip) :: lorder, ndat, i, j, ierr
    real(wp) :: current_x, chisqr

    ndat = size(x)
    if ( ndat<order+1 ) then
      lorder = ndat-1
    else
      lorder =  order
    end if

!    print*,'calc_smooth_derivs : ', ndat, lorder
!    print*,'dx = ', dx
!    print*,'t = ', x
!    print*,'y = ', y

    do i = 1, ndat 
      current_x = x(i)
      do j = 1, lorder+1
        this%xm(j,i) = current_x**(j-1)
      end do
      this%yv(i) = y(i)
    end do

    call multifit_linear ( ndat, lorder+1, this%xm(:,1:ndat), &
                           this%yv(1:ndat), this%cv, &
                           this%covm, chisqr, ierr )

!    print*,'coeffs = ', this%cv

    calc_smooth_derivs = (/ this%cv(lorder+1), lorder*this%cv(lorder+1), &
                          lorder*(lorder-1)*this%cv(lorder+1) /)

    do i = lorder, 1, -1
      calc_smooth_derivs(1) = calc_smooth_derivs(1)*dx + this%cv(i)
      if (i>1) then
        calc_smooth_derivs(2) = calc_smooth_derivs(2)*dx + (i-1)*this%cv(i)
      end if
      if (i>2) then
        calc_smooth_derivs(3) = calc_smooth_derivs(3)*dx &
                              + (i-1)*(i-2)*this%cv(i)
      end if
    end do
            
  end procedure calc_smooth_derivs

end submodule acceleration_history_implementation
