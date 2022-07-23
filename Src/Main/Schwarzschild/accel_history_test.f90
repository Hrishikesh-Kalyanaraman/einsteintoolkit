program accel_history_test

  use self_force_base
  use acceleration_history

  real(wp), dimension(4) :: acc
  real(wp), dimension(4,3) :: eacc

  ah = accel_history ( 5, 2 )

  call sf%set_accel ( 1.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp )
  call sf%get_accel ( acc )
  call ah%cycle_timelevels ( acc, 0.5_wp )

  call sf%set_accel ( 1.75_wp, 1.25_wp, 0.0_wp, 0.75_wp, 0.0_wp )
  call sf%get_accel ( acc )
  call ah%cycle_timelevels ( acc, 1.0_wp )

  call sf%set_accel ( 4.75_wp, 0.25_wp, 0.0_wp, 1.75_wp, 0.0_wp )
  call sf%get_accel ( acc )
  call ah%cycle_timelevels ( acc, 1.5_wp )

  call sf%set_accel ( 13.0_wp, -5.0_wp, 0.0_wp, 7.0_wp, 0.0_wp )
  call sf%get_accel ( acc )
  
  call ah%cycle_timelevels ( acc, 1.0_wp )

  call sf%set_accel ( 21.0_wp, -11.0_wp, 0.0_wp, 13.0_wp, 0.0_wp )
  call sf%get_accel ( acc )
  
  call ah%cycle_timelevels ( acc, 1.0_wp )

  print*,'nlevels = ', ah%nlevels
  print*,'nfilled = ', ah%nfilled
  print*,'extrap_order = ', ah%extrap_order

  print*,'at = ', ah%at
  print*,'ar = ', ah%ar
  print*,'aphi = ', ah%aphi
  print*,'dt = ', ah%dt
  print*,'t_rel = ', ah%t_rel

  eacc = ah%extrapolate ( 0.5_wp )

  print*,'at(1.0) = ', eacc(1,:)
  print*,'ar(1.0) = ', eacc(2,:)
  print*,'aphi(1.0) = ', eacc(4,:)
  
end program accel_history_test
