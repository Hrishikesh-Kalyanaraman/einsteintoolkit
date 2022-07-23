program accel_test

  use parameters
  use time_info
  use all_integrators
  use osculating_schwarzschild
  use geodesic_schwarzschild
  use numerics
  use analytic_circular_orbit

  implicit none

  class(equation_pointer), dimension(:), allocatable :: eqs
  type(geod_schw), target :: geod_orbit
  type(osc_schw), target :: osc_orbit
  type(circular_orbit), target :: driver
  integer(ip) :: k

  call read_parameters ( )

  call init_time ( t_initial )

  allocate ( eqs(2) )
  eqs(1)%p => driver
  if (use_osculating_orbit) then
    eqs(2)%p => osc_orbit
  else
    eqs(2)%p => geod_orbit
  endif

  call driver%init ( )

  call eqs(2)%p%init ( )

  call choose_integrator ( )

  call mol_int%init(eqs)

  call set_dtime ( 1.0e-2_wp )

  call driver%output ( )
  call eqs(2)%p%output ( )

  k=0
  do
    if ( get_current_time ( ) > t_final ) exit
    k = k+1

    call mol_int%step ( )

    if (out0d_every>0 .and. mod(k,out0d_every)==0) then
      call driver%output ( )
      call eqs(2)%p%output ( )
    end if
  end do

  print*, 'Program ended after ', k, ' iterations.'
end program accel_test
