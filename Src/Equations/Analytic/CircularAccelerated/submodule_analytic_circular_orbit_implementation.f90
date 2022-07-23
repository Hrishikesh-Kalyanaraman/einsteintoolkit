submodule(analytic_circular_orbit) analytic_circular_orbit_implementation

  implicit none

contains

  module procedure co_init

    use parameters, only : p_orb, mass, accel_amp, accel_sigma, accel_t0, &
                           omega_ratio, phi_initial

    implicit none

    real(wp) :: f, ut, uphi
    integer :: allocation_status

    print*,'Initializing circular orbit'

    if ( .not. allocated(this%ename) ) then
      allocate ( character(14) :: this%ename, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating scalar_schw equation name'
        stop
      end if
    end if
    this%ename = 'Circular orbit'

    this%r = p_orb*mass
    this%omega = omega_ratio*sqrt(mass/this%r**3)
    this%amp = accel_amp
    this%sigma = accel_sigma
    this%t0 = accel_t0
    f = 1.0_wp-2.0_wp*mass/this%r
    ut = 1.0_wp/sqrt(f-this%r**2*this%omega**2)
    uphi = this%omega*ut
    this%En = f*ut
    this%Lz = this%r**2*uphi
    this%chiomega = mass*sqrt(this%r/mass-6.0_wp)/this%r**2
    this%phi_initial = phi_initial

    this%io_id = -1
    call this%save_globals_1 ( )

  end procedure co_init


! This is empty as nothing needs to be done.
  module procedure co_rhs

  end procedure co_rhs


! This is empty as nothing needs to be done.
  module procedure co_set_to_zero

  end procedure co_set_to_zero


! This is empty as nothing needs to be done.
  module procedure co_update_vars

  end procedure co_update_vars


! Save the current orbit in the orbit_base and self_force_base variables.
  module procedure co_save_globals_1

    use orbit_base, only : orbit_info
    use self_force_base, only : sf
    use time_info, only : get_current_qtime
    use parameters, only : use_constant_acceleration, mass, accel_t0, &
                           accel_sigma, accel_amp, use_gaussian_acceleration, &
                           q_mass, q_charge
    use accelerated_circular_orbit, only : circ_accel

    implicit none

    real(qp) :: phi_tmp, chi_tmp
    real(wp) :: En, Lz, accel_amp_actual
    real(wp) :: phi, chi, ar, aphi, dardt, daphidt, d2ardt2, d2aphidt2
    real(qp), parameter :: piq2 = 2.0_qp*acos(-1.0_qp)
    real(wp), parameter :: sqrteo2 = sqrt(exp(1.0_wp)*0.5_wp)
    real(wp) :: ffac, gurr, gutt, guphiphi, at, ut, ur, uphi, ft, fr, fphi
    real(wp) :: at2, ar2, aphi2 

    if ( use_constant_acceleration ) then
      phi_tmp = this%omega * get_current_qtime ( ) + this%phi_initial
      phi_tmp = mod(phi_tmp, piq2)
      phi = real(phi_tmp,wp)
      chi_tmp = this%chiomega * get_current_qtime ( )
      chi_tmp = mod(chi_tmp, piq2)
      chi = real(chi_tmp,wp)
      ar = acc_r ( this%omega, mass, this%r )

      call orbit_info%set_orbit ( this%r, phi, 0.0_wp, this%En, this%Lz, chi )
      call sf%set_accel ( 0.0_wp, ar, 0.0_wp, 0.0_wp, 0.0_wp )
      call sf%set_daccel_dt ( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )
      call sf%set_d2accel_dt2 ( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )
    else
      if ( use_gaussian_acceleration ) then
        accel_amp_actual = accel_amp
      else
        accel_amp_actual = accel_amp*sqrteo2*sqrt(mass/this%r**3)*accel_sigma
      end if
      call circ_accel ( get_current_qtime ( ), accel_t0, this%r, mass, &
                        accel_sigma, accel_amp_actual, phi, this%En, this%Lz, &
                        ar, aphi, dardt, daphidt, d2ardt2, d2aphidt2 )
      ffac = 1.0_wp - 2.0_wp*mass/this%r
      ut = this%En/ffac
      ur = 0.0_wp
      uphi = this%Lz/this%r**2

      gurr = ffac
      gutt = -1.0_wp/gurr
      guphiphi = 1.0_wp/this%r**2

      at = -gutt*(aphi*gurr*uphi+ar*guphiphi*ur)/(guphiphi*gurr*ut)

      call orbit_info%set_orbit ( this%r, phi, 0.0_wp, this%En, this%Lz, &
                                  0.0_wp )
      call sf%set_accel ( at, ar, 0.0_wp, aphi, 0.0_wp )
      call sf%set_daccel_dt ( 0.0_wp, dardt, 0.0_wp, daphidt )
      call sf%set_d2accel_dt2 ( 0.0_wp, d2ardt2, 0.0_wp, d2aphidt2 )

      ft = q_mass/q_charge*(ut+at)/gutt
      fr = q_mass/q_charge*(ur+ar)/gurr
      fphi = q_mass/q_charge*(uphi+aphi)/guphiphi
      call sf%set_force ( ft, fr, 0.0_wp, fphi )

!      at2 = q_charge/q_mass*((gutt+ut**2)*ft+ut*(ur*fr+uphi*fphi))
!      ar2 = q_charge/q_mass*((gurr+ur**2)*fr+ur*(ut*ft+uphi*fphi))
!      aphi2 = q_charge/q_mass*((guphiphi+uphi**2)*fphi+uphi*(ut*ft+ur*fr))
!      print*,'Acc2 = ', at2, ar2, aphi2
    end if
  
  end procedure co_save_globals_1


! This procedure is empty as there is no dependence on other saved variables.
  module procedure co_save_globals_2

  end procedure co_save_globals_2


! This is empty as nothing needs to be done.
  module procedure co_load_globals

  end procedure co_load_globals


  module procedure co_output

    use output_base
    use time_info

    implicit none

    character(len=14) :: filename = 'circ_orbit.asc'
    integer(ip) :: io_id
    real(wp) :: time

    io_id = this%io_id
    if (io_id < 0 ) then
      io_id = next_available_io_id ()
      this%io_id = io_id
      print*,'Opening ', filename, ' with id ', io_id
      open(io_id, file=filename, status='replace', action='write')
      write(io_id,*) '#1: time, 2: r, 3: phi, 4: chi, 5: En, 6: Lz'
    end if

    time = get_current_time ( )
    write(io_id,'(*(es23.15e3,1x))') time, this%r, &
                                     time*this%omega+this%phi_initial, &
                                     time*this%chiomega, this%En, this%Lz
  end procedure co_output


  module procedure co_print_data

  end procedure co_print_data


  function acc_r ( omega, mass, r )
    real(wp), intent(in) :: omega, mass, r
    real(wp) :: acc_r

    acc_r = (r - 2.0_wp*mass)*(mass - r**3*omega**2) &
            /(r**2*(r - 2.0_wp*mass - r**3*omega**2))
  end function acc_r

end submodule analytic_circular_orbit_implementation
