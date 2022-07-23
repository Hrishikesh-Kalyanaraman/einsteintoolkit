submodule(osculating_schwarzschild) osculating_schwarzschild_implementation
!! The implementation of the interface defined in [[osculating_schwarzschild]].

  implicit none

contains

  module procedure osc_schw_init

    use parameters, only : p_orb, ecc, mass, phi_initial, q_mass, &
                           use_analytic_driver
    use all_integrators, only : mol_ntmp
    use orbit_base, only : orbit_info

    implicit none

    real(wp) :: r, phi, ur
    integer :: allocation_status

    this%nvars = 6
    if ( .not. allocated(this%ename) ) then
      allocate ( character(24) :: this%ename, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating osc_schw equation name'
        stop
      end if
    end if
    this%ename = 'Schwarzschild_osculating'

    if ( .not. allocated(this%var_data) ) then
      allocate ( this%var_data(this%nvars), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating osc_schw variables'
        stop
      end if
    end if

    if ( .not. allocated(this%rhs_data) ) then
      allocate ( this%rhs_data(this%nvars), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating osc_schw rhs-variables'
        stop
      end if
    end if

    this%ntmp = mol_ntmp ( )

    if ( .not. allocated(this%tmp_data) ) then
      allocate ( this%tmp_data(this%nvars,this%ntmp), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating osc_schw tmp-variables'
        stop
      end if
    end if

    this%mass = mass
    associate ( p => this%var_data(3), &
                e => this%e, &
                En => this%En, &
                Lz => this%Lz, &
                r => this%r, &
                chi => this%var_data(1), &
                phi => this%var_data(2), &
                alpha => this%var_data(4), &
                beta => this%var_data(5), &
                qmass => this%var_data(6), &
                w => this%w )
      p = p_orb / mass
      e = ecc
      print*,'Initializing orbit with p = ', p, ' and e = ', e
      if (use_analytic_driver) then
        call orbit_info%get_orbit ( r, phi, ur, this%En, this%Lz )
      else
        En = sqrt((p -2.0_wp*(1.0_wp+e))*(p-2.0_wp*(1.0_wp-e)) &
                   /(p*(p-3.0_wp-e**2)))
        Lz = p*mass/sqrt(p-3.0_wp-e**2)
      end if
      chi = acos(-1.0_wp) ! apastron
      phi = phi_initial
      w = rzero
      alpha = e*sin(w)
      beta = e*cos(w)
      qmass = q_mass
    end associate

    this%io_id = -1

    call this%save_globals_1 ( )

  end procedure osc_schw_init


  module procedure osc_schw_rhs

    implicit none

    real(wp) :: psi, omega, betap, alphap, pp, tp, phip, sinchi, coschi, dtaudt

    associate ( chi => this%var_data(1), &
                p => this%var_data(3), &
                alpha => this%var_data(4), &
                beta => this%var_data(5), &
                dchidt => this%rhs_data(1), &
                d2chidt2 => this%d2chidt2, &
                dphidt => this%rhs_data(2), &
                dpdt => this%rhs_data(3), &
                dalphadt => this%rhs_data(4), &
                dbetadt => this%rhs_data(5), &
                dqmassdt => this%rhs_data(6), &
                En => this%En, &
                accel => this%accel, &
                udota => this%udota, &
                mass => this%mass, &
                r => this%r, &
                drdt => this%drdt, &
                d2rdt2 => this%d2rdt2 )
      sinchi = sin(chi)
      coschi = cos(chi)
      psi = alpha*sinchi
      omega = beta*coschi

      betap = (mass*p**2*(-3 - alpha**2 - beta**2 + p)*accel(2)* &
           (2*alpha*(1 + omega) + (-2*(3 + beta**2) + p)*Sin(chi)))/ &
         ((-4*(-9 + alpha**2 + beta**2) - 12*p + p**2)*(1 + omega + psi)**2) &
         + (mass**2*p**2.5_wp*(-3 - alpha**2 - beta**2 + p)*accel(4)* &
           (beta*(12 + 4*(alpha**2 + beta**2) - 10*p + p**2) + &
             (2*(-3 + p) + (-6 + p)*(omega + psi) - 2*(omega + psi)**2)* &
              (-2*beta*(omega + psi) + (-6 + p)*Cos(chi)) - &
             4*alpha*(alpha*beta*Cos(2*chi) + &
                ((alpha**2 - beta**2)*Sin(2*chi))/2.0_wp)))/ &
         ((-4*(alpha**2 + beta**2) + (-6 + p)**2)*(1 + omega + psi)**4* &
           Sqrt(-6 + p - 2*(omega + psi)))
      alphap = -((mass*p**2*(-3 - alpha**2 - beta**2 + p)* &
             (2*beta*(1 + psi) + (-2*(3 + alpha**2) + p)*Cos(chi))*accel(2))/ &
           ((-4*(-9 + alpha**2 + beta**2) - 12*p + p**2)*(1 + omega + psi)**2) &
           ) + (mass**2*p**2.5_wp*(-3 - alpha**2 - beta**2 + p)*accel(4)* &
           (alpha*(12 + 4*(alpha**2 + beta**2) - 10*p + p**2) +  &
             (2*(-3 + p) + (-6 + p)*(omega + psi) - 2*(omega + psi)**2)* &
              (-2*alpha*(omega + psi) + (-6 + p)*Sin(chi)) +  &
             4*beta*(alpha*beta*Cos(2*chi) +  &
                ((alpha**2 - beta**2)*Sin(2*chi))/2.0_wp)))/ &
         ((-4*(alpha**2 + beta**2) + (-6 + p)**2)*(1 + omega + psi)**4* &
           Sqrt(-6 + p - 2*(omega + psi)))
      pp = (2*mass**2*p**3.5_wp*(-3 - alpha**2 - beta**2 + p)* &
           Sqrt(-6 + p - 2*(omega + psi))*(-3 + p - (omega + psi)**2)*accel(4))/ &
         ((-4*(alpha**2 + beta**2) + (-6 + p)**2)*(1 + omega + psi)**4) +  &
        (2*mass*p**3*(-3 - alpha**2 - beta**2 + p)*accel(2)* &
           (alpha*Cos(chi) - beta*Sin(chi)))/ &
         ((-4*(-9 + alpha**2 + beta**2) - 12*p + p**2)*(1 + omega + psi)**2)
      tp = (mass*p**2*Sqrt((-4*(alpha**2 + beta**2) + (-2 + p)**2)/ &
            (p - 2*(3 + omega + psi))))/ &
        ((1 + omega + psi)**2*(p - 2*(1 + omega + psi)))
! Comment this out for now.
!      if (use_constant_acceleration) then
!        tp = dtdchifac * tp
!      end if
      phip = Sqrt(p/(p - 2*(3 + omega + psi)))

      dchidt = 1.0_wp/tp
      dbetadt = betap*dchidt
      dalphadt = alphap*dchidt
      dpdt = pp*dchidt
      dphidt = phip*dchidt
      r = (mass*p)/(1.0_wp + beta*Cos(chi) + alpha*Sin(chi))
      dtaudt = (1.0_wp - 2.0_wp*mass/r) / this%En
      dqmassdt = -dtaudt*udota

      d2chidt2 = ((1 + beta*Cos(chi) + alpha*Sin(chi))* &
          (-((-8*(alpha*dalphadt + beta*dbetadt) + 2*dpdt*(-2 + p))*p* &
               (1 + beta*Cos(chi) + alpha*Sin(chi))* &
               (p - 2*(1 + beta*Cos(chi) + alpha*Sin(chi)))* &
               (p - 2*(3 + beta*Cos(chi) + alpha*Sin(chi)))) - &
            4*dpdt*(-4*(-1 + alpha**2 + beta**2) - 4*p + p**2)* &
             (1 + beta*Cos(chi) + alpha*Sin(chi))* &
             (p - 2*(1 + beta*Cos(chi) + alpha*Sin(chi)))* &
             (p - 2*(3 + beta*Cos(chi) + alpha*Sin(chi))) +  &
            4*p*(-4*(-1 + alpha**2 + beta**2) - 4*p + p**2)* &
             (p - 2*(1 + beta*Cos(chi) + alpha*Sin(chi)))* &
             (p - 2*(3 + beta*Cos(chi) + alpha*Sin(chi)))* &
             (dbetadt*Cos(chi) + dalphadt*Sin(chi) +  &
               dchidt*(alpha*Cos(chi) - beta*Sin(chi))) + &
            p*(-4*(-1 + alpha**2 + beta**2) - 4*p + p**2)* &
             (1 + beta*Cos(chi) + alpha*Sin(chi))* &
             (p - 2*(1 + beta*Cos(chi) + alpha*Sin(chi)))* &
             (dpdt - 2*(dbetadt*Cos(chi) + dalphadt*Sin(chi) +  &
                  dchidt*(alpha*Cos(chi) - beta*Sin(chi)))) +  &
            2*p*(-4*(-1 + alpha**2 + beta**2) - 4*p + p**2)* &
             (1 + beta*Cos(chi) + alpha*Sin(chi))* &
             (p - 2*(3 + beta*Cos(chi) + alpha*Sin(chi)))* &
             (dpdt - 2*(dbetadt*Cos(chi) + dalphadt*Sin(chi) +  &
                  dchidt*(alpha*Cos(chi) - beta*Sin(chi))))))/ &
        (2*mass*p**3*(-4*(-1 + alpha**2 + beta**2) - 4*p + p**2)** &
           1.5_wp*Sqrt(p - 2*(3 + beta*Cos(chi) + alpha*Sin(chi))))

      drdt = -((dchidt*mass*p*(alpha*Cos(chi) - beta*Sin(chi)))/ &
          (1 + beta*Cos(chi) + alpha*Sin(chi))**2)
      d2rdt2 = -((d2chidt2*mass*p*(alpha*Cos(chi) - beta*Sin(chi)))/ &
           (1 + beta*Cos(chi) + alpha*Sin(chi))**2) -  &
        (dchidt*mass*p*(dalphadt*Cos(chi) - beta*dchidt*Cos(chi) -  &
             dbetadt*Sin(chi) - alpha*dchidt*Sin(chi)))/ &
         (1 + beta*Cos(chi) + alpha*Sin(chi))**2 +  &
        (2*dchidt*mass*p*(alpha*Cos(chi) - beta*Sin(chi))* &
           (dbetadt*Cos(chi) + alpha*dchidt*Cos(chi) + dalphadt*Sin(chi) -  &
             beta*dchidt*Sin(chi)))/(1 + beta*Cos(chi) + alpha*Sin(chi))**3
    end associate
  end procedure osc_schw_rhs


! Save the current state of the orbit in the orbit_base variables.
  module procedure osc_schw_save_globals_1

    use time_info, only: get_current_time
    use orbit_base

!    print*,'Osc Schw save globals called at time: ', get_current_time()
    call this%rhs ( )
    call this%calc_dependent ( )
    call orbit_info%set_orbit ( this%r, this%var_data(2), &
                                this%ur, this%En, this%Lz, this%var_data(1) )

    call tdc_info%set_tdc ( this%r, this%drdt, this%d2rdt2 )
  end procedure osc_schw_save_globals_1


! This procedure is empty since no saved variables are needed.
  module procedure osc_schw_save_globals_2

  end procedure osc_schw_save_globals_2


  module procedure osc_schw_load_globals

    use parameters, only : mass, q_charge, q_mass, use_analytic_driver, &
                           a_timelevels
    use time_info, only : get_current_time, get_current_dtime
    use self_force_base, only : sf
    use acceleration_history, only : ah
    use all_integrators, only : mol_int

    implicit none

    real(wp), dimension(4) :: xp, up
    real(wp) :: ffac, gurr, gutt, guthth, guphiphi, udota
    real(wp), parameter :: pi = acos (-1.0_wp)
    real(qp) :: Enq, Lzq, rq, phiq, ffacq, chiq, alphaq, betaq, pq
    real(qp), dimension(4) :: upq

!    print*,'Osc Schw load globals called at time: ', get_current_time()
! The time window factor for the force is applied in get_force, so the
! acceleration that is stored later in this routine will have the same
! time window factor applied.
    call sf%get_force ( this%force )

    associate ( r => this%r, &
                phi => this%var_data(2), &
                En => this%En, &
                Lz => this%Lz, &
                chi => this%var_data(1), &
                p => this%var_data(3), &
                alpha => this%var_data(4), &
                beta => this%var_data(5), &
                qmass => this%var_data(6), &
                force => this%force, &
                accel => this%accel, &
                udota => this%udota )
      xp(1) = get_current_time ( )
      xp(2) = r
      xp(3) = 0.5_wp*pi
      xp(4) = phi

      Enq = real(En,qp)
      Lzq = real(Lz,qp)
      rq = real(r,qp)
      phiq = real(phi,qp)
      ffacq = 1.0_qp - 2.0_qp*real(mass,qp)/rq
      chiq = real(chi,qp)
      alphaq = real(alpha,qp)
      betaq = real(beta,qp)
      pq = real(p,qp)

      upq(1) = Enq/ffacq
      upq(2) = (betaq*sin(chiq)-alphaq*cos(chiq))*sqrt(pq-6-2*(alphaq &
                *sin(chiq)+betaq*cos(chiq)))/sqrt(pq*(pq-3-alphaq**2-betaq**2))
      upq(3) = 0.0_qp
      upq(4) = Lzq/rq**2

      gurr = real(ffacq,wp)
      gutt = -1.0_wp/gurr
      guthth = 1.0_wp/r**2
      guphiphi = guthth

      up = real(upq,wp)

      accel(1) = (gutt + up(1)**2)*force(1) + up(1)*( up(2)*force(2) &
                                                    + up(4)*force(4) )
      accel(2) = (gurr + up(2)**2)*force(2) + up(2)*( up(1)*force(1) &
                                                    + up(4)*force(4) )
      accel(3) = 0.0_wp
      accel(4) = (guphiphi + up(4)**2)*force(4) &
                           + up(4)*( up(1)*force(1) + up(2)*force(2) )

      accel = q_charge/q_mass*accel

      if (use_analytic_driver) then
        udota = 0.0_wp
      else
        udota = q_charge*dot_product(up, force)
      end if
!      print*,'udota = ', udota

      call sf%set_accel ( accel(1), accel(2), accel(3), accel(4), udota )

      if ( a_timelevels > 0 .and. mol_int%complete_step ( ) ) then
        call ah%cycle_timelevels ( accel, get_current_dtime ( ) )
      end if
    end associate

  end procedure osc_schw_load_globals


  module procedure osc_schw_output

    use output_base
    use time_info

    implicit none

    character(len=13) :: filename = 'osc_orbit.asc'
    integer(ip) :: io_id
    real(wp) :: time

    io_id = this%io_id
    if (io_id < 0 ) then
      io_id = next_available_io_id ()
      this%io_id = io_id
      print*,'Opening ', filename, ' with id ', io_id
      open(io_id, file=filename, status='replace', action='write')
      write(io_id,*) '#1: time, 2: chi, 3: phi, 4: p, 5: alpha, 6: beta, 7: qmass, 8: dchidt, 9:dphidt, 10: dpdt, 11: dalphadt, 12: dbetadt, 13: dqmassdt, 14: e, 15: r, 16: drdt, 17: d2rdt2, 18: En, 19: Lz'
    end if

    time = get_current_time ( )
    call this%rhs ( )
    call this%calc_dependent ( )
    write(io_id,'(*(es23.15e3,1x))') time, this%var_data, this%rhs_data, &
                                     this%e, this%r, this%drdt, this%d2rdt2, &
                                     this%En, this%Lz
  end procedure osc_schw_output


  module procedure close_osc_schw

    deallocate ( this%var_data, this%rhs_data, this%tmp_data )
    if ( this%io_id>=0 ) then
      close (this%io_id)
    end if

  end procedure close_osc_schw


  module procedure calc_dependent

    use parameters, only : use_analytic_driver
    use orbit_base, only : orbit_info

    implicit none

    ! We do this in quad precision to avoid unnecessary noise in ur.
    real(qp) :: eq, pq, rposq, Enq, Lzq, urq, massq, alphaq, betaq, chiq, wq

    chiq = real(this%var_data(1),qp)
    pq = real(this%var_data(3),qp)
    alphaq = real(this%var_data(4),qp)
    betaq = real(this%var_data(5),qp)
    massq = real(this%mass,qp)

    wq = atan2(alphaq,betaq) 
    eq = sqrt(alphaq**2+betaq**2)
    Enq = sqrt((pq -2.0_qp*(1.0_qp+eq))*(pq-2.0_qp*(1.0_qp-eq)) &
               /(pq*(pq-3.0_qp-eq**2)))
    Lzq = pq*massq/sqrt(pq-3.0_qp-eq**2)
    rposq = pq*massq/(1.0_qp+alphaq*sin(chiq)+betaq*cos(chiq))
    urq = sign(1.0_qp,sin(chiq-wq)) &
              *sqrt(abs(Enq**2 - (1.0_qp + (Lzq/rposq)**2) &
                                *(1.0_qp - 2.0_qp*massq/rposq)))
    
    this%w = real(wq,wp)
    this%e = real(eq,wp)
    if ( .not. use_analytic_driver ) then
      this%En = real(Enq,wp)
      this%Lz = real(Lzq,wp)
      this%ur = real(urq,wp)
    else
      call orbit_info%get_orbit(this%r, this%var_data(2), this%ur, this%En, &
                             this%Lz)
    end if

  end procedure calc_dependent

end submodule osculating_schwarzschild_implementation
