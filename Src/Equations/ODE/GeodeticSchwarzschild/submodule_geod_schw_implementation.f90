submodule(geodesic_schwarzschild) geodesic_schwarzschild_implementation
!! The implementation of the interfaces defined in [[geodesic_schwarzschild]].

  implicit none

contains

  module procedure geod_schw_init

    use parameters, only : p_orb, ecc, mass, phi_initial, q_mass, &
                           use_analytic_driver
    use all_integrators, only : mol_ntmp
    use orbit_base, only : orbit_info

    implicit none

    integer :: allocation_status

    this%nvars = 7
    if ( .not. allocated(this%ename) ) then
      allocate ( character(22) :: this%ename, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating geod_schw equation name'
        stop
      end if
    end if
    this%ename = 'Schwarzschild_geodesic'

    if ( .not. allocated(this%var_data) ) then
      allocate ( this%var_data(this%nvars), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating geod_schw variables'
        stop
      end if
    end if

    if ( .not. allocated(this%rhs_data) ) then
      allocate ( this%rhs_data(this%nvars), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating geod_schw rhs-variables'
        stop
      end if
    end if

    this%ntmp = mol_ntmp ( )

    if ( .not. allocated(this%tmp_data) ) then
      allocate ( this%tmp_data(this%nvars,this%ntmp), stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating geod_schw tmp-variables'
        stop
      end if
    end if

    associate ( p => this%p, &
                e => this%e, &
                En => this%En, &
                Lz => this%Lz, &
                w => this%w, &
                r => this%var_data(1), &
                phi => this%var_data(2), &
                ut => this%var_data(3), &
                ur => this%var_data(4), &
                uphi => this%var_data(5), &
                qmass => this%var_data(6), &
                chi => this%var_data(7) )
      p = p_orb / mass
      e = ecc

      if (use_analytic_driver) then
        call orbit_info%get_orbit ( r, phi, ur, this%En, this%Lz )
      else
        En = sqrt((p -2.0_wp*(1.0_wp+e))*(p-2.0_wp*(1.0_wp-e)) &
                    /(p*(p-3.0_wp-e**2)))
        Lz = p*mass/sqrt(p-3.0_wp-e**2)
      end if
      r = p*mass/(1.0_wp-e)
      phi = phi_initial
      chi = acos(-1.0_wp) ! pi
      w = 0.0_wp
      ut = (En*r)/(r - 2*mass )
      uphi = Lz/r**2
! Currently starts at apastron
      ur = 0.0_wp
!      ur = sqrt(((2.0_wp*mass - r)*(r + r**3*uphi**2 + 2.0_wp*mass*ut**2 &
!                                   - r*ut**2))/r**2)
      qmass = q_mass
    end associate
    
    this%io_id = -1

    call this%save_globals_1 ( )

  end procedure geod_schw_init


  module procedure geod_schw_rhs

    use parameters, only : mass

    implicit none

    real(wp) :: utp, urp, uphip, utinv, urpp, ecoschi

    associate ( r => this%var_data(1), &
                phi => this%var_data(2), &
                ut => this%var_data(3), &
                ur => this%var_data(4), &
                uphi => this%var_data(5), &
                chi => this%var_data(7), &
                drdt => this%rhs_data(1), &    
                dphidt => this%rhs_data(2), &    
                dutdt => this%rhs_data(3), &    
                durdt => this%rhs_data(4), &    
                duphidt => this%rhs_data(5), &
                dqmassdt => this%rhs_data(6), &
                dchidt => this%rhs_data(7), &
                En => this%En, &
                Lz => this%Lz, &
                p => this%p, &
                e => this%e, &
                w => this%w, &
                accel => this%accel, &
                udota => this%udota, &
                d2rdt2 => this%d2rdt2 )
      utp = (2.0_wp*mass*ur*ut)/(2.0_wp*mass*r - r**2) + accel(1)
      urp = (-(mass*r**2*ur**2) + (-2*mass + r)**2* &
           (-(r**3*uphi**2) + mass*ut**2))/((2*mass - r)*r**3) + accel(2)
      uphip = (-2.0_wp*uphi*ur)/r + accel(4)

      utinv = 1.0_wp/ut

      drdt = ur*utinv
      dphidt = uphi*utinv
      dutdt = utp*utinv
      durdt = urp*utinv
      duphidt = uphip*utinv
      dqmassdt = -utinv*udota
      ecoschi = e*cos(chi-w)
      dchidt = (p-2.0_wp-2.0_wp*ecoschi)*(1.0_wp+ecoschi)**2 &
              *sqrt((p-6.0_wp-2.0_wp*ecoschi) &
                   /(p**4*mass**2*(p-2.0_wp-2.0_wp*e)*(p-2.0_wp+2.0_wp*e)))

      d2rdt2 = utinv**2*urp - utinv**3*utp*ur
      En = (1.0_wp - 2.0_wp*mass/r)*ut
      Lz = r**2*uphi
    end associate 

  end procedure geod_schw_rhs


! Save the current orbital information in the orbit_base variables.
  module procedure geod_schw_save_globals_1

    use orbit_base
    use parameters, only : mass

    implicit none

    real(wp) :: acosarg, w1, w2
    real(wp), save :: w_old
    real(wp), parameter :: pi = acos(-1.0_wp)

    call this%rhs ( )
    call invert_pe ( this%En, this%Lz, this%p, this%e )
    associate ( w => this%w, &
                chi => this%var_data(7), &
                p => this%p, &
                e => this%e, &
                r => this%var_data(1) )
      w_old = w
      acosarg = (p*mass/r - 1.0_wp)/e
      if (acosarg<-1.0_wp) acosarg = -1.0_wp
      if (acosarg>1.0_wp) acosarg = 1.0_wp
      w1 = chi - acos(acosarg)
      w2 = chi + acos(acosarg) - 2.0_wp*pi
!      print*,'w1 = ', w1
!      print*,'w2 = ', w2
      if (abs(w1-w_old)<=abs(w2-w_old)) then
        w = w1
      else
        w = w2
      end if
!      print*,' acos argument = ', (p*mass/r - 1.0_wp)/e
!      print*,'w = ', w
    end associate
    call orbit_info%set_orbit ( this%var_data(1), this%var_data(2), &
                                this%var_data(4), this%En, this%Lz, &
                                this%var_data(7) )

    call tdc_info%set_tdc ( this%var_data(1), this%rhs_data(1), this%d2rdt2 )


  end procedure geod_schw_save_globals_1


! This procedure is empty as there is no dependence on saved variables.
  module procedure geod_schw_save_globals_2

  end procedure geod_schw_save_globals_2


  module procedure geod_schw_load_globals

    use orbit_base
    use self_force_base, only : sf
    use parameters, only : mass, q_charge, q_mass, use_analytic_driver, &
                           a_timelevels
    use time_info, only : get_current_dtime
    use all_integrators, only : mol_int
    use acceleration_history, only : ah


    implicit none

    real(wp) :: gurr, gutt, guthth, guphiphi, updota

    call sf%get_force ( this%force )

    associate ( r => this%var_data(1), &
                phi => this%var_data(2), &
                ut => this%var_data(3), &
                ur => this%var_data(4), &
                uphi => this%var_data(5), &
                drdt => this%rhs_data(1), &
                dphidt => this%rhs_data(2), &
                dutdt => this%rhs_data(3), &
                durdt => this%rhs_data(4), &
                duphidt => this%rhs_data(5), &
                force => this%force, &
                accel => this%accel, &
                udota => this%udota )
      gurr = 1.0_wp - 2.0_wp*mass/r
      gutt = -1.0_wp/gurr
      guthth = 1.0_wp/r**2
      guphiphi = guthth

      accel(1) = (gutt + ut**2)*force(1) + ut*( ur*force(2) &
                                              + uphi*force(4) )
      accel(2) = (gurr + ur**2)*force(2) + ur*( ut*force(1) &
                                              + uphi*force(4) )
      accel(3) = 0.0_wp
      accel(4) = (guphiphi + uphi**2)*force(4) &
                           + uphi*( ut*force(1) + ur*force(2) )

      accel = q_charge/q_mass*accel

      if (use_analytic_driver) then
        udota = 0.0_wp
      else
        udota = q_charge*dot_product( (/ ut, ur, 0.0_wp, uphi /), force)
      end if

      call sf%set_accel ( accel(1), accel(2), accel(3), accel(4), udota )

      if ( a_timelevels > 0 .and. mol_int%complete_step ( ) ) then
        call ah%cycle_timelevels ( accel, get_current_dtime ( ) )
      end if
    end associate
  end procedure geod_schw_load_globals


  subroutine invert_pe ( En, Lz, p, e )
  !! Routine to convert from Energy and Angular momentum per unit mass to
  !! semi-latus rectum and eccentricity.
  !!
  !! The equations in this routine was obtained using Mathematica by
  !! inverting the relationships
  !! \[
  !!    E^2 = \frac{(p-2-2e)(p-2+2e)}{p(p-3-e^2)},
  !! \]
  !! \[
  !!    L_z^2 = \frac{p^2 M^2}{p-3-e^2}.
  !! \]
  !! Note that this routine only works outside the separatrix and will give
  !! wrong results inside.
    use parameters, only : mass

    implicit none

    real(wp), intent(in) :: En
    !! The energy per unit mass, \(E\).
    real(wp), intent(in) :: Lz
    !! The angular momentum per unit mass, \(L_z\).
    real(wp), intent(out) :: p
    !! On output the semi-latus rectum, \(p\).
    real(wp), intent(out) :: e
    !! On output the eccentricity.
    complex(wp), parameter :: zi = (0.0_wp ,1.0_wp)
    complex(wp), parameter :: onepisqrt3 = 1.0_wp+zi*sqrt(3.0_wp)
    complex(wp), parameter :: onemisqrt3 = 1.0_wp-zi*sqrt(3.0_wp)
    real(wp), parameter :: threesqrt3 = 3.0_wp*sqrt(3.0_wp)
    real(wp), parameter :: third = 1.0_wp/3.0_wp
    complex(wp) :: pc, ec

    real(wp) :: E2, Lz2, M2, M4, M6, E22, Lz22, Lz23, onem3E2, E22Lz2, Lz2M2
    real(wp) :: onemE2, fac, twoLz2p4M2, sixE2M2inv, thirtysixE22M2inv
    complex(wp) :: sqrtfac, sqrtfac2, cuberootfac, cuberootfacinv, mainfac

    M2 = mass**2
    M4 = M2**2
    M6 = M2*M4

    E2 = En**2
    E22 = E2**2
    Lz2 = Lz**2
    Lz22 = Lz2**2
    Lz23 = Lz2*Lz22
    onem3E2 = 1.0_wp - 3.0_wp*E2
    E22Lz2 = E22*Lz2
    Lz2M2 = Lz2*M2
    onemE2 = 1.0_wp - E2

    sqrtfac = sqrt(cmplx(E22Lz2*M6*((8.0_wp-36.0*E2+27.0*E22)*Lz2M2 &
                                +16.0_wp*M4+Lz22*onemE2),0.0_wp,wp))
    sqrtfac2 = 8.0_wp*(-8.0_wp*M6+threesqrt3*sqrtfac)
    cuberootfac = (-Lz23-24.0_wp*(2.0_wp-6.0_wp*E2+9.0_wp*E22)*Lz2*M4 &
                   -12.0_wp*Lz22*M2*onem3E2+sqrtfac2)**third
    cuberootfacinv = 1.0_wp/cuberootfac
    fac = Lz22+16.0_wp*M4+8.0_wp*Lz2M2*onem3E2
    twoLz2p4M2 = 2.0_wp*(Lz2+4.0_wp*M2)
    mainfac = onemisqrt3*cuberootfac+onepisqrt3*cuberootfacinv*fac+twoLz2p4M2
    sixE2M2inv = 1.0_wp/(6.0_wp*E2*M2)
    thirtysixE22M2inv = sixE2M2inv/(6.0_wp*E2)

    pc = mainfac*sixe2M2inv
    ec = sqrt(-3.0_wp+mainfac*sixE2M2inv-mainfac**2*thirtysixE22M2inv/Lz2)

    p = real(pc,wp)
    e = real(ec,wp)
  end subroutine invert_pe


  module procedure geod_schw_output

    use output_base
    use time_info

    implicit none

    character(len=9) :: filename = 'orbit.asc'
    integer(ip) :: io_id
    real(wp) :: time

    io_id = this%io_id 
    if (io_id < 0 ) then
      io_id = next_available_io_id ()
      this%io_id = io_id
      print*,'Opening ', filename, ' with id ', io_id
      open(io_id, file=filename, status='replace', action='write')
      write(io_id,*) '#1: time, 2: r, 3: phi, 4: ut, 5: ur, 6: uphi, 7: qmass, 8: chi, 9: drdt, 10: dphidt, 11: dutdt, 12: durdt, 13: duphidt, 14: dqmassdt, 15: dchidt,  16: d2rdt2, 17: En, 18: Lz, 19: p, 20: e, 21: w'
    end if

    time = get_current_time ( )
    write(io_id,'(*(es23.15e3,1x))') time, this%var_data, this%rhs_data, &
                                    this%d2rdt2, this%En, this%Lz, this%p, &
                                    this%e, this%w
  end procedure geod_schw_output


  module procedure close_geod_schw

    deallocate ( this%var_data, this%rhs_data, this%tmp_data )
    if ( this%io_id>=0 ) then
      close (this%io_id)
    end if

  end procedure close_geod_schw

end submodule geodesic_schwarzschild_implementation
