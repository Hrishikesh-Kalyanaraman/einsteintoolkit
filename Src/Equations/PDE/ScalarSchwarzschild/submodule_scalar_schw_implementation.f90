submodule(scalar_schw) scalar_schw_implementation
!! The implementation of the interfaces defined in [[scalar_schw]].

    implicit none

contains

  module procedure scal_schw_init

    use grid
    use parameters, only : n_elems, order, mass, lmin, lmax, r_center, &
                           gaussian_center, sigma, amplitude, use_particle, &
                           use_generic_orbit, use_exact_initial_data, &
                           use_filter, filter_order, filter_cutoff
    use numerics, only : transition, invert_tortoise
    use all_integrators, only : mol_ntmp
    use acceleration_history, only : ah

    implicit none
 
    integer(ip) :: i, j, k
    integer(ip), parameter :: ncoeffs = 4
    character(len=3), dimension(3) :: var_names = (/ 'psi', 'rho', 'pi ' /)
    character(len=7), dimension(3) :: rhs_names = (/ 'psi_rhs', 'rho_rhs', &
                                                     'pi_rhs ' /)
    character(len=7), dimension(3) :: tmp_names = (/ 'psi_tmp', 'rho_tmp', &
                                                     'pi_tmp ' /)
    character(len=8), dimension(2) :: flux_names = (/ 'rho_flux', 'pi_flux ' /)
    character(len=3), dimension(4) :: coeff_names = (/ 'crr', 'ctr', &
                                                       'cr ', 'ct ' /)
    character(len=12), dimension(2) :: lambda_names = (/ 'lambda_minus', &
                                                        'lambda_plus ' /)
    character(len=3), dimension(2,2) :: s_names = reshape ( &
                                        (/ 'S11', 'S12', 'S21', 'S22' /), &
                                        (/ 2, 2 /) )
    character(len=7), dimension(2,2) :: sinv_names = reshape ( &
                                        (/ 'S11_inv', 'S12_inv', &
                                           'S21_inv', 'S22_inv' /), &
                                        (/ 2, 2 /) )
    real(wp) :: fT, fTp, fTpp, Omega, Omegap, eL, eLp, H, Hp, rm2M, &
                term1, term2, C02, C11, sqrtfac, sqrtfacinv
    integer(ip) :: lmode, nmodes, nvars
    integer :: allocation_status

    print*,'scalar schwarzschild init called'

! Set the number of variables per mode.
    this%nvars = 3
! Set the equation name.
    if ( .not. allocated(this%ename) ) then
      allocate ( character(18) :: this%ename, stat=allocation_status )
      if ( allocation_status > 0 ) then
        print*,'Error while allocating scalar_schw equation name'
        stop
      end if
    end if
    this%ename = 'Schwarzschild_wave'

! Setup mode info.
    nmodes = nmodes_of_l ( lmin, lmax )
    this%nmodes = nmodes
    allocate ( this%ll(nmodes), stat=allocation_status )
    if (allocation_status>0) then
      print*,'Allocation of ll failed'
      stop
    end if
    allocate ( this%mm(nmodes), stat=allocation_status )
    if (allocation_status>0) then
      print*,'Allocation of mm failed'
      stop
    end if
    call set_lm_mode_info ( this, lmin, lmax )

! Setup grid functions for the evolution data and it's rhs.
    nvars = this%nvars
    allocate ( this%eq_data(nvars,nmodes) )
    allocate ( this%eq_rhs_data(nvars,nmodes) )
    do i = 1, nvars
      do j = 1, nmodes
        this%eq_data(i,j) = cgf ( n_elems, order, &
                                  convert_var_name(var_names(i),j) )
        this%eq_rhs_data(i,j) = cgf ( n_elems, order, &
                                      convert_var_name(rhs_names(i),j) )
      end do
    end do

! Setup up temporary storage required by the time integrator.
    this%ntmp = mol_ntmp ( )
    print*,'Allocating ', this%ntmp, ' grid functions for temporary storage'
    allocate ( this%eq_tmp_data(nvars,nmodes,this%ntmp) )
    do k = 1, this%ntmp
      do i = 1, nvars
        do j = 1, nmodes
          this%eq_tmp_data(i,j,k) = cgf ( n_elems, order, &
                                        convert_var_name( &
                                        convert_var_name(tmp_names(i),j), k) )
        end do
      end do
    end do

! Setup grid functions for the flux data.
    allocate ( this%eq_flux_data(nvars,nmodes) )
    do i = 1, 2
      do j = 1, nmodes
        this%eq_flux_data(i,j) = cgf ( n_elems, order, &
                                       convert_var_name(flux_names(i),j) )
      end do
    end do

! Setup grid functions for the equation coefficients.
    allocate ( this%eq_coeffs(ncoeffs) )
    do i = 1, ncoeffs
      this%eq_coeffs(i) = rgf ( n_elems, order, coeff_names(i) )
    end do
    allocate ( this%eq_lcoeffs(nmodes) )
    do i = 1, nmodes
      this%eq_lcoeffs(i) = rgf ( n_elems, order, convert_var_name('c',i) )
    end do

! Setup boundary grid functions for the characteristic data.
    allocate ( this%eq_lambda(2), this%eq_s(2,2), this%eq_sinv(2,2) )
    do i = 1, 2
      this%eq_lambda(i) = rgfb ( n_elems, lambda_names(i) )
      do j = 1, 2
        this%eq_s(i,j) = rgfb ( n_elems, s_names(i,j) )
        this%eq_sinv(i,j) = rgfb ( n_elems, sinv_names(i,j) )
      end do
    end do

! Setup grid functions for various radial coordinates.
    this%r_schw = rgf ( n_elems, order, 'r_schw' )
    this%r_star = rgf ( n_elems, order, 'r_star' )

! Setup arrays for storage of derivatives. This is done in this way in order
! to allow for OpenMP parallelization.
    allocate ( drho(order+1,nmodes), dpi(order+1,nmodes) )
    if (.not. allocated(flux_result) ) then
      allocate ( flux_result(order+1,2,nmodes) )
    end if

! Setup data pointers as required by the MoL integrator.
    allocate ( this%data_pointer(nvars,nmodes,-1:this%ntmp) )
    do i = 1, nvars
      do j = 1, nmodes
        this%data_pointer(i,j,-1)%p => this%eq_rhs_data(i,j)
        this%data_pointer(i,j,0)%p => this%eq_data(i,j)
        do k = 1, this%ntmp
          this%data_pointer(i,j,k)%p => this%eq_tmp_data(i,j,k)
        end do
      end do
    end do

! Initialize the reference element.
    if (use_filter) then
      this%refelem = ref_element( order, filter_order, filter_cutoff )
    else
      this%refelem = ref_element( order )
    end if

! Compute coefficients in the inner hyberloidal region.
    do i = 1, Tminus_ind(1)
      do j = 1, order+1
        if ( (i==1) .and. (j==1) ) then
          this%r_star%elems(i)%var(j) = maxval ( empty )
          this%r_schw%elems(i)%var(j) = 2.0_wp*mass
          this%eq_coeffs(1)%elems(i)%var(j) = 0.0_wp
          this%eq_coeffs(2)%elems(i)%var(j) = 1.0_wp
          this%eq_coeffs(3)%elems(i)%var(j) = 0.0_wp
          this%eq_coeffs(4)%elems(i)%var(j) = 0.0_wp
          do k = 1, this%nmodes
            this%eq_lcoeffs(k)%elems(i)%var(j) = 0.0_wp
          end do
        else
          associate ( rho => rho%elems(i)%var(j), &
                      rstar => this%r_star%elems(i)%var(j), &
                      rschw => this%r_schw%elems(i)%var(j) )
            call transition ( rho, Tminus, rho_min, fT, fTp, fTpp )
            Omega = 1.0_wp - rho / rho_min * fT
            Omegap = -(fT+rho*fTp)/rho_min
            eL = 1.0_wp + rho**2 * fTp / rho_min
            eLp = rho*(2.0_wp*fTp+rho*fTpp) / rho_min
            H = -1.0_wp + Omega**2/eL
            Hp = (2.0_wp*Omega*Omegap*eL-Omega**2*eLp)/eL**2
            rstar = rho/Omega
            rm2M = invert_tortoise ( rstar, mass )
            rschw = 2.0_wp * mass + rm2M
            this%eq_coeffs(1)%elems(i)%var(j) = (1.0_wp+H)/(1.0_wp-H)
            this%eq_coeffs(2)%elems(i)%var(j) = -2.0_wp*H/(1.0_wp-H)
            this%eq_coeffs(3)%elems(i)%var(j) = Hp/(1.0_wp-H)
            this%eq_coeffs(4)%elems(i)%var(j) = -Hp/(1.0_wp-H)
            term1 = rm2M/rschw**3
            term2 = 2.0_wp*mass/rschw
            do k = 1, this%nmodes
              lmode = this%ll(k)
              this%eq_lcoeffs(k)%elems(i)%var(j) = -1.0_wp/(1.0_wp-H**2)*term1 &
                                                   *( lmode*(lmode+1) + term2 )
            end do
          end associate
        end if
      end do
    end do 

! Compute coefficients in the source region.
    do i = Tminus_ind(2), Tplus_ind(1)
      do j = 1, order+1
        associate ( rho => rho%elems(i)%var(j), &
                    rstar => this%r_star%elems(i)%var(j), &
                    rschw => this%r_schw%elems(i)%var(j) )
          rstar = rho
          rm2M = invert_tortoise ( rstar, mass )
          rschw = 2.0_wp * mass + rm2M
          this%eq_coeffs(1)%elems(i)%var(j) = 1.0_wp
          this%eq_coeffs(2)%elems(i)%var(j) = 0.0_wp
          this%eq_coeffs(3)%elems(i)%var(j) = 0.0_wp
          this%eq_coeffs(4)%elems(i)%var(j) = 0.0_wp
          term1 = rm2M/rschw**3
          term2 = 2.0_wp*mass/rschw
          do k = 1, this%nmodes
            lmode = this%ll(k)
            this%eq_lcoeffs(k)%elems(i)%var(j) = -term1 &
                                                 *( lmode*(lmode+1) + term2 )
          end do

        end associate
      end do
    end do 

! Compute coefficients in the outer hyperboloidal region.
    do i = Tplus_ind(2), n_elems
      do j = 1, order+1
        if ( (i==n_elems) .and. (j==order+1) ) then
          this%r_star%elems(i)%var(j) = -maxval ( empty )
          this%r_schw%elems(i)%var(j) = -maxval ( empty )
          this%eq_coeffs(1)%elems(i)%var(j) = 0.0_wp
          this%eq_coeffs(2)%elems(i)%var(j) = -1.0_wp
          this%eq_coeffs(3)%elems(i)%var(j) = 0.0_wp
          this%eq_coeffs(4)%elems(i)%var(j) = 0.0_wp
          do k = 1, this%nmodes
            lmode = this%ll(k)
            this%eq_lcoeffs(k)%elems(i)%var(j) = -lmode*(lmode+1) &
                                                /(2.0_wp*rho_max**2)
          end do
        else
          associate ( rho => rho%elems(i)%var(j), &
                      rstar => this%r_star%elems(i)%var(j), &
                      rschw => this%r_schw%elems(i)%var(j) )
            call transition ( rho, Tplus, rho_max, fT, fTp, fTpp )
            Omega = 1.0_wp - rho / rho_max * fT
            Omegap = -(fT+rho*fTp)/rho_max
            eL = 1.0_wp + rho**2 * fTp / rho_max
            eLp = rho*(2.0_wp*fTp+rho*fTpp) / rho_max
            H = 1.0_wp - Omega**2/eL
            Hp = -(2.0_wp*Omega*Omegap*eL-Omega**2*eLp)/eL**2
            rstar = rho/Omega
            rm2M = invert_tortoise ( rstar, mass )
            rschw = 2.0_wp * mass + rm2M
            this%eq_coeffs(1)%elems(i)%var(j) = (1.0_wp-H)/(1.0_wp+H)
            this%eq_coeffs(2)%elems(i)%var(j) = -2.0_wp*H/(1.0_wp+H)
            this%eq_coeffs(3)%elems(i)%var(j) = -Hp/(1.0_wp+H)
            this%eq_coeffs(4)%elems(i)%var(j) = -Hp/(1.0_wp+H)
            term1 = rm2M/rschw**3
            term2 = 2.0_wp*mass/rschw
            do k = 1, this%nmodes
              lmode = this%ll(k)
              this%eq_lcoeffs(k)%elems(i)%var(j) = -1.0_wp/(1.0_wp-H**2)*term1 &
                                                   *( lmode*(lmode+1) + term2 )
            end do
          end associate
          
        end if
      end do
    end do 

! Store information at element boundaries for characteristic fluxes.
    do i = 1, n_elems
      do j = 1, 2
        if (j==1) then
          C02 = this%eq_coeffs(1)%elems(i)%var(1)
          C11 = this%eq_coeffs(2)%elems(i)%var(1)
        else
          C02 = this%eq_coeffs(1)%elems(i)%var(order+1)
          C11 = this%eq_coeffs(2)%elems(i)%var(order+1)
        end if
        sqrtfac = sqrt(4.0_wp*C02+C11**2)
        sqrtfacinv = 1.0_wp/sqrt(4.0_wp*C02+C11**2)
        this%eq_lambda(1)%elems(i)%bvar(j) = 0.5_wp*(-C11-sqrtfac)
        this%eq_lambda(2)%elems(i)%bvar(j) = 0.5_wp*(-C11+sqrtfac)
        this%eq_s(1,1)%elems(i)%bvar(j) = 0.5_wp*(C11+sqrtfac)
        this%eq_s(1,2)%elems(i)%bvar(j) = 0.5_wp*(C11-sqrtfac)
        this%eq_s(2,1)%elems(i)%bvar(j) = 1.0_wp
        this%eq_s(2,2)%elems(i)%bvar(j) = 1.0_wp
        this%eq_sinv(1,1)%elems(i)%bvar(j) = sqrtfacinv
        this%eq_sinv(1,2)%elems(i)%bvar(j) = 0.5_wp*(1.0_wp-C11*sqrtfacinv)
        this%eq_sinv(2,1)%elems(i)%bvar(j) = -sqrtfacinv
        this%eq_sinv(2,2)%elems(i)%bvar(j) = 0.5_wp*(1.0_wp+C11*sqrtfacinv)
      end do

! When not using a particle set up a gaussian initial data profile.
      if (.not. use_particle) then
        do k = 1, this%nmodes
          do j = 1, rho%elems(i)%order+1
            associate ( rho => rho%elems(i)%var(j), &
                        dpsidt => this%eq_data(2,k)%elems(i)%var(j) )
              dpsidt = amplitude * exp ( -0.5_wp*((rho-gaussian_center) &
                                         /sigma)**2 )
            end associate
          end do
        end do
      end if
    end do

! When using a particle initialize the effective source.
    if (use_particle) then
      print*,'nmodes = ', this%nmodes
      print*,'nvars = ', this%nvars
      print*,'Initializing effective source'
      call this%effs%init ( this%nmodes, this%nvars, this%ll, this%mm, mass )
    end if

! When using a generic orbit initialize the time dependent coordinate
! transformation. This requires the orbit evolution to have already been
! initialized
    if (use_generic_orbit ) then
      call this%time_dep_coord%init ( )
      call this%time_dep_coord%set_coefficients ( this%eq_coeffs, &
                                                  this%eq_lcoeffs, &
                                                  this%eq_lambda, &
                                                  this%eq_s, this%eq_sinv, &
                                                  rho, this%r_star, &
                                                  this%r_schw, this%ll )
    end if

    print*,'Scal_Schw initialized'
  end procedure scal_schw_init


  module procedure scal_schw_rhs

    use grid, only : rho, drdrho
    use parameters, only : n_elems, use_particle, tsigma, torder, &
                           exact_initial_data_lmax, use_generic_orbit, &
                           q_charge, order, turn_on_source_smoothly, &
                           a_timelevels, use_adot, use_addot
    use numerics, only : time_window
    use time_info, only : get_current_time, get_current_sub_dtime
    use world_tube, only : wt
    use self_force_base, only : sf
!    use all_integrators, only : mol_int
    use acceleration_history, only : ah

    implicit none

    integer(ip) :: i, j, k, nmodes
    real(wp) :: r, phi, ur, En, Lz
    real(wp) :: tfac, dtfac_dt, dt2fac_dt2
    real(wp), dimension(4,3) :: accel

! If we are using a particle.
    if (use_particle) then
! Get the current orbit information.
      call orbit_info%get_orbit ( r, phi, ur, En, Lz )
! Calculate the time window values.
      call time_window ( get_current_time (), tsigma, torder, tfac, &
                         dtfac_dt, dt2fac_dt2 )
! Tell the effective source about the time window.
      if ( turn_on_source_smoothly ) then
!        print*,'Setting time window: ', get_current_time ( ), tfac, dtfac_dt, dt2fac_dt2
        call this%effs%set_time_window ( tfac, dtfac_dt, dt2fac_dt2, -1 )
      else
        call this%effs%set_time_window ( tfac, dtfac_dt, dt2fac_dt2, &
                                         exact_initial_data_lmax )
      end if

      if ( a_timelevels < 0 ) then
! Get the current acceleration.
        call sf%get_accel ( accel(:,1) )
        call sf%get_daccel_dt ( accel(:,2) )
        call sf%get_d2accel_dt2 ( accel(:,3) )
      else
        accel = ah%extrapolate ( get_current_sub_dtime ( ) )
      end if

      if ( .not. (use_adot .or. use_addot) ) then
        accel(:,2) = 0.0_wp
      end if
      if ( .not. use_addot ) then
        accel(:,3) = 0.0_wp
      end if

! Set the particle position information. We only need the r and phi
! components of the acceleration and it's derivatives.
      call this%effs%set_particle_pos ( r, phi, ur, En, Lz, accel(2,1), &
                                        accel(4,1), accel(2,2), accel(4,2), &
                                        accel(2,3), accel(4,3) )
!      print*,'time = ', get_current_time ()
!      print*,'set_paticle_pos arguments = ', get_current_time (), &
!                                        r, phi, ur, En, Lz, accel(2,1), &
!                                        accel(4,1), accel(2,2), accel(4,2), &
!                                        accel(2,3), accel(4,3)

! Evaluate the source at the current time.
      if ( wt%wsize > 0 ) then
        call this%effs%evaluate_source ( this%r_schw, wt )
      end if
    end if

    nmodes = this%nmodes

! Calculate the current fluxes at all element boundaries.
    call this%flux()

! Calculate the right hand sides for all modes and all elements.
!$OMP PARALLEL DO private(k,i) shared(this,drdrho,drho,dpi,nmodes)
    do k = 1, nmodes
      do i = 1, n_elems
! Calculate the radial derivatives (by multiplying with the derivative matrix
! and the Jacobian) and store them. Storing them this way is done to avoid
! segmentation faults when OpenMP parallelizing.
        drho(:,k) = drdrho%elems(i)%var * &
                    matmul ( this%refelem%dr, &
                             this%eq_data(2,k)%elems(i)%var )
        dpi(:,k) = drdrho%elems(i)%var * &
                    matmul ( this%refelem%dr, &
                             this%eq_data(3,k)%elems(i)%var )
        do j = 1, order+1
          associate ( psi => this%eq_data(1,k)%elems(i)%var(j), &
                      rho => this%eq_data(2,k)%elems(i)%var(j), &
                      pi => this%eq_data(3,k)%elems(i)%var(j), &
                      rho_flux => this%eq_flux_data(1,k)%elems(i)%var(j), &
                      pi_flux => this%eq_flux_data(2,k)%elems(i)%var(j), &
                      dpsidt => this%eq_rhs_data(1,k)%elems(i)%var(j), &
                      drhodt => this%eq_rhs_data(2,k)%elems(i)%var(j), &
                      dpidt => this%eq_rhs_data(3,k)%elems(i)%var(j), &
                      crr => this%eq_coeffs(1)%elems(i)%var(j), &
                      crt => this%eq_coeffs(2)%elems(i)%var(j), &
                      cr => this%eq_coeffs(3)%elems(i)%var(j), &
                      ct => this%eq_coeffs(4)%elems(i)%var(j), &
                      c => this%eq_lcoeffs(k)%elems(i)%var(j) )
!                      effs => this%effs%source(1,k)%elems(i)%var(j) )

            dpsidt = rho
            drhodt = crr*dpi(j,k) + crt*drho(j,k) + cr*pi + ct*rho &
                    + c*psi + rho_flux
            dpidt = drho(j,k) + pi_flux
            if ( use_particle .and. wt%wsize > 0 ) then
! Add the effective source if necessary.
              drhodt = drhodt + q_charge*this%effs%source(1,k)%elems(i)%var(j)
            end if
          end associate
        end do
      end do
    end do
!$OMP END PARALLEL DO

  end procedure scal_schw_rhs


! This procedure is empty as extracting the self-force requires updated
! values for the orbit.
  module procedure scal_schw_save_globals_1

  end procedure scal_schw_save_globals_1


! Extract the self-force (using the self-force observer) and store the values
! in self_force_base variables to be used by the orbit evolution.
  module procedure scal_schw_save_globals_2

    use self_force_observer
    use self_force_base
    use parameters, only : lmax, use_particle, fit_high_l
    use time_info, only : get_current_time
    use world_tube, only : wt

    implicit none

    real(wp) :: ft, fr, ftheta, fphi

!     print*,'Scal_Schw save globals called at time: ', get_current_time()
    if (use_particle) then
      call self_force_extract ( this%effs, wt%wsize )

      if (fit_high_l) then
! Store the already summed and extrapolated values for the self-force.
        ft = sfobs%ftl(lmax+2)
        fr = sfobs%frl(lmax+2)
        ftheta = 0.0_wp
        fphi = sfobs%fphil(lmax+2)
      else
! Sum over m and store in the self_force_base class.
        ft = sum(sfobs%ftl)
        fr = sum(sfobs%frl)
        ftheta = 0.0_wp
        fphi = sum(sfobs%fphil)
      end if

      call sf%set_force ( ft, fr, ftheta, fphi )
    end if
    
  end procedure scal_schw_save_globals_2


! This routine was introduced in order to be able to point to the effective
! source object.
  subroutine self_force_extract ( effs, wsize )

    use self_force_observer, only : sfobs

    implicit none

    class(eff_source), target, intent(in) :: effs
    integer(ip), intent(in) :: wsize
    class(eff_source), pointer :: effsp

    if ( wsize > 0 ) then
      call sfobs%extract ( )
    else
      effsp => effs
      call sfobs%extract ( effsp )
    end if

  end subroutine self_force_extract


! Based on the current orbit, set the time dependent evolution equation
! coefficients and adjust the timestep.
  module procedure scal_schw_load_globals

    use grid, only : rho, delta_rho_min
    use parameters, only : use_generic_orbit, coufac, tsigma, t_initial
    use time_info, only : set_dtime, short_timesteps_active, get_current_time

    implicit none

!    print*,'Scal_Schw load globals called at time: ', get_current_time()
    if (use_generic_orbit ) then
      call this%time_dep_coord%set_coefficients ( this%eq_coeffs, &
                                                  this%eq_lcoeffs, &
                                                  this%eq_lambda, &
                                                  this%eq_s, this%eq_sinv, &
                                                  rho, this%r_star, &
                                                  this%r_schw, this%ll )
      if (short_timesteps_active) then
        if ( get_current_time() >=4.0_wp*tsigma+t_initial ) then
          print*,'Source turned on, switching back to normal timesteps'
          call set_dtime ( coufac*delta_rho_min/this%time_dep_coord%maxspeed )
          short_timesteps_active = .false.
        end if
      else
        call set_dtime ( coufac*delta_rho_min/this%time_dep_coord%maxspeed )
      end if
    end if

  end procedure scal_schw_load_globals


  module procedure scal_schw_apply_filter

    use parameters, only : n_elems

    implicit none

    integer(ip) :: i, j, k, nmodes, nvars

    nmodes = this%nmodes
    nvars = this%nvars

    if (this%refelem%have_filter) then
!$OMP PARALLEL DO private(i,j,k) shared(this,nmodes,nvars,n_elems)
      do k = 1, nmodes
        do j = 1, nvars 
          do i = 1, n_elems
            this%eq_data(j,k)%elems(i)%var = matmul ( this%refelem%filter, &
                                               this%eq_data(j,k)%elems(i)%var )
          end do
        end do
      end do
!$OMP END PARALLEL DO
    end if
  end procedure scal_schw_apply_filter


! Calcuate the characteristic fluxes. This has all been moved to a subroutine
! in order to avoid segmentation fault issues with OpenMP parallelization.
  module procedure scal_schw_flux

    implicit none

    integer(ip) :: k, nmodes

    nmodes = this%nmodes

!$OMP PARALLEL DO private(k) shared(nmodes,this)
    do k = 1, nmodes
      call get_elem_flux ( this%eq_data(2,k), this%eq_data(3,k), nmodes, k, &
                           this%r_schw, this%effs, this%time_dep_coord, &
                           this%eq_lambda, this%eq_s, this%eq_sinv, &
                           this%eq_coeffs, this%refelem, &
                           this%eq_flux_data(1,k), this%eq_flux_data(2,k) )
    end do
!$OMP END PARALLEL DO

  end procedure scal_schw_flux


! Read in frequency domain data in Schwarzschild time and tortoise coordinates
! and transform to the computational coordinates.
  module procedure read_all_modes

    use parameters, only : exact_initial_data_lmax, q_charge, mass, &
                           use_world_tube, use_generic_orbit
    use output_base
    use grid, only : rho, Tminus_ind, Tplus_ind
    use world_tube, only : wt

    implicit none

    integer(ip) :: i, j, jhigh, k, l, n, np, io_id=-1
    character(len=1024) :: filein
    real(wp) :: rho_r, Psire, Psiim, drstarPsire, drstarPsiim, dtPsire, dtPsiim
    real(wp) :: radius
    integer(ip) :: sgn
    complex(wp), dimension(1) :: ss, sst, ssr
    real(wp), parameter :: eps = 1.0e-12_wp

! All grid functions have the same number of elements, so we can get that value
! from any of them.
    n = this%r_schw%n

! Loop over modes.
    do k = 1, this%nmodes
! Only read in data when required
      if ( this%ll(k) <= exact_initial_data_lmax ) then
        filein = construct_filename ( this%ll(k), this%mm(k) )
        print*,'Overwriting mode ', k,' with analytic data'
! Use the features of output_base to get the next file unit number.
        if (io_id < 0) then
          io_id = next_available_io_id ( )
        end if
! Open the file.
        open(io_id, file = filein, status='old', action='read')

! Loop over all elements and nodes, read in retarded data, scale it with the
! charge and r_schw
        do i = 1, n
          np = rho%elems(i)%order+1
          jhigh = np
          if (i==n) jhigh = np - 1
          do j = 1, jhigh
            read(io_id,*) rho_r, Psire, Psiim, drstarPsire, drstarPsiim, &
                          dtPsire, dtPsiim
            if ( abs(rho%elems(i)%var(j)-rho_r)>1e-12) then
              print*,'Coordinates of initial data inconsistent with grid for (i,j) = ', i, j
              print*,'rho%elems(i)%var(j) = ', rho%elems(i)%var(j)
              print*,'rho = ', rho_r
              stop
            end if

            associate ( psi => this%eq_data(1,k)%elems(i)%var(j), &
                        dpsidt => this%eq_data(2,k)%elems(i)%var(j), &
                        dpsidr => this%eq_data(3,k)%elems(i)%var(j), &
                        rschw => this%r_schw%elems(i)%var(j) )
              psi = q_charge*rschw*cmplx(Psire,Psiim,wp)
              dpsidt = q_charge*rschw*cmplx(dtPsire,dtPsiim,wp)
! We are evolving r_schw*Psi but need the radial derivative with respect to
! r_star.
              dpsidr = q_charge*(rschw*cmplx(drstarPsire,drstarPsiim,wp) &
                      +(1.0_wp-2.0_wp*mass/rschw)*cmplx(Psire, Psiim, wp))
            end associate
          end do
        end do

! When using a world tube, we have to subtract the singular field inside it.
        if (use_world_tube .and. wt%wsize>0 ) then
          do i = wt%windex1, wt%windex2
            np = rho%elems(i)%order+1
            ! The following assumes the particle is in the middle of the
            ! domain.
            if ( i<=n/2 ) then
              ! If we are to the left of the particle sgn=-1
              sgn = -1
            else 
              ! Otherwise sgn = 1
              sgn = 1
            end if
            do j = 1, np
              radius = this%r_schw%elems(i)%var(j)
              call this%effs%get_singular ( radius, k, ss )
              call this%effs%get_dsingular_dt ( radius, sgn, k, sst )
              call this%effs%get_dsingular_dr ( radius, sgn, k, ssr )
              associate ( psi => this%eq_data(1,k)%elems(i)%var(j), &
                          dpsidt => this%eq_data(2,k)%elems(i)%var(j), &
                          dpsidr => this%eq_data(3,k)%elems(i)%var(j) )
! Scale with the charge before subtracting.
                psi = psi - q_charge*ss(1)
                dpsidt = dpsidt - q_charge*sst(1)
                dpsidr = dpsidr - q_charge*ssr(1)
              end associate
            end do
          end do
        end if

! Transform to hyperboloidal coordinates in the inner hyperboloidal region.
        do i = 1, Tminus_ind(1)
          np = rho%elems(i)%order+1
          do j = 1, np
            associate ( dpsidt => this%eq_data(2,k)%elems(i)%var(j), &
                        dpsidr => this%eq_data(3,k)%elems(i)%var(j) )
              call this%tortoise_to_hyperboloidal ( i, j, .true., &
                                                    dpsidt, dpsidr )
            end associate
          end do
        end do

! For generic orbits transform to time dependent coordinates in the source
! region.
        if ( use_generic_orbit ) then
          do i = Tminus_ind(2), Tplus_ind(1)
            np = rho%elems(i)%order+1
            do j = 1, np
              associate ( dpsidt => this%eq_data(2,k)%elems(i)%var(j), &
                          dpsidr => this%eq_data(3,k)%elems(i)%var(j) )
                call this%time_dep_coord%tortoise_to_tdc ( i, j, &
                                                           dpsidt, dpsidr )
              end associate
            end do
          end do
        end if

! Transform to hyperboloidal coordinates in the outer hyperboloidal region.
        do i = Tplus_ind(2), n
          np = rho%elems(i)%order+1
          do j = 1, np
            associate ( dpsidt => this%eq_data(2,k)%elems(i)%var(j), &
                        dpsidr => this%eq_data(3,k)%elems(i)%var(j) )
              call this%tortoise_to_hyperboloidal ( i, j, .false., &
                                                    dpsidt, dpsidr )
            end associate
          end do
        end do

! Fill in the data at the horizon and scri+ by linear extrapolation.
        np = rho%elems(n)%order+1
        do i = 1, 3
          this%eq_data(i,k)%elems(1)%var(1) = &
            linear_extrapolate ( rho%elems(1)%var(1:3), &
                                 this%eq_data(i,k)%elems(1)%var(2:3) )
          this%eq_data(i,k)%elems(n)%var(np) = &
            linear_extrapolate ( rho%elems(n)%var(np:np-2:-1), &
                                 this%eq_data(i,k)%elems(n)%var(np-1:np-2:-1) )
        end do
! Close this input file and release the file unit.
        call release_io_id ( io_id )
      end if
    end do
  end procedure read_all_modes


! Output of the computational coordinates to be used by the frequency domain
! code to calculate correct retarded data at the gridpoints used.
  module procedure output_coords

    use parameters, only : mass
    use grid, only : rho, Tminus_ind, Tplus_ind
    use output_base
    use time_info, only : get_current_time
    use numerics, only : invert_tortoise

    implicit none

    integer(ip) :: i, j, io_id=-1
    character(len=10) :: filename = 'coords.asc'
    real(wp) :: diff

! Use output_base to get the next file unit.
    if (io_id < 0 ) then 
      io_id = next_available_io_id ( )
    end if
! Open the file and write the header.
    print*,'Opening ', filename, ' with id ', io_id
    open(io_id, file=filename, status='replace', action='write')
    write(io_id,*)
    write(io_id,*)
    write(io_id,*) '#time = ', get_current_time ( )

! Write the radial coordinate data.
    do i = 1, rho%n
      do j = 1, rho%elems(i)%order+1
        associate ( rho => rho%elems(i)%var(j), &
                    rstar => this%r_star%elems(i)%var(j) )
! In the inner hyperboloidal region we need the difference between the
! computational coordinate and rstar.
          if ( i<=Tminus_ind(1) ) then
            diff = rho - rstar
! In the source region the computational coordinate is rstar.
          else if ( Tminus_ind(2) <= i .and. i<=Tplus_ind(1) ) then
            diff = 0.0_wp
! In the outer hyperboloidal region we need the difference between rstar
! and the computational coordinate.
          else if ( i>=Tplus_ind(2) ) then
            diff = rstar - rho
          end if
! Write the computational coordinate, rstar, r_schw-2M and the difference.
          write(io_id,'(*(es23.15e3,1x))') rho, rstar, &
                                           invert_tortoise ( rstar, mass ), &
                                           diff
        end associate
      end do
    end do

! Close the file and release the file unit.
    call release_io_id ( io_id )

  end procedure output_coords


! Transform from tortoise coordinates to hyperboloidal coordinates.
  module procedure tortoise_to_hyperboloidal

    implicit none

    associate ( c1 => this%eq_coeffs(1)%elems(elem)%var(node) )
      if ( inner ) then
! H = (c1-1)/(c1+1)
! dpsi/dtau = dpsi/dt, dpsi/drho = H/(1+H) dpsi/dt + 1/(1+H) dpsi/drstar
! H/(1+H) = (c1-1)/(2 c1), 1/(1+H) = (c1+1)/(2 c1)
       dpsidr = ( (c1-1.0_wp)*dpsidt + (c1+1.0_wp)*dpsidr )/(2.0_wp*c1)
      else
! H = (1-c1)/(1+c1)
! dpsi/dtau = dpsi/dt, dpsi/drho = H/(1-H) dpsi/dt + 1/(1-H) dpsi/drstar
! H/(1-H) = (1-c1)/(2 c1), 1/(1-H) = (1+c1)/(2 c1)
       dpsidr = ( (1.0_wp-c1)*dpsidt + (1.0_wp+c1)*dpsidr )/(2.0_wp*c1)
      end if
    end associate
  end procedure tortoise_to_hyperboloidal


  function count_digits ( l )
  !! Utility routine to count the number of digits in an integer. Used to
  !! construct filenames.

    implicit none

    integer(ip), intent(in) :: l
    !! The integer to count the number of digits in.
    integer(ip) :: count_digits
    integer(ip) :: ltmp

    count_digits = 1
    ltmp = l/10

    count_the_digits: do
      if (ltmp<1) exit count_the_digits
      count_digits = count_digits+1
      ltmp = ltmp/10
    end do count_the_digits
  end function count_digits


  function construct_filename ( l, m )
  !! Construct unique filenames using a base name and the mode numbers l and m
  !! for reading in external inital data. The base name is constructed based
  !! on the run time parameters [[parameters:input_directory]] and
  !! [[parameters:input_basename]].

    use parameters, only : input_directory, input_basename

    implicit none

    integer(ip), intent(in) :: l
    !! The \(\ell\)-mode of this mode.
    integer(ip), intent(in) :: m
    !! The \(m\)-mode of this mode.
    character(len=1024) :: construct_filename
    !! The return value is the concatenation of the
    !! [[parameters:input_directory]] and [[parameters:input_basename]] and the
    !! mode information.
    character(len=6) :: lstr, mstr
    character(len=4), dimension(6) :: form = (/ "(i1)", "(i2)", "(i3)", &
                                                "(i4)", "(i5)", "(i6)" /)

    write(lstr,form(count_digits(l))) l
    write(mstr,form(count_digits(m))) m

    construct_filename = trim(input_directory)//'/'//trim(input_basename) &
                         //'_l'//trim(lstr)//'m'//trim(mstr)//'.dat'
  end function construct_filename


  function flux ( u, C11, C02 )
  !! Calculate a flux from the derivative data and the equation coefficient.

    use kinds

    implicit none

    complex(wp), dimension(2), intent(in) :: u
    !! The input state vector \((\Upsilon,\Pi)\).
    real(wp), intent(in) :: C11
    !! The coefficient of the second radial derivative term, \(c_{rr}\).
    real(wp), intent(in) :: C02
    !! The coefficient of the mixed second derivative term, \(c_{tr}\).
    complex(wp), dimension(2) :: flux
    !! On return the flux.
    
    flux(1) = -C11 * u(1) - C02 * u(2)
    flux(2) = -u(1)

    return

  end function flux

! Add a mode number to a variable name. Used for constructing unique filenames.
  module procedure convert_var_name

    implicit none

    character(len=4) :: tmp
    character(len=7), dimension(4) :: form = (/ "(a1,i1)", "(a1,i2)", &
                                                "(a1,i3)", "(a1,i4)" /)
    integer(ip) :: ndigits

! Should probably use count_digits from above as it does not involve floating
! point arithmetic.
    ndigits = floor(log(real(mode,wp))/log(10.0_wp)+1.0_wp)
    if (ndigits>4) then
      print*,'convert_var_name called with a mode number with more than 4 digits'
      stop
    end if
    write(tmp,form(ndigits)) '.', mode
    convert_var_name = trim(var_name)//trim(tmp)
  end procedure convert_var_name


  function n_of_l ( l )
  !! Calculate the number of m-modes for a given value of l. For the scalar
  !! charge in orbit around a Schwarzschild black hole we only need to evolve
  !! modes with \(\ell+m\) even due to symmetries.

    implicit none

    integer(ip), intent(in) :: l
    !! The \(\ell\)-value.
    integer(ip) :: n_of_l
    !! The return value is the number of \(m\)-modes that is evolved for this
    !! \(\ell\).

    if ( mod(l,2) == 0 ) then
      n_of_l = ( l + 2 ) / 2
    else
      n_of_l = ( l + 1 ) / 2
    end if
  end function n_of_l


  function nmodes_of_l ( lmin, lmax )
  !! Calculate the total number of modes from lmin to lmax.
  !!
  !! Uses the function n_of_l.

    implicit none

    integer(ip), intent(in) :: lmin
    !! The minimum value of \(\ell\) to include in the sum,
    !! \(\ell_{\mathrm{min}}\).
    integer(ip), intent(in) :: lmax
    !! The maximum value of \(\ell\) to include in the sum.
    !! \(\ell_{\mathrm{max}}\).
    integer(ip) :: nmodes_of_l
    !! The return value is to total number of \(m\)-modes from all modes from
    !! \(\ell_{\mathrm{min}}\) to \(\ell_{\mathrm{max}}\).
    integer(ip) :: l

    nmodes_of_l = 0

    do l = lmin, lmax
      nmodes_of_l = nmodes_of_l + n_of_l ( l )
    end do
  end function nmodes_of_l


  subroutine set_lm_mode_info ( this, lmin, lmax )
  !! Loop over all mode and store it's corresponding \(\ell\) and m-values.
  !! This should maybe made into a type bound procedure.
  !!
  !! This sets this%[[scal_schw:ll]] and this%[[scal_schw:mm]].

    implicit none

    type(scal_schw), intent(inout) :: this
    !! The mode information will be set for this [[scal_schw]] object.
    integer(ip), intent(in) :: lmin
    !! The mininum \(\ell\)-value.
    integer(ip), intent(in) :: lmax
    !! The mininum \(\ell\)-value.
    integer(ip) :: l, m, nn

    nn = 0
    do l = lmin, lmax
      do m = mod(l,2), l, 2
        nn = nn+1
        this%ll(nn) = l
        this%mm(nn) = m
      end do
    end do
  end subroutine set_lm_mode_info

  subroutine get_elem_flux ( rho, pi, nmodes, mode, r, effs, time_dc, &
                             e_lambda, e_s, e_sinv, e_coeffs, r_elem, &
                             rho_flux, pi_flux )
  !! Calculate the characteristic fluxes for all elements for a single mode.
  !!
  !! Having all of these calculations in a subroutine made it easier to ensure
  !! that all local variables are private (for OpenMP parallelization) and to
  !! avoid segmentation faults.

    use grid, only : drdrho, Tminus_ind, Tplus_ind
    use grid_function
    use parameters, only : order, use_world_tube, use_particle, &
                           use_generic_orbit, q_charge
    use world_tube, only : wt

    implicit none

    type(cgf), intent(in) :: rho
    !! A complex grid function containing the time derivative variable,
    !! \(\Upsilon\).
    type(cgf), intent(in) :: pi
    !! A complex grid function containing the radial derivative variable,
    !! \(\Pi\).
    integer(ip), intent(in) :: nmodes
    !! The total number of modes.
    integer(ip), intent(in) :: mode
    !! The mode for which the flux is calculated.
    type(rgf), intent(in) :: r
    !! A real grid function that contains the Scharzschild radial coordinate.
    type(scal_schw_eff), intent(inout) :: effs
    !! The effective source object. Needed to add and subtract the derivatives
    !! of the singular field when calculating fluxes at the world-tube
    !! boundary.
    type(tdc) :: time_dc
    !! The time dependent coordinate object. Needed to transform the
    !! derivatives of the field when calculating fluxes at the boundary of
    !! the time dependent coordinate region.
    type(rgfb), dimension(:), intent(in) :: e_lambda
    !! A 1d array of size (2) of real boundary grid functions that  contains
    !! the characteristic speeds at the element boundary.
    type(rgfb), dimension(:,:), intent(in) :: e_s
    !! A 2d array of size (2,2) of real boundary grid functions that contains
    !! the matrix used to convert from characteristic to evolved variables.
    type(rgfb), dimension(:,:), intent(in) :: e_sinv
    !! A 2d array of size (2,2) of real boundary grid functions that contains
    !! the matrix used to convert from evolved to characteristic variables.
    type(rgf), dimension(:), intent(in) :: e_coeffs
    !! A 1d array of real grid functions that contains the coefficients for the
    !! wave equation.
    type(ref_element), intent(in) :: r_elem
    !! The reference element.
    type(cgf), intent(inout) :: rho_flux
    !! On output contains the \(\Upsilon\) component of the characteristic
    !! flux.
    type(cgf), intent(inout) :: pi_flux
    !! On output contains the \(\Pi\) component of the characteristic flux.

    integer(ip), dimension(2) :: ind
    integer(ip) :: i, j, l, ne, np, sgn
    complex(wp), dimension(2,2) :: uint, uext
    real(wp), dimension(2,2) :: lambda
    real(wp), dimension(2,2,2) :: s, sinv
    complex(wp), dimension(2,2) :: flux_int
    complex(wp), dimension(1,2) :: psi_sing, psi_sing2
    complex(wp) :: psi_t, psi_rstar, psi_lambda, psi_xi
    logical :: debug_output = .false.

    ne = rho%n
    ind(1) = 1
! Loop over all elements.
    do i = 1, ne
      np = rho%elems(i)%order + 1
      ind(2) = np

! Get the interior data.
      uint(1,1) = rho%elems(i)%var(1)
      uint(2,1) = rho%elems(i)%var(np)
      uint(1,2) = pi%elems(i)%var(1)
      uint(2,2) = pi%elems(i)%var(np)
! Get the exterior data
      if (i>1) then
        uext(1,1) = rho%elems(i-1)%var(np)
        uext(1,2) = pi%elems(i-1)%var(np)
      else
! Zero exterior data at the horizon.
        uext(1,:) = 0.0_wp
      end if
      if (i<ne) then
        uext(2,1) = rho%elems(i+1)%var(1)
        uext(2,2) = pi%elems(i+1)%var(1)
      else
! Zero exterior data at Scri+.`
        uext(2,:) = 0.0_wp
      endif

      ! The following assumes the particle is in the middle of the domain.
      if ( i<=ne/2 ) then
        ! If we are to the left of the particle sgn=-1
        sgn = -1
      else
        ! Otherwise sgn=1
        sgn = 1
      end if
! If we have a world tube and use an effective source we need to subtract
! and/or add the singular field appropriately.
      if ( use_world_tube .and. use_particle ) then
        associate ( r1 => r%elems(i)%var(1), &
                    r2 => r%elems(i)%var(np), &
                    b1 => wt%boundary_info%elems(i)%bvar(1), &
                    b2 => wt%boundary_info%elems(i)%bvar(2) )
! Use the worldtube object's boundary information to determine if we are
! at a worldtube boundary.
          if ( b1 /= izero ) then
            if ( wt%wsize > 0 ) then
              call effs%get_dsingular_dt ( r1, sgn, mode, psi_sing(1,1:) )
              call effs%get_dsingular_dr ( r1, sgn, mode, psi_sing(1,2:) )
            else
              call effs%get_dsingular_dt ( r1, -b1, mode, psi_sing(1,1:) )
              call effs%get_dsingular_dr ( r1, -b1, mode, psi_sing(1,2:) )
              call effs%get_dsingular_dt ( r1, b1, mode, psi_sing2(1,1:) )
              call effs%get_dsingular_dr ( r1, b1, mode, psi_sing2(1,2:) )
              psi_sing(1,:) = psi_sing2(1,:)-psi_sing(1,:)
            end if
! If we are on a generic orbit transform the derivatives of the singular
! field to time dependent coordinates before subtracting or adding (b1 is
! either -1 or 1 depending on whether the current element is inside or
! outside of the worldtube).
            if ( use_generic_orbit ) then
              psi_t = psi_sing(1,1)
              psi_rstar = psi_sing(1,2)
              call time_dc%tortoise_to_tdc ( i, 1, psi_t, &
                                                   psi_rstar, &
                                                   psi_lambda, psi_xi )
              psi_sing(1,1) = psi_lambda
              psi_sing(1,2) = psi_xi
            end if
            uext(1,:) = uext(1,:) + q_charge*b1*psi_sing(1,:)
          end if

! Use the worldtube object's boundary information to determine if we are
! at a worldtube boundary.
          if ( b2 /= izero ) then
            if ( wt%wsize > 0 ) then
              call effs%get_dsingular_dt ( r2, sgn, mode, psi_sing(1,1:) )
              call effs%get_dsingular_dr ( r2, sgn, mode, psi_sing(1,2:) )
            else
              call effs%get_dsingular_dt ( r2, b2, mode, psi_sing(1,1:) )
              call effs%get_dsingular_dr ( r2, b2, mode, psi_sing(1,2:) )
              call effs%get_dsingular_dt ( r2, -b2, mode, psi_sing2(1,1:) )
              call effs%get_dsingular_dr ( r2, -b2, mode, psi_sing2(1,2:) )
              psi_sing(1,:) = psi_sing2(1,:)-psi_sing(1,:)
            end if
! If we are on a generic orbit transform the derivatives of the singular
! field to time dependent coordinates before subtracting or adding (b1 is
! either -1 or 1 depending on whether the current element is inside or
! outside of the worldtube).
            if ( use_generic_orbit ) then
              psi_t = psi_sing(1,1)
              psi_rstar = psi_sing(1,2)
              call time_dc%tortoise_to_tdc ( i, 2, psi_t, &
                                                   psi_rstar, &
                                                   psi_lambda, psi_xi )
              psi_sing(1,1) = psi_lambda
              psi_sing(1,2) = psi_xi
            end if
            uext(2,:) = uext(2,:) + q_charge*b2*psi_sing(1,:)
          end if
        end associate
      end if

! If we are on the boundary of the time dependent coordinate region, we
! have to transform the exterior data (either to or from time dependent
! coordinates) to match the interior data.
      if ( use_generic_orbit ) then
        if ( i == Tminus_ind(1) ) then
          call time_dc%tdc_to_tortoise ( i+1, 1, uext(2,1), &
                                                 uext(2,2), psi_t, &
                                                 psi_rstar )
          uext(2,1) = psi_t
          uext(2,2) = psi_rstar
        end if
        if ( i == Tminus_ind(2) ) then
          call time_dc%tortoise_to_tdc ( i, 1, uext(1,1), &
                                               uext(1,2), psi_lambda, &
                                               psi_xi )
          uext(1,1) = psi_lambda
          uext(1,2) = psi_xi
        end if
        if ( i == Tplus_ind(1) ) then
          call time_dc%tortoise_to_tdc ( i, 2, uext(2,1), &
                                               uext(2,2), psi_lambda, &
                                               psi_xi )
          uext(2,1) = psi_lambda
          uext(2,2) = psi_xi
        end if
        if ( i == Tplus_ind(2) ) then
          call time_dc%tdc_to_tortoise ( i-1, 2, uext(1,1), &
                                                 uext(1,2), psi_t, &
                                                 psi_rstar )
          uext(1,1) = psi_t
          uext(1,2) = psi_rstar
        end if
      end if

! Store the eigen-data in local variables.
      do j = 1, 2
        lambda(:,j) = e_lambda(j)%elems(i)%bvar(:)
        do l = 1, 2
          s(:,l,j) = e_s(l,j)%elems(i)%bvar(:)
          sinv(:,l,j) = e_sinv(l,j)%elems(i)%bvar(:)
        end do
! Calculate the interior flux.
        flux_int(j,:) = flux ( uint(j,:), &
                               e_coeffs(2)%elems(i)%var(ind(j)), &
                               e_coeffs(1)%elems(i)%var(ind(j)) )
      end do

! Pass the information to the characteristic flux routine from the
! reference element object and get the numerical flux back (the lift matrix
! has already been applied.
      flux_result(:,:,mode) = r_elem%characteristic_flux ( 2, order, uint, uext, &
                                                 flux_int, lambda, s, sinv, &
                                                 debug_output )

! The only thing left to do is to store the result after multiplying with the
! Jacobian.
      do j = 1, np
        rho_flux%elems(i)%var(j) = drdrho%elems(i)%var(j)*flux_result(j,1,mode)
        pi_flux%elems(i)%var(j) = drdrho%elems(i)%var(j)*flux_result(j,2,mode)
      end do
    end do

  end subroutine get_elem_flux


  function linear_extrapolate ( x, y )
  !! Utility function for linear extrapolation of data to the horizon and
  !! Scri+.
  !!
  !! Used to import external initial data from files that do not contain valid
  !! data at those locations.

    implicit none

    real(wp), dimension(3), intent(in) :: x
    !! A 1d array of size (3) of reals that contains the 3 coordinate locations
    !! involved in the extrapolation. The first element is the coordinate
    !! where the extrapolation is needed, that last two elements are the
    !! locations where data is known.
    complex(wp), dimension(2:3), intent(in) :: y
    !! A 1d array of size(2) that contains the known data.
    complex(wp) :: linear_extrapolate
    !! The return value is the extrapolated data at the required location.

    linear_extrapolate = y(2) - (y(3)-y(2))*(x(2)-x(1))/(x(3)-x(2))

  end function linear_extrapolate
end submodule scalar_schw_implementation
