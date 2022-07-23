submodule(scalar_schw_eff) scalar_schw_eff_implementation
!! The implementation of the interfaces provided in [[scalar_schw_eff]].

  implicit none

  interface
    subroutine initialize_source ( nmodes, l, m, mass ) bind(c, name='init_source')
    !! Interface to the C++ routine init_source for initializing the C++
    !! effective source class.
      use iso_c_binding
      integer(c_int) :: nmodes
      !! The number of modes.
      integer(c_int), dimension(nmodes) :: l
      !! A 1d array that contains the \(\ell\)-values.
      integer(c_int), dimension(nmodes) :: m
      !! A 1d array that contains the \(m\)-values.
      real(c_double) :: mass
      !! The mass of the black hole.
    end subroutine initialize_source

    subroutine set_time_window_coeffs ( tfac, dtfac_dt, d2tfac_dt2, &
                                        do_smooth_after_lmax, nmodes ) &
                                        bind(c, name='set_time_window')
    !! Interface to the C++ routine set_time_window that sets the time window
    !! information in the C++ effective source class.
      use iso_c_binding
      real(c_double) :: tfac
      !! The time window factor.
      real(c_double) :: dtfac_dt
      !! The first time derivative of the time window factor.
      real(c_double) :: d2tfac_dt2
      !! The second time derivative of the time window factor.
      integer(c_int) :: nmodes
      !! The number of modes.
      integer(c_int) :: do_smooth_after_lmax
      !! Only set the time window for
      !! \(\ell>\)do_smooth_after_lmax.
    end subroutine set_time_window_coeffs

    subroutine set_particle_pos ( r, phi, ur, En, Lz, ar, aphi, dardt, &
                                  daphidt, d2ardt2, d2aphidt2, &
                                  nmodes ) bind(c, name='set_particle')
    !! Interface to the C++ routine set_particle that sets the particle
    !! information in the C++ effective source class.
      use iso_c_binding
      real(c_double) :: r
      !! The radial coordinate of the particle (in Schwarzschild coordinates),
      !! \(r\).
      real(c_double) :: phi
      !! The azimuthal angle of the particle, \(\phi\).
      real(c_double) :: ur
      !! The radial component of the 4-velocity, \(u^r\).
      real(c_double) :: En
      !! The energy per unit mass of the particle, \(E\).
      real(c_double) :: Lz
      !! The angular momentum per unit mass of the particle \(L_z\).
      real(c_double) :: ar
      !! The radial component of the 4-acceleration, \(a^r\).
      real(c_double) :: aphi
      !! The \(\phi\)-component of the 4-acceleration, \(a^{\phi}\).
      real(c_double) :: dardt
      !! The derivative of \(a^r\) with respect to Schwarzschild coordinate
      !! time.
      real(c_double) :: daphidt
      !! The derivative of \(a^{\phi}\) with respect to Schwarzschild
      !! coordinate time.
      real(c_double) :: d2ardt2
      !! The second derivative of \(a^r\) with respect to Schwarzschild
      !! coordinate time.
      real(c_double) :: d2aphidt2
      !! The second derivative of \(a^{\phi}\) with respect to Schwarzschild
      !! coordinate time.
      integer(c_int) :: nmodes
      !! The number of modes.
    end subroutine set_particle_pos

    subroutine evaluate_source_all ( mode, n, r, win, dwin, d2win, sre, sim ) &
                                      bind(c, name='eval_source_all')
    !! Interface to the C++ routine eval_source_all that evaluates the
    !! effective source on the locations provided in a 1d input array.
      use iso_c_binding
      integer(c_int) :: mode
      !! The mode.
      integer(c_int) :: n
      !! The number of elements in the input radial coordinate array, \(r\).
      real(c_double) :: r(*)
      !! The radial coordinates in Schwarzschild coordinates, \(r_i\).
      real(c_double) :: win(*)
      !! The window function, \(W_i\).
      real(c_double) :: dwin(*)
      !! The radial derivative of the window function,
      !! \( \left .\frac{dW}{dr}\right |_i \).
      real(c_double) :: d2win(*)
      !! The second radial derivative of the window function,
      !! \(\left .\frac{d^2W}{dr^2}\right |_i\).
      real(c_double) :: sre(*)
      !! The real part of the effective source,
      !! \(\mathcal{R}(S_{\mathrm{eff}})\).
      real(c_double) :: sim(*)
      !! The imaginary part of the effective source,
      !! \(\mathcal{I}(S_{\mathrm{eff}})\).
    end subroutine  evaluate_source_all

    subroutine get_phi ( mode, r, phire, phiim ) bind(c, name='Phi')
    !! Interface to the C++ routine Phi that evaluates the singuler field
    !! at a given radial coordinate in Schwarzschild coordinates.
      use iso_c_binding
      integer(c_int) :: mode
      !! The mode.
      real(c_double) :: r
      !! The radial coordinate, \(r\), in Schwarzschild coordinates.
      real(c_double) :: phire
      !! On output the real part of the singular field \(\cal{R}(\Phi(r))\).
      real(c_double) :: phiim
      !! On output the imaginary part of the singular field
      !! \(\cal{I}(\Phi(r))\).
    end subroutine get_phi

    subroutine get_dphidt ( mode, r, sgn, dphidtre, dphidtim ) &
                            bind(c, name='dPhi_dt')
    !! Interface to the C++ routine dPhi_dt that evaluates the time derivative
    !! of the singuler field at a given radial coordinate in Schwarzschild
    !! coordinates.
      use iso_c_binding
      integer(c_int) :: mode
      !! The mode.
      real(c_double) :: r
      !! The radial coordinate, \(r\), in Schwarzschild coordinates.
      integer(c_int) :: sgn
      !! The value should be 1 when the derivative to the right of the
      !! particle is needed and -1 when the derivative to the left of the
      !! particle is needed.
      real(c_double) :: dphidtre
      !! On output the real part of the time derivative of the singular field
      !! \(\cal{R}\left (\left .\frac{d\Phi}{dt}\right |_r\right )\).
      real(c_double) :: dphidtim
      !! On output the imaginary part of the time derivative of the singular
      !! field \(\cal{R}\left (\left .\frac{d\Phi}{dt}\right |_r\right ) \).
    end subroutine get_dphidt

    subroutine get_dphidr ( mode, r, sgn, dphidrre, dphidrim ) &
                            bind(c, name='dPhi_dr')
    !! Interface to the C++ routine dPhi_dt that evaluates the radial
    !! derivative of the singuler field at a given radial coordinate in
    !! Schwarzschild coordinates.
      use iso_c_binding
      integer(c_int) :: mode
      !! The mode.
      real(c_double) :: r
      !! The radial coordinate, \(r\), in Schwarzschild coordinates.
      integer(c_int) :: sgn
      !! The value should be 1 when the derivative to the right of the
      !! particle is needed and -1 when the derivative to the left of the
      !! particle is needed.
      real(c_double) :: dphidrre
      !! On output the real part of the radial derivative of the singular field
      !! \(\cal{R}\left (\left .\frac{d\Phi}{dr}\right |_r\right )\).
      real(c_double) :: dphidrim
      !! On output the imaginary part of the radial derivative of the singular
      !! field \(\cal{R}\left (\left .\frac{d\Phi}{dr}\right |_r\right )\).
    end subroutine get_dphidr
  end interface

contains

  module procedure scal_schw_eff_init

    use iso_c_binding
    use parameters, only : n_elems, order
    use scalar_schw, only : convert_var_name

    implicit none

    integer :: allocation_status
    integer(c_int) :: nmodes_c
    real(c_double) :: mass_c
    integer(ip) :: np, i, j
    character(len=6) :: source_name = 'source'

    this%nvars = nvars/3
    this%nmodes = nmodes
!    print*,'nvars = ', this%nvars
!    print*,'nmodes = ', this%nmodes
!    print*,'l = ', l
!    print*,'m = ', m
    allocate ( this%source(this%nvars,nmodes), stat=allocation_status )
    if ( allocation_status > 0 ) then
      print*,'Allocation of Schwarzschild scalar source failed.'
      stop
    end if

    do i = 1, this%nmodes
      do j = 1, this%nvars
        this%source(j,i) = cgf ( n_elems, order, &
                                 convert_var_name(source_name,i) )
      end do
    end do

    nmodes_c = this%nmodes
    mass_c = mass

    call initialize_source ( nmodes_c, l, m, mass_c )

    np = this%source(1,1)%elems(1)%order+1
    
    allocate ( this%sre(np,nmodes), this%sim(np,nmodes) )

  end procedure scal_schw_eff_init


  module procedure scal_schw_eff_set_particle_pos

    use iso_c_binding

    implicit none

    real(c_double) :: r_c, phi_c, ur_c, En_c, Lz_c, ar_c, aphi_c, &
                      dardt_c, daphidt_c, d2ardt2_c, d2aphidt2_c

    r_c = r
    phi_c = phi
    ur_c = ur
    En_c = En
    Lz_c = Lz
    ar_c = ar
    aphi_c = aphi
    dardt_c = dardt
    daphidt_c = daphidt
    d2ardt2_c = d2ardt2
    d2aphidt2_c = d2aphidt2

    call set_particle_pos ( r_c, phi_c, ur_c, En_c, Lz_c, ar_c, aphi_c, &
                            dardt_c, daphidt_c, d2ardt2_c, d2aphidt2_c, &
                            this%nmodes ) 

  end procedure scal_schw_eff_set_particle_pos


  module procedure scal_schw_eff_set_time_window

    use iso_c_binding

    implicit none

    real(c_double) :: tfac_c, dtfac_dt_c, d2tfac_dt2_c
    integer(c_int) :: do_smooth_after_lmax_c

    tfac_c = tfac
    dtfac_dt_c = dtfac_dt
    d2tfac_dt2_c = d2tfac_dt2
    do_smooth_after_lmax_c = do_smooth_after_lmax

    call set_time_window_coeffs ( tfac_c, dtfac_dt_c, d2tfac_dt2_c, &
                                  do_smooth_after_lmax_c, this%nmodes )

  end procedure scal_schw_eff_set_time_window


  module procedure scal_schw_eff_evaluate_source

    use iso_c_binding

    implicit none

    integer(ip) :: i, j, k, l, nmodes, nvars, nelems
    integer(c_int) :: i_c, np

    nmodes = this%nmodes
    nvars = this%nvars 
    nelems = this%source(1,1)%n
    np = this%source(1,1)%elems(1)%order+1

!$OMP PARALLEL DO private(i,j,k,i_c,l) shared(r,wt,this)
    do i = 1, nmodes
      i_c = i
      do j = 1, nvars
        do k = wt%windex1, wt%windex2
          call evaluate_source_all ( i-1, np, r%elems(k)%var, &
                                     wt%win%elems(k)%var, &
                                     wt%dwin%elems(k)%var, &
                                     wt%d2win%elems(k)%var, &
                                     this%sre(:,i), this%sim(:,i) )
          do l = 1, np
            this%source(j,i)%elems(k)%var(l) = cmplx ( this%sre(l,i), &
                                                       this%sim(l,i), wp )
          end do
        end do

        do k = 1, wt%windex1-1
          this%source(j,i)%elems(k)%var = czero
        end do
        do k = wt%windex2+1, nelems
          this%source(j,i)%elems(k)%var = czero
        end do
      end do
    end do
!$OMP END PARALLEL DO

  end procedure scal_schw_eff_evaluate_source


  module procedure scal_schw_eff_get_singular

    use iso_c_binding

    implicit none

    integer(ip) :: i, j, k, l, nmodes, nvars, nelems
    real(wp) :: phire, phiim

    call get_phi ( mode-1, r, phire, phiim ) 
    psi(1) = cmplx ( phire, phiim, wp )

  end procedure scal_schw_eff_get_singular


  module procedure scal_schw_eff_get_dsingular_dt

    use iso_c_binding

    implicit none

    integer(ip) :: i, j, k, l, nmodes, nvars, nelems
    real(wp) :: dphidtre, dphidtim

    call get_dphidt ( mode-1, r, sgn, dphidtre, dphidtim ) 
    dpsidt(1) = cmplx ( dphidtre, dphidtim, wp )

  end procedure scal_schw_eff_get_dsingular_dt


  module procedure scal_schw_eff_get_dsingular_dr

    use parameters, only : mass
    use iso_c_binding

    implicit none

    integer(ip) :: i, j, k, l, nmodes, nvars, nelems
    real(wp) :: dphidrre, dphidrim, rfac

!    print*,'mode = ', mode
!    print*,'r = ', r
    call get_dphidr ( mode-1, r, sgn, dphidrre, dphidrim ) 
!    print*,'dphidr = ', dphidrre, dphidrim
! Convert to r_star coordinates
    rfac = (r-2.0_wp*mass)/r
    dpsidr(1) = rfac*cmplx ( dphidrre, dphidrim, wp )

  end procedure scal_schw_eff_get_dsingular_dr

end submodule scalar_schw_eff_implementation
