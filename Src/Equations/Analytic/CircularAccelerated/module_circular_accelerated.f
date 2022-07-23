! The expressions gets filled in using
! /home/diener/Mathematica_new/Non-uniform_circular_worldline_corrected_dt.nb
! Don't execute the whole notebook as it takes way too long to simplify some
! expressions that are not needed.
      module accelerated_circular_orbit

      use kinds
      use parameters, only : use_gaussian_acceleration

      implicit none

      contains

      subroutine circ_accel ( t, t0, r0, mass, sigma, amp, phi, En, Lz,
     -                        ar, aphi, dardt, daphidt, 
     -                        d2ardt2, d2aphidt2 )

        use kinds
        
        real(qp), intent(in) :: t
        real(wp), intent(in) :: t0, r0, mass, sigma, amp
        real(wp), intent(out) :: phi, En, Lz, ar, aphi
        real(wp), intent(out) :: dardt, daphidt, d2ardt2, d2aphidt2
        real(wp) :: dt, omega0, f0, omega, omegap, omegapp, omegappp
        real(qp) :: phitmp
        real(qp), parameter :: twopiq = 2.0_qp*acos(-1.0_qp)
        real(wp), parameter :: pi = acos(-1.0_wp)

        dt = t - t0
        f0 = 1 - (2*mass)/r0
        omega0 = Sqrt(mass/r0**3)
        if (use_gaussian_acceleration) then
          omega = omega0 + amp*omega0*exp(-(dt**2/sigma**2))
          omegap = (-2*amp*dt*omega0*exp(-(dt**2/sigma**2)))/sigma**2
          omegapp =         (-2*amp*omega0*(-2*dt**2 + sigma**2)*
     -    exp(-(dt**2/sigma**2)))/sigma**4
          omegappp =         (4*amp*omega0*(-2*dt**3 + 3*dt*sigma**2)*
     -    exp(-(dt**2/sigma**2)))/sigma**6
          phitmp =         omega0*t + (amp*omega0*Sqrt(Pi)*sigma*
     -     Erf(dt/sigma))/2. + 
     -  (amp*omega0*Sqrt(Pi)*sigma*Erf(t0/sigma))/2.
        else
          omega = omega0 - (2*amp*dt*exp(-(dt**2/sigma**2)))/sigma**2
          omegap =         (-2*amp*(-2*dt**2 + sigma**2)*
     -    exp(-(dt**2/sigma**2)))/sigma**4
          omegapp =         (-4*amp*dt*(2*dt**2 - 3*sigma**2)*
     -    exp(-(dt**2/sigma**2)))/sigma**6
          omegappp =         (4*amp*(4*dt**4 - 12*dt**2*sigma**2 + 
     -      3*sigma**4)*exp(-(dt**2/sigma**2)))/sigma**8
          phitmp =         omega0*t + amp*exp(-(dt**2/sigma**2)) - 
     -  amp*exp(-(t0**2/sigma**2))
        end if
        phi = real(mod(phitmp,twopiq),wp)
        En = f0/Sqrt(f0 - omega**2*r0**2)
        Lz = (omega*r0**2)/Sqrt(f0 - omega**2*r0**2)
        ar =         (f0*(mass - omega**2*r0**3))/
     -  (r0**2*(f0 - omega**2*r0**2))
        aphi = (f0*omegap)/(f0 - omega**2*r0**2)**2
        dardt =         (2*f0*omega*omegap*(mass - f0*r0))/
     -  (f0 - omega**2*r0**2)**2
        daphidt =         (f0*(4*omega*omegap**2*r0**2 + 
     -      omegapp*(f0 - omega**2*r0**2)))/
     -  (f0 - omega**2*r0**2)**3
        d2ardt2 =         (2*f0*(mass - f0*r0)*
     -    (omega*omegapp*(f0 - omega**2*r0**2) + 
     -      omegap**2*(f0 + 3*omega**2*r0**2)))/
     -  (f0 - omega**2*r0**2)**3
        d2aphidt2 =         (f0*omegappp*(f0 - omega**2*r0**2)**2 + 
     -    4*f0*omegap*r0**2*
     -     (3*omega*omegapp*(f0 - omega**2*r0**2) + 
     -       omegap**2*(f0 + 5*omega**2*r0**2)))/
     -  (f0 - omega**2*r0**2)**4
        
      end subroutine circ_accel

      end module accelerated_circular_orbit
