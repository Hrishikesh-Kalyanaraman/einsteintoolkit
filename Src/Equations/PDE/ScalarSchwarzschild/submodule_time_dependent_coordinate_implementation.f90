submodule(time_dependent_coordinate) time_dependent_coordinate_implementation
!! Implementation of the interfaces defined in [[time_dependent_coordinate]].

  implicit none

contains

  module procedure tdc_init

    use parameters, only : n_elems, order

    implicit none

    this%dxdlambda = rgf ( n_elems, order, 'dxdlambda' )
    this%dxdxi = rgf ( n_elems, order, 'dxdxi' )
    this%d2xdlambda2 = rgf ( n_elems, order, 'd2xdlambda2' )
    this%d2xdxi2 = rgf ( n_elems, order, 'd2xdxi2' )
    this%d2xdlambdadxi = rgf ( n_elems, order, 'd2xdlambdadxi' )
    this%dxdlambda_b = rgfb ( n_elems, 'dxdlambda_b' )
    this%dxdxi_b = rgfb ( n_elems, 'dxdxi_b' )
    this%rm2m = rgf ( n_elems, order, 'rm2m' )
  end procedure tdc_init


  module procedure tdc_set_coefficients

    use parameters, only : mass
    use grid, only : Tminus_ind, Tplus_ind, Tminus, Tplus, rho_particle
    use orbit_base, only : tdc_info
    use numerics, only : invert_tortoise
   
    implicit none

    integer(ip) :: i, j, k, nnodes, nmodes
    real(wp) :: r, drdt, d2rdt2
    real(wp) :: dxdxiinv, dxdxiinv2, dxdxiinv3, dxdlambda2, rfac

    this%maxspeed = 1.0_wp

    call tdc_info%get_tdc ( r, drdt, d2rdt2 )
    call convert_rschw_to_rstar ( r, drdt, d2rdt2 )

    do i = Tminus_ind(2), Tplus_ind(1)
      nnodes = this%dxdlambda%elems(i)%order + 1
      associate ( dxdlambda => this%dxdlambda%elems(i)%var, &
                  dxdxi => this%dxdxi%elems(i)%var, &
                  dxdlambda_b => this%dxdlambda_b%elems(i)%bvar, &
                  dxdxi_b => this%dxdxi_b%elems(i)%bvar, &
                  d2xdlambda2 => this%d2xdlambda2%elems(i)%var, &
                  d2xdxi2 => this%d2xdxi2%elems(i)%var, &
                  d2xdlambdadxi => this%d2xdlambdadxi%elems(i)%var, &
                  cxixi => coeffs(1)%elems(i)%var, &
                  clambdaxi => coeffs(2)%elems(i)%var, &
                  cxi => coeffs(3)%elems(i)%var, &
                  clambda => coeffs(4)%elems(i)%var, &
                  rho => rho%elems(i)%var, &
                  rstar => rstar%elems(i)%var, &
                  rschw => rschw%elems(i)%var, &
                  rm2m => this%rm2m%elems(i)%var, &
                  maxspeed => this%maxspeed )
        call coord_trans ( Tminus, Tplus, r, rho_particle, drdt, d2rdt2, &
                           rho, rstar, dxdlambda, dxdxi, d2xdlambda2, &
                           d2xdxi2, d2xdlambdadxi ) 
        dxdlambda_b(1) = dxdlambda(1)
        dxdlambda_b(2) = dxdlambda(nnodes)
        dxdxi_b(1) = dxdxi(1)
        dxdxi_b(2) = dxdxi(nnodes)

        do j = 1, nnodes
          dxdxiinv = 1.0_wp / dxdxi(j)
          dxdxiinv2 = dxdxiinv*dxdxiinv
          dxdxiinv3 = dxdxiinv2*dxdxiinv
          dxdlambda2 = dxdlambda(j)*dxdlambda(j)
          
          cxixi(j) = ( 1.0_wp - dxdlambda2 ) * dxdxiinv2
          clambdaxi(j) = 2.0_wp * dxdlambda(j) * dxdxiinv
          cxi(j) = ( d2xdxi2(j) * ( dxdlambda2 - 1.0_wp ) &
                 +dxdxi(j) * ( -2.0_wp*dxdlambda(j)*d2xdlambdadxi(j) &
                               +dxdxi(j)*d2xdlambda2(j) ) ) * dxdxiinv3
          clambda(j) = 0.0_wp
          if (j==1) then
            k = 1
          else if (j==nnodes) then
            k = 2
          end if
          if (j==1 .or. j==nnodes) then
            lambda(1)%elems(i)%bvar(k) = -(1.0_wp+dxdlambda(j))*dxdxiinv
            lambda(2)%elems(i)%bvar(k) = (1.0_wp-dxdlambda(j))*dxdxiinv
            maxspeed = max ( abs(lambda(1)%elems(i)%bvar(k)), &
                             abs(lambda(2)%elems(i)%bvar(k)), &
                             maxspeed )
            s(1,1)%elems(i)%bvar(k) = (1.0_wp+dxdlambda(j))*dxdxiinv
            s(1,2)%elems(i)%bvar(k) = (-1.0_wp+dxdlambda(j))*dxdxiinv
            s(2,1)%elems(i)%bvar(k) = 1.0_wp
            s(2,2)%elems(i)%bvar(k) = 1.0_wp
            sinv(1,1)%elems(i)%bvar(k) = 0.5_wp*dxdxi(j)
            sinv(1,2)%elems(i)%bvar(k) = 0.5_wp*(1.0_wp-dxdlambda(j))
            sinv(2,1)%elems(i)%bvar(k) = -0.5_wp*dxdxi(j)
            sinv(2,2)%elems(i)%bvar(k) = 0.5_wp*(1.0_wp+dxdlambda(j))
          end if

          rm2m(j) = invert_tortoise ( rstar(j), mass )

          rschw(j) = 2.0_wp*mass + rm2m(j)
        end do
      end associate
    end do    

    nmodes = size(ll)

!$OMP PARALLEL DO private(k,i,nnodes,j,rfac) shared(nmodes,TMinus,Tplus,this)
    do k = 1, nmodes
      do i = Tminus_ind(2), Tplus_ind(1)
        nnodes = this%dxdlambda%elems(i)%order + 1
        do j = 1, nnodes
          associate ( rm2m => this%rm2m%elems(i)%var(j), &
                      rschw => rschw%elems(i)%var(j), &
                      lcoeffs => lcoeffs(k)%elems(i)%var(j) )
            rfac = rm2m / rschw**4 
            lcoeffs = -(ll(k)*(ll(k)+1)*rschw+2.0_wp*mass)*rfac
          end associate
        end do
      end do
    end do
!$OMP END PARALLEL DO
  end procedure tdc_set_coefficients


  module procedure tdc_tdc_to_tortoise_cvec

    implicit none

    associate ( drstardxi => this%dxdxi_b%elems(ielem)%bvar(dir), &
                drstardlambda => this%dxdlambda_b%elems(ielem)%bvar(dir) )
      dudrstar = dudxi/drstardxi
      dudt = dudlambda - drstardlambda*dudrstar
    end associate
    
  end procedure tdc_tdc_to_tortoise_cvec


  module procedure tdc_tdc_to_tortoise_cscal

    implicit none

    associate ( drstardxi => this%dxdxi_b%elems(ielem)%bvar(dir), &
                drstardlambda => this%dxdlambda_b%elems(ielem)%bvar(dir) )
      dudrstar = dudxi/drstardxi
      dudt = dudlambda - drstardlambda*dudrstar
    end associate
    
  end procedure tdc_tdc_to_tortoise_cscal


  module procedure tdc_tdc_to_tortoise_rscal

    implicit none

    associate ( drstardxi => this%dxdxi_b%elems(ielem)%bvar(dir), &
                drstardlambda => this%dxdlambda_b%elems(ielem)%bvar(dir) )
      dudrstar = dudxi/drstardxi
      dudt = dudlambda - drstardlambda*dudrstar
    end associate
    
  end procedure tdc_tdc_to_tortoise_rscal


  module procedure tdc_tortoise_to_tdc_cvecb

    implicit none

    associate ( drstardxi => this%dxdxi_b%elems(ielem)%bvar(dir), &
                drstardlambda => this%dxdlambda_b%elems(ielem)%bvar(dir) )
      dudlambda = dudt + drstardlambda*dudrstar
      dudxi = drstardxi*dudrstar
    end associate
    
  end procedure tdc_tortoise_to_tdc_cvecb


  module procedure tdc_tortoise_to_tdc_cscalb

    implicit none

    associate ( drstardxi => this%dxdxi_b%elems(ielem)%bvar(dir), &
                drstardlambda => this%dxdlambda_b%elems(ielem)%bvar(dir) )
      dudlambda = dudt + drstardlambda*dudrstar
      dudxi = drstardxi*dudrstar
    end associate
    
  end procedure tdc_tortoise_to_tdc_cscalb


  module procedure tdc_tortoise_to_tdc_cscal

    implicit none

    associate ( dxdlambda => this%dxdlambda%elems(elem)%var(node), &
                dxdxi => this%dxdxi%elems(elem)%var(node) )
      dpsidt = dpsidt + dxdlambda*dpsidr
      dpsidr = dxdxi*dpsidr
    end associate 

  end procedure tdc_tortoise_to_tdc_cscal


  subroutine coord_trans ( a, b, xp, xip, dxpdt, d2xpdt2, xi, &
                           x, dxdt, dxdxi, d2xdt2, d2xdxi2, d2xdtdxi )
  !! Routine that calculates the Tortoise coordinates \((t,r_*)\) from time
  !! dependent coordinates \((\lambda,\xi)\) where the particle is kept at a
  !! fixed coordinate location as well as the informtation needed to transform
  !! the wave equation to time dependent coordinates.

    implicit none

    real(wp), intent(in) :: a
    !! The lower boundary location of the time dependent coordinate region.
    real(wp), intent(in) :: b
    !! The upper boundary location of the time dependent coordinate region.
    real(wp), intent(in) :: xp
    !! The current particle location in Tortoise coordinates,
    !! \(r^{\mathrm{p}}_*\).
    real(wp), intent(in) :: xip
    !! The constant particle location in time dependent coordinates,
    !! \(\xi^{\mathrm{p}}\).
    real(wp), intent(in) :: dxpdt
    !! The current time derivative of the particle location in Tortoise
    !! coordinates, \(\frac{d r_*^{\mathrm{p}}}{dt}\).
    real(wp), intent(in) :: d2xpdt2
    !! The current second time derivative of the particle location in Tortoise
    !! coordinates, \(\frac{d^2 r_*^{\mathrm{p}}}{dt^2}\).
    real(wp), dimension(:), intent(in) :: xi
    !! A 1d array containing coordinate values \(\xi_i\).
    real(wp), dimension(:), intent(out) :: x
    !! A 1d array that on output contains the Tortoise coordinates \((r_*)_i\).
    real(wp), dimension(:), intent(out) :: dxdt
    !! A 1d array that on outout contains \(\frac{d r_*}{dt}\).
    real(wp), dimension(:), intent(out) :: dxdxi
    !! A 1d array that on outout contains \(\frac{d r_*}{d\xi}\).
    real(wp), dimension(:), intent(out) :: d2xdt2
    !! A 1d array that on outout contains \(\frac{d^2 r_*}{dt^2}\).
    real(wp), dimension(:), intent(out) :: d2xdxi2
    !! A 1d array that on outout contains \(\frac{d^2 r_*}{d\xi^2}\).
    real(wp), dimension(:), intent(out) :: d2xdtdxi
    !! A 1d array that on outout contains \(\frac{d^2 r_*}{dt d\xi}\).
    integer(ip) :: i
    real(wp) :: xpma, xipma, xima, bmxp, bmxip, bmxi, bma, ximxip, xipmxp
    real(wp) :: xipmainv, xipmamulbmxipinv, dtfac

    do i = 1, size(xi)
      xpma = xp - a
      xipma = xip - a
      xima = xi(i) - a
      bmxp = b - xp
      bmxi = b - xi(i)
      bmxip = b - xip
      bma = b - a
      ximxip = xi(i) - xip
      xipmainv = 1.0_wp/xipma
      xipmamulbmxipinv = xipmainv/bmxip
      xipmxp = xip - xp
      x(i) = a + xpma*xipmainv*xima &
           +( bmxp*xipma-xpma*bmxip )*xipmamulbmxipinv/bma*xima*ximxip
      dtfac = xima*bmxi*xipmamulbmxipinv
      dxdt(i) = dtfac*dxpdt
      dxdxi(i) = ( (2.0_wp*xi(i)-xip-a)*xipmxp + xpma*bmxip ) &
                 * xipmamulbmxipinv
      d2xdt2(i) = dtfac*d2xpdt2
      d2xdxi2(i) = 2.0_wp*xipmxp*xipmamulbmxipinv
      d2xdtdxi(i) = ( a + b - 2.0_wp * xi(i) ) * xipmamulbmxipinv * dxpdt
    end do
  end subroutine coord_trans


  subroutine convert_rschw_to_rstar ( rp, drpdt, d2rpdt2 )
  !! Routine that converts the particle location and time derivatives in
  !! Schwarszschild coordinates to Tortoise coordinates.

    use numerics, only : rstar_of_r
    use parameters, only : mass

    implicit none

    real(wp), intent(inout) :: rp
    !! On input the particle location in Schwarschild coordinates,
    !! \(r^{\mathrm{p}}\). On output the particle location in Tortoise
    !! coordinates, \(r_*^{\mathrm{p}}\).
    real(wp), intent(inout) :: drpdt
    !! On input the time derivative of the particle location in Schwarschild
    !! coordinates, \(\frac{dr^{\mathrm{p}}}{dt}\). On output the time
    !! derivative of the particle location in Tortoise coordinates,
    !! \(\frac{dr_*^{\mathrm{p}}}{dt}\).
    real(wp), intent(inout) :: d2rpdt2
    !! On input the second time derivative of the particle location in
    !! Schwarschild coordinates, \(\frac{d^2r^{\mathrm{p}}}{dt^2}\). On output
    !! the second time derivative of the particle location in Tortoise
    !! coordinates, \(\frac{d^2r_*^{\mathrm{p}}}{dt^2}\).

    d2rpdt2 = ( -2.0_wp * mass * drpdt**2 &
                + rp * ( rp - 2.0_wp * mass ) * d2rpdt2 ) &
              / ( rp - 2.0_wp * mass )**2
    drpdt = rp / ( rp - 2.0_wp * mass ) * drpdt
    rp = rstar_of_r ( rp, mass )
  end subroutine
end submodule time_dependent_coordinate_implementation
