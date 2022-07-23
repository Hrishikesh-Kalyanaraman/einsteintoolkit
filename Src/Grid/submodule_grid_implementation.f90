submodule(grid) grid_implementation
!! The implementation of the interfaces defined in [[grid]].

contains

  module procedure init_grid_coordinates

    use DG_structures
    use parameters
    use numerics

    implicit none

    integer(ip) :: i, j
    real(wp) :: nk
    real(wp), dimension(2) :: vx
    real(wp) :: rstar_center, delta_rstar, rmax, rmin

    rho_min = Sminus
    print*,'rho_min = ', rho_min
    if (use_particle) then
      rmin = p_orb*mass/(1.0_wp+ecc)
      print*,'rmin = ', rmin
      rmax = p_orb*mass/(1.0_wp-ecc)
      print*,'rmax = ', rmax
      rstar_center = rstar_of_r ( 0.5_wp*(rmin+rmax), mass )
      print*,'rstar_center = ', rstar_center
      rho_particle = rstar_center
    else
      rstar_center = rstar_of_r ( r_center, mass )
      print*,'rstar_center = ', rstar_center
    end if
    delta_rstar = ( rstar_center - Sminus )*2.0_wp / n_elems
    print*,'delta_rstar = ', delta_rstar
    rho_max = rstar_center + nint(0.5_wp*n_elems)*delta_rstar
    print*,'rho_max = ', rho_max
    Tminus = rstar_center - t_size*delta_rstar
    Tminus_ind(1) = n_elems/2 - t_size
    Tminus_ind(2) = n_elems/2 - t_size+1
    print*,'Tminus = ', Tminus
    print*,'Tminus_ind = ', Tminus_ind
    Tplus = rstar_center + t_size*delta_rstar
    Tplus_ind(1) = n_elems/2 + t_size
    Tplus_ind(2) = n_elems/2 + t_size+1
    print*,'Tplus = ', Tplus
    print*,'Tplus_ind = ', Tplus_ind
    nk = real(n_elems,wp)
      
    rho = rgf ( n_elems, order, 'rho' )
    drhodr = rgf ( n_elems, order, 'drhodr' )
    drdrho = rgf ( n_elems, order, 'drdrho' )

    if (use_particle) then
      particle_element(1) = n_elems/2
      particle_element(2) = n_elems/2+1
      particle_node(1) = order+1
      particle_node(2) = 1
    else
      particle_element = -1
      particle_node = -1
    end if

!$OMP PARALLEL DO private(i,vx) shared(rho_min,rho_max,nk,relem, &
!$OMP                               order,rho,drhodr,drdrho)
    do i = 1, n_elems
      vx(1) = (rho_max-rho_min)*real(i-1,wp)/nk + rho_min
      vx(2) = (rho_max-rho_min)*real(i,wp)/nk + rho_min

      call Jacobian ( order, vx, relem%r, relem%dr, rho%elems(i)%var, &
                      drhodr%elems(i)%var, drdrho%elems(i)%var )
    end do
!$OMP END PARALLEL DO

    delta_rho_min = rho%elems(1)%var(2) - rho%elems(1)%var(1)

    print*,'rho(Tminus_ind(1)) = ', rho%elems(Tminus_ind(1))%var(order+1)
    print*,'rho(Tminus_ind(2)) = ', rho%elems(Tminus_ind(2))%var(1)
    print*,'rho(Tplus_ind(1)) = ', rho%elems(Tplus_ind(1))%var(order+1)
    print*,'rho(Tplus_ind(2)) = ', rho%elems(Tplus_ind(2))%var(1)
    print*,'delta_rho_min = ', delta_rho_min

  end procedure init_grid_coordinates


  subroutine Jacobian ( n, vx, r, dr, x, xr, rx )
  !! Routine that sets up the physical coordinates as well as the Jacobian and
  !! inverse Jacobian to convert derivatives between reference element and
  !! physical element.
  !!
  !! This is almost exactly the same as the routine of the same name in
  !! [[submodule_DG_implementation.f90]] except it also returns xr. Should
  !! probably be combined. I don't think xr is used anywhere so maybe this is
  !! not needed.

    integer(ip), intent(in) :: n
    !! The order of the element.
    real(wp), dimension(2), intent(in) :: vx
    !! The physical coordinates of the left and right boundary of the element.
    real(wp), dimension(n+1), intent(in) :: r
    !! The node locations inside the reference element, \(r_i\).
    real(wp), dimension(n+1,n+1), intent(in) :: dr
    !! The derivative matrix for the reference element, \(D_{ij}\).
    real(wp), dimension(n+1), intent(out) :: x
    !! On output the physical coordinates of the nodes of the element, \(x_i\).
    real(wp), dimension(n+1), intent(out) :: xr
    !! On output contains \((xr)_i=\sum_{i,j}^{n+1} D_{ij} x_j\).
    real(wp), dimension(n+1), intent(out) :: rx
    !! On output contains \((rx)_i=\frac{1}{(xr)_i}\).
    integer(ip) :: i,j

    do i = 1, n+1
      x(i) = vx(1)+0.5_wp*(1.0_wp+r(i))*(vx(2)-vx(1))
    end do
    xr = 0.0_wp
    do i = 1, n+1
      do j = 1, n+1
        xr(i) = xr(i) + dr(i,j) * x(j)
      end do
    end do
    rx = 1.0_wp / xr
    
  end subroutine Jacobian
  
end submodule grid_implementation
