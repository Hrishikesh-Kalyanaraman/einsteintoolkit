submodule(self_force_observer) self_force_observer_implementation
!! The implementation of the interface in [[self_force_observer]].

  implicit none

contains

  module procedure sf_init

    use parameters, only : lmin, lmax, fit_high_l, first_l, last_l
    implicit none

    integer(ip) :: n_unique
    select type(object)
    type is (scal_schw)
      this%p => object
    class default
      print*,'Self-force observer initialized with invalid type'
      stop
    end select

    allocate ( this%vname, source = this%p%ename )

    if (size(rad)/=1) then
      print*,'Self-force observer should be initialized with 1 extraction radius'
      stop
    end if
    this%nradii = size(rad)
    this%radii = rad
    this%ioo_id = -1

    call find_indices ( rad, coord, this%elem_index, this%node_index )

    n_unique = n_unique_values ( this%p%ll )
    if (.not. fit_high_l) then
      allocate ( this%fl(lmin:lmax), this%ftl(lmin:lmax), &
                 this%fphil(lmin:lmax), this%frl(lmin:lmax) )
    else
      allocate ( this%fl(lmin:lmax+3), this%ftl(lmin:lmax+3), &
                 this%fphil(lmin:lmax+3), this%frl(lmin:lmax+3) )
    end if
    if ( fit_high_l ) then
      if ( lmin>0 ) then
        print*,'fit_high_l should only be .true. when lmin=0'
        stop
      end if
      if ( first_l < lmin ) then
        print*,'Error: first_l < lmin'
        stop
      end if
      if ( last_l > lmax ) then
        print*,'Error: last_l < lmax'
        stop
      end if
    end if

  end procedure sf_init


  module procedure sf_extract

    use orbit_base
    use gsl_interface, only : legendre_sphPlm
    use parameters, only : mass, use_generic_orbit, lmax, fit_high_l, &
                           first_l, last_l, fit_order, nfit, use_particle, &
                           q_charge
    use numerics, only : correct_for_higher_modes
    use time_info, only : get_current_time
    use world_tube, only : wt

    implicit none

    real(wp) :: r, phi, ur, En, Lz, m_fold_factor, y_lm, ftl_new, frl_new
    real(wp) :: flm, ftlm, ftlml, ftlmr, frlm, frlml, frlmr, fphilm
    complex(wp) :: phifactor, phase, psip, psitpr, psitpl, psirpr, psirpl
    complex(wp) :: psit_new, psir_new
    integer(ip) :: i, j, k
    complex(wp), dimension(1) :: psis, dpsisdrl, dpsisdrr, &
                                 dpsisdtl, dpsisdtr, dpsisdphi

    i = this%elem_index(1)
    j = this%node_index(1)
    this%fl = 0.0_wp
    this%ftl = 0.0_wp
    this%fphil = 0.0_wp
    this%frl = 0.0_wp

    call orbit_info%get_orbit ( r, phi, ur, En, Lz )    

    associate ( nmodes => this%p%nmodes, &
                mm => this%p%mm, &
                ll => this%p%ll, &
                var => this%p%eq_data, &
                fl => this%fl, &
                ftl => this%ftl, &
                fphil => this%fphil, &
                frl => this%frl )
      do k = 1, nmodes
        if ( mm(k) == 0 ) then
          m_fold_factor = 1.0_wp
        else
          m_fold_factor = 2.0_wp
        end if

        phifactor = zi*mm(k)
        phase = cmplx ( cos(mm(k)*phi), sin(mm(k)*phi), wp )
        y_lm = legendre_sphPlm ( ll(k), mm(k), 0.0_wp )
        psip = 0.5_wp * ( var(1,k)%elems(i)%var(j) + &
                          var(1,k)%elems(i+1)%var(1) )
        psitpl = var(2,k)%elems(i)%var(j)
        psitpr = var(2,k)%elems(i+1)%var(1)
        psirpl = var(3,k)%elems(i)%var(j)
        psirpr = var(3,k)%elems(i+1)%var(1)

        if (use_generic_orbit) then
          call this%p%time_dep_coord%tdc_to_tortoise ( i, 2, psitpl, psirpl, &
                                                       psit_new, psir_new )
          psitpl = psit_new
          psirpl = psir_new
          call this%p%time_dep_coord%tdc_to_tortoise ( i, 2, psitpr, psirpr, &
                                                       psit_new, psir_new )
          psitpr = psit_new
          psirpr = psir_new
        end if

        if ( present(effs) ) then
          call effs%get_singular ( r, k, psis )         
          call effs%get_dsingular_dt ( r, -1, k, dpsisdtl )         
          call effs%get_dsingular_dt ( r, 1, k, dpsisdtr )         
          call effs%get_dsingular_dr ( r, -1, k, dpsisdrl )         
          call effs%get_dsingular_dr ( r, 1, k, dpsisdrr )         

          psip = psip - q_charge*psis(1)
          psitpl = psitpl - q_charge*dpsisdtl(1)
          psitpr = psitpr - q_charge*dpsisdtr(1)
          psirpl = psirpl - q_charge*dpsisdrl(1)
          psirpr = psirpr - q_charge*dpsisdrr(1)
        end if

        flm = m_fold_factor*y_lm*real(phase*psip,wp)
        ftlml = m_fold_factor*y_lm*real(phase*psitpl,wp)
        ftlmr = m_fold_factor*y_lm*real(phase*psitpr,wp)
        frlml = m_fold_factor*y_lm*real(phase*psirpl,wp)
        frlmr = m_fold_factor*y_lm*real(phase*psirpr,wp)
        fphilm = m_fold_factor*y_lm*real(phifactor*phase*psip,wp)

        ftlm = 0.5_wp * ( ftlml + ftlmr )
        frlm = 0.5_wp * ( frlml + frlmr )

        fl(ll(k)) = fl(ll(k)) + flm
        ftl(ll(k)) = ftl(ll(k)) + ftlm
        frl(ll(k)) = frl(ll(k)) + frlm
        fphil(ll(k)) = fphil(ll(k)) + fphilm
      end do 

      frl = frl/(r-2.0_wp*mass)-fl/r**2
      fl = fl/r
      ftl = ftl/r
      fphil = fphil/r

      if (fit_high_l) then
        fl(lmax+1:lmax+3) = correct_for_higher_modes ( fl(0:lmax), &
                                    lmax+1, first_l, last_l, fit_order, nfit )
        ftl(lmax+1:lmax+3) = correct_for_higher_modes ( ftl(0:lmax), &
                                    lmax+1, first_l, last_l, fit_order, nfit )
        frl(lmax+1:lmax+3) = correct_for_higher_modes ( frl(0:lmax), &
                                    lmax+1, first_l, last_l, fit_order, nfit )
        fphil(lmax+1:lmax+3) = correct_for_higher_modes ( fphil(0:lmax), &
                                    lmax+1, first_l, last_l, fit_order, nfit )
      end if
    end associate

  end procedure sf_extract


  module procedure sf_output

    use output_base
    use time_info

    implicit none

    integer(ip) :: ioo_id, tmp_id
    character(len=len_trim(this%vname)+10) :: filename

    ioo_id = this%ioo_id
    if (ioo_id<0) then
      filename = trim(this%vname) // '.fl.asc'
      ioo_id = next_available_io_id ()
      this%ioo_id = ioo_id
      print*,'Opening ', filename, ' with id ', ioo_id
      open(ioo_id, file=filename, status='replace', action='write')
      filename = trim(this%vname) // '.ftl.asc'
      tmp_id = next_available_io_id ()
      print*,'Opening ', filename, ' with id ', tmp_id
      open(tmp_id, file=filename, status='replace', action='write')
      filename = trim(this%vname) // '.fphil.asc'
      tmp_id = next_available_io_id ()
      print*,'Opening ', filename, ' with id ', tmp_id
      open(tmp_id, file=filename, status='replace', action='write')
      filename = trim(this%vname) // '.frl.asc'
      tmp_id = next_available_io_id ()
      print*,'Opening ', filename, ' with id ', tmp_id
      open(tmp_id, file=filename, status='replace', action='write')
    end if
    write(ioo_id,'(*(es23.15e3,1x))') get_current_time ( ), this%fl(:) 
    write(ioo_id+1,'(*(es23.15e3,1x))') get_current_time ( ), this%ftl(:) 
    write(ioo_id+2,'(*(es23.15e3,1x))') get_current_time ( ), this%fphil(:) 
    write(ioo_id+3,'(*(es23.15e3,1x))') get_current_time ( ), this%frl(:) 
    
  end procedure sf_output


  module procedure close_sf_observer
  end procedure close_sf_observer


  function n_unique_values ( var ) result(n)
  !! Helper function that finds the number of unique values in a 1d integer
  !! array.
  !!
  !! This function does not seem to be used, so can probably safely be
  !! removed.

    use iso_c_binding

    integer(c_int), dimension(:), intent(in) :: var
    !! An 1d integer arrays.
    integer(ip) :: n
    !! The return value is the number of unique values in the input array.
    integer(c_int), dimension(size(var)) :: tmp
    integer(ip) :: i

    tmp(1) = var(1)
    n = 1

    do i=2, size(var)
      if ( .not. any ( tmp(1:n) == var(i) ) ) then
        n = n+1
        tmp(n) = var(i)
      end if
    end do
   
  end function n_unique_values

end submodule self_force_observer_implementation 
