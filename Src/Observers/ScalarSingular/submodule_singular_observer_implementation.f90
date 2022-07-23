submodule(singular_observer) singuler_observer_implementation

  implicit none

contains

  module procedure sobs_init

    use parameters, only: lmin, lmax
    use grid, only: particle_element, particle_node, rho_particle
    use effective_source, only: eff_source

    implicit none

    integer :: nl, nvars
    select type(object)
    class is (eff_source)
      this%p => object
    class default
      print*,'Singular observer initialized with invalid type'
      stop
    end select


    nl = lmax - lmin + 1
    nvars = this%p%nvars
    this%nl = nl
    allocate ( this%vname, source = 'singular_observer' )
    allocate ( this%psi(nvars,nl), this%dpsidt(nvars,nl,2), &
               this%dpsidr(nvars,nl,2) )
    this%nradii = 2
    this%ioo_id = -1
    allocate ( this%elem_index(2), this%node_index(2), this%radii(2) )
    this%elem_index = particle_element
    this%node_index = particle_node
    this%radii = rho_particle

  end procedure sobs_init


  module procedure sobs_extract

    implicit none

  end procedure sobs_extract


  module procedure sobs_output

    use output_base
    use time_info

    implicit none

    character(len=16) :: filename = 'psi_singular.asc'
    integer(ip) :: ioo_id
    character(len=92) :: key = &
    '#1: time, 2,3: psi-, 4,5: dpsidtr-, 6,7: dpsidr-, 8,9: psi+, 10,11: dpsidt+, 12,13: dpsidr'

    ioo_id = this%ioo_id
    if ( ioo_id<0) then
      ioo_id = next_available_io_id ()
      this%ioo_id = ioo_id
      print*,'Opening ', filename, ' with id ', ioo_id
      open(ioo_id, file=filename, status='replace', action='write')
      write(ioo_id,*) key
    end if

!    write(ioo_id,'(*(es32.15e3,1x))') get_current_time ( ), &
!                                     real(this%psi(1),wp), &
!                                     aimag(this%psi(1)), &
!                                     real(this%dpsidt(1),wp), &
!                                     aimag(this%dpsidt(1)), &
!                                     real(this%dpsidr(1),wp), &
!                                     aimag(this%dpsidr(1)), &
!                                     real(this%psi(2),wp), &
!                                     aimag(this%psi(2)), &
!                                     real(this%dpsidt(2),wp), &
!                                     aimag(this%dpsidt(2)), &
!                                     real(this%dpsidr(2),wp), &
!                                     aimag(this%dpsidr(2))
  end procedure sobs_output


  module procedure close_sing_observer

    implicit none

    deallocate ( this%psi, this%dpsidt, this%dpsidr )
    nullify ( this%p )

  end procedure close_sing_observer

end submodule singuler_observer_implementation
