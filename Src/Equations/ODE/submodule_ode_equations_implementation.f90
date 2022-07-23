submodule(ode_equations) ode_equations_implementation
!! The implementation of the interfaces defined in [[ode_equations]].

  implicit none

contains

  module procedure ode_set_to_zero

    implicit none

    if (dest<-1 .or. dest>this%ntmp) then
      print*,'Error: geod_schw_set_to_zero called with invalid destination'
      stop
    end if

    if ( dest == -1 ) then
      this%rhs_data = rzero
    else if ( dest == 0 ) then
      this%var_data = rzero
    else
      this%tmp_data(:,dest) = rzero
    end if

  end procedure ode_set_to_zero


  module procedure ode_update_vars

    implicit none

    integer(ip) :: i, ntmp
    real(wp) :: alpha, beta

!    print*,'ODE equation update vars called'
!    print*,'Data = ', this%var_data
!    print*,'RHS = ', this%rhs_data
    ntmp = this%ntmp

    if (source<=-2 .or. source>ntmp) then
      print*,'Error: update_vars called with incorrect source argument'
      stop
    end if
    if (dest<=-1 .or. dest>ntmp) then
      print*,'Error: update_vars called with incorrect dest argument'
      stop
    end if
    if ( present(source2) .and. ( source2<=-2 .or. source2>ntmp) ) then
      print*,'Error: update_vars called with incorrect source2 argument'
      stop
    end if
    if ( present(source2) .and. (source2==dest) ) then
      print*,'Error: when source2 is present in update_vars it has to be different than dest'
      stop
    end if

    if ( present(scalar) ) then
      alpha = scalar
    else
      alpha = 1.0_wp
    end if

    if ( present(scalar2) ) then
      beta = scalar2
    else
      beta = 1.0_wp
    end if

    if ( dest==0 ) then
      select case (source)
      case (-1)
        this%var_data = alpha*this%rhs_data
      case (0)
        this%var_data = alpha*this%var_data
      case (1:)
        this%var_data = alpha*this%tmp_data(:,source)
      end select
      if ( present(source2) ) then
        select case (source2)
        case (-1)
          this%var_data = this%var_data+beta*this%rhs_data
        case (0)
          this%var_data = this%var_data+beta*this%var_data
        case (1:)
          this%var_data = this%var_data+beta*this%tmp_data(:,source2)
        end select
      end if
    end if
    if ( dest>0 ) then
      select case (source)
      case (-1)
        this%tmp_data(:,dest) = alpha*this%rhs_data
      case (0)
        this%tmp_data(:,dest) = alpha*this%var_data
      case (1:)
        this%tmp_data(:,dest) = alpha*this%tmp_data(:,source)
      end select
      if ( present(source2) ) then
        select case (source2)
        case (-1)
          this%tmp_data(:,dest) = this%tmp_data(:,dest)+beta*this%rhs_data
        case (0)
          this%tmp_data(:,dest) = this%tmp_data(:,dest)+beta*this%var_data
        case (1:)
          this%tmp_data(:,dest) = this%tmp_data(:,dest) &
                                  +beta*this%tmp_data(:,source2)
        end select
      end if
    end if
!    print*,'Data = ', this%var_data
!    print*,'RHS = ', this%rhs_data
!    print*,'ODE equation update vars exited'
  end procedure ode_update_vars


  module procedure ode_print_data

    implicit none

    integer(ip) :: j

    print*,'ode data = ', this%var_data(1:3)
    print*,'rhs data = ', this%rhs_data(1:3)
    do j = 1, this%ntmp
      print*,'temp data(', j, ') = ', this%tmp_data(1:3,j)
    end do
  end procedure ode_print_data

end submodule ode_equations_implementation
