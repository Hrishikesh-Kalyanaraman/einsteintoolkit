module all_integrators
!! Module that provides a single view of all integrators. The implementation
!! is provided in the module itself.

  use rk4_integrator
  use rk5_integrator
  use abmv5_integrator

  implicit none

  type(rk4), target :: rk4_int
  !! The RK4 integrator.
  type(rk5), target :: rk5_int
  !! The RK5 integrator.
  type(abmv5), target :: abmv5_int
  !! The ABMV5 integrator.
  class(integrator), pointer :: mol_int
  !! Pointer that can point to any kind of integrator.

  contains

    subroutine choose_integrator ( )
    !! Routine that chooses an integrator by pointing to one of the
    !! available ones.

      use parameters, only : mol_integrator

      implicit none

      select case ( mol_integrator )
      case ( 'rk4' )
        mol_int => rk4_int
      case ( 'rk5' )
        mol_int => rk5_int
      case ( 'abmv5' )
        mol_int => abmv5_int
      case default
        print*,'Invalid integrator requested'
        stop
      end select

    end subroutine choose_integrator


    function mol_ntmp ( ) result ( ntemp )
    !! Function that interfaces with the ntemp functions provided by the
    !! individual integrators. This function can be called before
    !! [[choose_integrator]] is called.

      use parameters, only : mol_integrator

      implicit none

      integer(ip) :: ntemp
      !! The number of temporary storage levels needed.

      select case ( mol_integrator )
      case ( 'rk4' )
        ntemp = rk4_int%ntemp ( )
      case ( 'rk5' )
        ntemp = rk5_int%ntemp ( )
      case ( 'abmv5' )
        ntemp = abmv5_int%ntemp ( )
      case default
        print*,'Invalid integrator requested'
        stop
      end select

    end function mol_ntmp

end module all_integrators
