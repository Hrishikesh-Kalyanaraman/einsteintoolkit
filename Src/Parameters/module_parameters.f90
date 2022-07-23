module parameters
!! Definition of all parameters.

  use kinds

  implicit none

  character(len=32) :: equation_name = 'scalar_schwarzschild'
  !! The name of the PDE equation to evolve.
  integer(ip) :: n_elems = 32
  !! The number of DG elements.
  integer(ip) :: order = 16
  !! The order of the DG elements.
  real(wp) :: Sminus = -20.0_wp
  !! The computational coordinate of the horizon.
  real(wp) :: r_center = 10.0_wp
  !! The Schwarzschild coordinate of the center of the domain. Not used if
  !! use_particle = .true.
  integer(ip) :: t_size = 5
  !! The radius of the region where time-dependent coordinates are used (in
  !! units if DG elements).
  logical :: use_particle = .false.
  !! If .true. use the effective source, if .false. setup an initial Gaussian
  !! initial data profile.
  real(wp) :: mass = 1.0_wp
  !! The mass of the black hole.
  integer(ip) :: lmin = 0
  !! The minimum l-value to evolve.
  integer(ip) :: lmax = 20
  !! The maximim l-value to evolve.
  real(wp) :: gaussian_center = 0.0_wp
  !! The tortoise coordinate of the center of the Gaussian initial data profile.
  real(wp) :: sigma = 1.0_wp
  !! The width of the Gaussian initial data profile.
  real(wp) :: amplitude = 1.0_wp
  !! The amplitude of the Gaussian initial data profile.
  real(wp) :: t_initial = 0.0_wp
  !! The initial time value.
  real(wp) :: t_final = 10.0_wp
  !! The final time value.
  real(wp) :: p_orb = 10.0_wp
  !! The semilatus rectum of the initial geodesic orbit.
  real(wp) :: ecc = 0.0_wp
  !! The eccentricity of the initial geodesic orbit.
  real(wp) :: phi_initial = 0.0_wp
  !! The initial value of the \(\phi\) coordinate (used to match externally
  !! provided initial data).
  real(wp) :: coufac = 0.5_wp
  !! The courant factor.
  logical :: use_world_tube = .false.
  !! If .true. use a world_tube around the particle. Should always be set when
  !! using a particle.
  integer(ip) :: world_tube_width = 1
  !! The width of the world_tube in units of the DG element.
  logical :: turn_on_source_smoothly = .false.
  !! If .true. turn the effective source on smoothly for all element. When
  !! .false. at least one mode requires external initial data.
  real(wp) :: tsigma = 10.0_wp
  !! The width of the 'torder'th order Gaussian used to turn on the effective
  !! source smoothly.
  integer(ip) :: torder = 4
  !! The order of the Gaussian used to turn on the effective source smoothly.
  integer(ip) :: out0d_every = -1
  !! Zero dimensional output frequency. -1 means no output.
  integer(ip) :: out1d_every = -1
  !! One dimensional output frequency. -1 means no output.
  real(wp) :: q_charge = 1.0_wp
  !! The charge \(q\) of the scalar charge.
  real(wp) :: q_mass = 1.0_wp
  !! The initial mass \(m\) of the scalar charge.
  logical :: use_osculating_orbit = .true.
  !! If .true., the osculating orbits framework is used for orbit evolution. If
  !! .false., the geodesic evolution equations are used for orbit evolution.
  logical :: use_generic_orbit = .false.
  !! If .true., time-dependent coordinates are used. If .false., the orbit has
  !! to be circular.
  logical :: evolve_orbit = .false.
  !! If .true., apply the self-force back-reaction. If .false., a geodesic
  !! orbit is maintained.
  logical :: turn_on_force_smoothly = .true.
  !! If .true., turn on the back-reaction force smoothly. If .false., switch
  !! to full force instantaneously.
  logical :: use_chi = .false.
  !! If .true., use \(\chi\) for turning on the source. If .false., use time.
  real(wp) :: evolve_after = 0.0_wp
  !! After this value of (\chi\) or time turn on the force.
  real(wp) :: force_sigma = acos(-1.0_wp)
  !! The width of the smooth turn on of the force.
  logical :: fit_high_l = .false.
  !! If .true., use some higher l-modes to fit the falloff and use that to
  !! improve the total self-force. If .false., use the sum of the evolved modes
  !! for the total self-force.
  integer(ip) :: first_l = 12
  !! Smallest l-mode to include in the fit.
  integer(ip) :: last_l = 20
  !! Largest l-mode to include in the fit.
  integer(ip) :: fit_order = 2
  !! The order of the lowest order fit term.
  integer(ip) :: nfit = 1
  !! The total number of terms to include in the fit.
  logical :: output_coords_for_exact = .false.
  !! When .true., the coordinates are output to file.
  logical :: use_exact_initial_data = .false.
  !! If .true., use external initial data for some modes.
  integer(ip) :: exact_initial_data_lmax = -1
  !! Use external initial data for all modes from lmin to this value.
  character(len=1024) :: input_directory
  !! Directory containing the external initial data.
  character(len=256) :: input_basename
  !! The base name for the input files.
  character(len=5) :: mol_integrator = 'rk4'
  !! The integrator to use. Currently allowed values are 'rk4', 'rk5' and
  !! 'abmv5'.
  logical :: use_constant_acceleration = .false.
  !! Use a constant accelerated geodesic orbit (either circular or eccentric).
  logical :: use_gaussian_acceleration = .false.
  !! Use a Gaussian acceleration profile (only for circuler orbits).
  real(wp) :: accel_amp = 0.0
  !! Amplitude of the Gaussian acceleration profile.
  real(wp) :: accel_sigma = 1.0
  !! Width (in time) of the Gaussian acceleration profile.
  real(wp) :: accel_t0 = 0.0
  !! The time at which the acceleration event is centered.
  real(wp) :: omega_ratio = 1.0
  !! For constant accelerated orbits, the ratio of the accelerated to geodesic
  !! angular velocity.
  logical :: use_filter = .false.
  !! Should a filter be applied to the evolved variables.
  integer :: filter_order
  !! The order of the filter.
  integer :: filter_cutoff
  !! The filter cutoff below which the low modes are left untouched.
  logical :: use_analytic_driver = .false.
  !! Should an analytic prescription for the acceleration be used to drive the
  !! evolution.
  integer(ip) :: a_timelevels = -1
  !! Number of timelevels to use for smooth acceleration derivatives when
  !! Runge-Kutta time integrators are used. A value of -1 means that these
  !! are disabled.
  integer(ip) :: a_smooth_derivative_order = 3
  !! The order of the smooth acceleration derivatives. The value should at most
  !! be [[a_timelevels]]-1.
  logical :: use_hermite_extrapolation = .false.
  !! Should hermite interpolation be used?
  integer(ip) :: hermite_order = 3
  !! The order of the hermite interpolation. Currently allowed values are 3 and
  !! 5.
  logical :: use_smooth_derivs_for_hermite = .false.
  !! Should smoothing derivatives be used to estimate the derivative values
  !! for Hermite interpolation.
  logical :: use_adot = .false.
  !! Use the first time derivative of the 4-acceleration in the effective
  !! source.
  logical :: use_addot = .false.
  !! Use the second time derivative of the 4-acceleration in the effective
  !! source.
  logical :: use_field_observer = .false.
  !! Use an observer to extract field data at the horizon and scri+.

  namelist /params/ equation_name, n_elems, order, Sminus, r_center, t_size, &
                    use_particle, mass, lmin, lmax, gaussian_center, &
                    sigma, amplitude, t_initial, t_final, p_orb, ecc, &
                    phi_initial, coufac, use_world_tube, world_tube_width, &
                    turn_on_source_smoothly, tsigma, torder, &
                    out0d_every, out1d_every, q_charge, q_mass, &
                    use_osculating_orbit, use_generic_orbit, evolve_orbit, &
                    turn_on_force_smoothly, use_chi, evolve_after, &
                    force_sigma, fit_high_l, first_l, last_l, fit_order, &
                    nfit, output_coords_for_exact, use_exact_initial_data, &
                    exact_initial_data_lmax, input_directory, input_basename, &
                    mol_integrator, use_constant_acceleration, &
                    use_gaussian_acceleration, accel_amp, accel_sigma, &
                    accel_t0, omega_ratio, use_filter, filter_order, &
                    filter_cutoff, use_analytic_driver, a_timelevels, &
                    a_smooth_derivative_order, use_hermite_extrapolation, &
                    hermite_order, use_smooth_derivs_for_hermite, use_adot, &
                    use_addot, use_field_observer

  contains
    subroutine read_parameters ()
    !! Read in the run-time parameters from a file. The name of the file
    !! is read from the first command line option. If no command line options
    !! are given, use 'input.par'.
    !!
    !! The parameters are read in one go using the namelist mechanism.

      implicit none
      integer :: argc, arglen, errno
      character(len=256) :: inputpar

      argc = command_argument_count()
      if (argc == 0) then
        inputpar = 'input.par'
      else
        call get_command_argument(1, inputpar, arglen, errno)
        if (errno .ne. 0) then
          write(*,*) "Error retrieving command line argument"
          stop
        endif
      endif

      open(1, file = inputpar, status='old', form = 'formatted', &
              action = 'read' )
      read(1, nml=params)
      close(1)
    end subroutine read_parameters

end module parameters
