# SelfForce-1D

SelfForce-1D is a code infrastructure for simulating Extreme Mass Ratio Inspirals using the effective source approach to the self-force problem. Currently, only a scalar charge in a Schwarzschild spacetime background has been implemented, but it is the hope that more systems will be added soon.

This is a complete rewrite of 1DScalarWaveDG. The new code is much more modular in order to make it easy to extend the code with other systems of equations.

## Requirements

The code is written in modern object oriented fortran and, hence, will not work with older compilers. The code has been successfully compiled with gfortran 7.3.0 and ifort 17.0.7. It may also work with some compilers older than that, but that is untested.

SCons is used for building (instead of make) and version 3.1.0 is required as previous versions contained a bug that prevented correct dependency determination when submodules and type bound procedures are used. A local version of SCons 3.1.1 is included for people who choose not to install a new enough version of SCons.

In terms of external libraries, BLAS/LAPACK and GSL are needed.

The C++ code that provides the effective source for a scalar charge in a Schwarzschild spacetime (by Barry Wardell) is now included in the checkout (before it had to be obtained separately). In order to extract the Wigner D include files, the `h5dump` command from `HDF5` is now needed.

## Installation

Checkout the code from bitbucket using: `git clone https://peterdiener@bitbucket.org/peterdiener/selfforce-1d.git SelfForce-1D`

Change to the directory: `cd SelfForce-1D`

Modify the `SConstruct` file appropriately to match your compilers and libraries. Example `SConstruct` files for the intel (`SConstruct.intel`) and gnu (`SConstruct.gnu`) compilers are provided. It should be fairly simple to figure out what changes needs to be made.

If there is an appropriate version of `SCons` on your system and `scons` is in your path, you can now build the code using: `scons`

In case a new enough version of SCons is not available and you choose not to install it yourself, the SCons-local version provided in the `SCons-local` directory can be used. To do this, the provided tarball needs to be extracted. In order not to clutter the source tree with untracked files, we suggest creating a directory elsewhere and copy the tarball there before extracting it, e.g. to create a `SCons-local-3.1.1` directory in your home directory do:

`mkdir ~/SCons-local-3.1.1`

`cp SCons-local/scons-local-3.1.1.tar.gz ~/SCons-local-3.1.1`

`cd ~/SCons-local-3.1.1`

`tar xvzf scons-local-3.1.1.tar.gz`

and then change back to the `SelfForce-1D` directory. Then build the code using: `python ~/SCons-local-3.1.1/scons.py`

During first time compilation you may see warnings about include directories not existing. These can be safely ignored as the directories are created during compilation.

Upon successful completion, a `Build` directory will have been created and inside that, directories for each of the directories in the `Src` directory. These contain the object and module files. Finally, an `Exe` directory will have been created. This directory contains three exectubles of which two (`accel_history_test.x` and `accel_test.x`) are simpler test codes. Despite it's name (and this should change soon) `test.x` is the main executable.

Due to an issue with scons, the first build will fail if it is done in parallel (with `scons -jN` where `N>1`). However subsequent builds can be done in parallel as long as the `Build` directory is not removed.

## Running the code

The code uses the Fortran namelist I/O feature to read in it's parameters from parameter files. A parameter file starts with `&params` (params is the internal name for the namelist) and ends with an `/`. In between, values can be assigned to the variables in the namelist. Everything following a `!` is interpreted as a comment. All valid parameters are defined in `Src/Parameters/module_parameters.f90`.

### Evolution of Gaussian initial data

An example parameter file is provided in `Par/scalargaussian_o8.par`. This parameter file sets up a Gaussian initial data pulse of a scalar field in a Schwarzschild spacetime and evolves it for 1000 M.

Create a directory somewhere (preferrably outside of the `SelfForce-1D` directory) and copy the example parameter file into it. The code is parallelized using OpenMP, so you may want to set the `OMP_NUM_THREADS` environment to something suitable before running. In this parameter file, only 4 modes are evolved. Since the parallelization is done over modes, I would recommend to at most using `OMP_NUM_THREADS=4`. Change to the directory and launch the executable:

`<path to SelfForce1D>/Build/Exe/test.x scalargaussian_o8.par`

To visualize some of the data produced in the run, copy the `Par/plot_tails_ScalarGaussianOrder8.gp` file to the directory that contains the data, start up `gnuplot` and, at the gnuplot command prompt, execute:

`load 'plot_tails_ScalarGaussianOrder8.gp'`

This gnuplot script plots data from the files `psi.[123].extract.asc` that contains the scalar field extracted at the horizon and future null infinity (Scri+) as a function of coordinate time in a log-log plot for three different modes (l=0,1 and 2). What you should see is the expected power law tail behavior. At the horizon, the exponent should be -(2l+3), while at Scri+ it should be -(l+2). The dotted lines indicates pure power law decay with exponents (-2,-3,-4 and -5). It is clear that the expected power law decay is seen in all cases except for l=2 at the horizon, where roundoff error is reached before the tail behavior sets in.

Another example visualization script is `Par/plot_qnm_ScalarGaussianOrder8.gp`. Copy that as well and execute:

`load 'plot_qnm_ScalarGaussianOrder8.gp'`

Here the early part of the waveform at the horizon is plotted in a log-plot and the exptected quasi-normal mode (QNM) oscillations (using Emanuele Berti's data from <https://pages.jh.edu/~eberti2/ringdown/>) are plotted on top. For l=0, there are not quite enough oscillations before the tail sets in to see a clear agreement. On the other hand, for l=1 and l=2 the agreement is excellent after the initial transient have died down and the tail (or roundoff error) sets in.

### Evolution of a particle on an eccentric geodesic

An example parameter file that use a particle is provided in `Par/scalar_p7.2_e0.5_o8.par`. This parameter file sets up a particle moving on an eccentric geodesic with semi-latus rectum p=7.2 and eccentricity e=0.5. It uses a world-tube of default size 1. The osculating orbit description is used and since this is a generic orbit (as opposed to a circular orbit with constant radius) the time dependent coordinate transformation is used in a region of 10 DG elements on either side of the partcle location. The effective source is turned on smoothly (starting from zero initial data) but quickly (timescale of 0.1 M). All modes from l=0 to l=10 are evolved.

As in the previous example, create a directory for this example, copy the parameter file to that directory and change to it. Then execute:

`<path to SelfForce1D>/Build/Exe/test.x scalar_p7.2_e0.5_o8.par`

This run will take a bit longer, since more modes are evolved, but will also be able to take advantage of more cores if they are available. To look at some data, copy `Par/plot_sf_ScalarP7.2E0.5_Order8.gp` to the directory, start up gnuplot, and execute:

`load 'plot_sf_ScalarP7.2E0.5_Order8.gp'`

You will see a timeseries plot of the radial component of the self-force extracted at the particle location. Note how there is an initial transient followed by periodic behavior (the radial period of the orbit). Also note that the amplitude of the modes (for l>4) decreases with l, as should be expected for a convergent mode-sum.

### Initial data

When starting an evolution with zero initial data and a smooth activation of the effective source, some initial transients are generated that need to propagate away before the correct self-force for the given orbit can be extracted. Due to the slow decay of the tails for low l-modes, it may be necessary to evolve for quite a long time before sufficient accuracy in the self-force is reached. If initial data consistent with the orbit is available, it is possible to start the evolution with the effective source turned on at full strength and the correct self-force at the location of the particle from the beginning. Niels Warburton has provided a python code (available in the `InitialData` directory) using computations in the frequency domain that is able to provide initial data for the full retarded field for both geodesic and accelerated eccentric orbits for any given set of coordinates. Before using this code, I suggest copying `InitialData` and it's subdirectories to a location outside of the `SelfForce-1D` directory tree in order to avoid cluttering it with files that are not tracked in the git repository. Generating your own initial data and using them is a multi-step process.

#### Step 1: obtaining a set of coordinates where the initial data should be provided

Write an input parameter file for the required orbit and set the parameter `output_coords_for_exact = .true.`. Example parameter files are provided in `Par/scalar_p8.0_e0.1_o<n>.par` where `n=8,16,32`. They generate files with the coordinates for the case of semilatus rectum `p=8` and eccentricity `e=0.1` for DG order `n=8,16, 32`. Note that `t_final = 0.0` as these files are only used to generate a file with the computational coordinates. After running any of the parameter files, a file named `coords.asc` should be created.

#### Step 2: running the python code

Copy the `coords.asc` file generated at the end of step 1 to the `InitialData/data/input` directory with the name based on the template `coords_p<p>_e<e>_n<n>_accel.dat`, where `<p>` should be replaced with the value for the semilatus rectum `p`, `<e>` by the value of the eccentricity `e` and `<n>` by the DG order `n`. As can be seen there are already example files provided for the 3 cases mentioned in step 1. In this particular case, execute:

`python SSF_init.py 8 0.1 1.0 8 0 2`

in the `InitialData` directory in order to run the python code. In general the usage is:

`python SSF_init.py p e kappa n lmin lmax`

where `p` is the semilatus rectum, `e` is the eccentricity, `kappa` indicates how accelerated the orbit is, `n` is the DG order used and `lmin` and `lmax` are the minimum and maximum l-values to calculate initial data for. In the examples provided, `p=8`, `e=0.1`, `kappa=1` (geodesic orbit), `n=8,16,32`, `lmin=0` and `lmax=2`. While running the python code it will print a line starting with `('phi value at chi=pi: '`. Note down the numerical value following that. That is the initial value for the azimuthal angle that should be used when starting the evolution in step 3. After running, initial data files are created in `InitialData/data/output` with filenames according to the template `SSF_init_data_p<p>_e<e>_n<n>_l<l>m<m>.dat`. Example initial data files are already provided for the 3 cases mentioned in step 1.

#### Step 3: Running the self-force code with the generated initial data.

The parameter file used to generate the coordinate file in step 1 now needs to be modified in order to be able to read in the initial data files produced in step 2. First, change the value of `t_final`, then change the value of `turn_on_source_smoothly` from `.true.` to `.false.`. In addition, values for 4 parameters have to be provided. These are `use_exact_initial_data = .true.`, `exact_initial_data_lmax`, `input_directory`, and `input_basename`. When using initial data files, the value of `use_exact_initial_data` should always be `.true.`. The parameter `exact_initial_data_lmax` specifies up to what value of l initial data should be used and obviously depends on the situation. If the value of `exact_initial_data_lmax` is less than `lmax`, the effective source is still turned on smoothly for the modes with `exact_initial_data_lmax<l<=lmax`. The value of `input_directory` is a string that specifies in which directory (either full path or path relative to the run directory) the initial data files are located. The value of `input_basename` specifies the base part of the input data name, excluding the last part that specifies the `l` and `m` values for a given mode. The parameter files `Par/scalar_p8.0_e0.1_o<n>_evolve.par` are the modified versions of `Par/scalar_p8.0_e0.1_o<n>.par` that will read in the example input files mentioned in step 2--in other words, the relative path to the initial data files with respect to the parameter files inside the `SelfForce-1D` directory tree is used.

## Running tests

Regression and portability testing is supported by a python script `Src/Test/run_tests.py` and a set of test parameter files locted in the `Test` directory. To run the test, simply execute:

`Src/Test/run_tests.py`

in the main `SelfForce-1D` directory. The script will provide information about passing and/or failing tests when run. More careful examination can be performed afterwards by comparing the newly genereated test output in the `RunTests` directory with the original data in the `Test` directory.

## Documentation

Documentation of the interfaces of all classes and type bound procedures defined and used in the code were produced using FORD and is currently available at <https://www.cct.lsu.edu/~diener/SelfForce1D/Doc/index.html>.
