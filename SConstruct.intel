import sys
import os
env = DefaultEnvironment(ENV = {'PATH' : os.environ['PATH'],
                                'LD_LIBRARY_PATH' : os.getenv('LD_LIBRARY_PATH',default='')},
                         LINK = 'ifort',
                         LINKFLAGS = '-g -mkl',
                         TOOLS= ['default', 'ifort'])

# pretty output
if ARGUMENTS.get('VERBOSE') != '1':
    if sys.stdout.isatty():
        try:
            from termcolor import colored
        except ImportError:
            print("The module termcolor could not be imported so turning colors off.")
            print("To activate colors, either get an administrator to install termcolor or install it in user space by:\n")
            print("[pip|pip3] install --user termcolor\n")
            print("where the appropriate pip should be used depending on your python version")
            def colored(x,_):
                return x
    else:
        def colored(x,_):
            return x

    env['CXXCOMSTR'] = colored("Compiling", "green") + " $TARGET"
    env['F90COMSTR'] = colored("Compiling", "green") + " $TARGET"
    env['FORTRANCOMSTR'] = colored("Compiling", "green") + " $TARGET"
    env['LINKCOMSTR'] = colored("Linking", "blue") + " $TARGET"
    env['HDF5COMSTR'] = colored("Generating", "magenta") +  " $TARGET"

# Build options
env['LIBS'] = ['stdc++', 'gsl', 'm', 'gslcblas']
env['LIBPATH']  = ['/usr/local/packages/gsl-2.4-gnu/lib']
env['F90FILESUFFIXES']=['.f90','.f']
env['F90']      = 'ifort'
env['FORTRAN']  = 'ifort'
env['F90FLAGS'] = ['-O3', '-g', '-qopenmp',  '-xHOST', '-qopt-report', '-qopt-report-phase=vec,openmp', '-mkl', '-no-wrap-margin', '-diag-disable=10397']
env['FORTRANFLAGS'] = ['-O3', '-g', '-qopenmp', '-xHOST', '-qopt-report', '-qopt-report-phase=vec,openmp', '-mkl', '-no-wrap-margin', '-diag-disable=10397']
#env['FORTRANMODDIRPREFIX'] = '-module '
env['FORTRANMODDIR'] = '${TARGET.dir}'
env['CPPPATH']  = ['/usr/local/packages/gsl-2.4-gnu/include', 'EffSource/ScalarSchwarzschild/EffectiveSource']
env['CXXFLAGS'] = ['-O3', '-g', '-qopenmp', '-std=c++11', '-xHOST', '-qopt-report', '-qopt-report-phase=vec,openmp', '-Wall', '-diag-disable=10397']
env['CXX'] = 'icpc'

Export('env')

SConscript('Src/SConscript', exports='env', variant_dir='Build', duplicate=0)

Decider('MD5-timestamp')
