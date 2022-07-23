import libSchw 
import boundary_conditions as bc
from math import *
import numpy as np
from scipy.integrate import ode
from scipy.special import sph_harm, lpn, lqn
import matplotlib.pyplot as plt
import sys, getopt, cmath, os

range1 = lambda start, end: range(start, end+1)


# Routine to merge and sort two arrays while keeping track of where the
# elements of the combined array are coming from. The routine returns
# a tuple containing the merged array and a mask with value False for
# elements that was originally in the first array and True for elements
# that was in the second array.
def merge_arrays(arr1, arr2):

    arr1m = np.zeros(len(arr1),dtype='bool')
    arr2m = np.ones(len(arr2),dtype='bool')

    wtype = np.dtype([('data',arr1.dtype),('mask',arr1m.dtype)])

    marr1 = np.empty(len(arr1),dtype=wtype)
    marr2 = np.empty(len(arr2),dtype=wtype)

    marr1['data'] = arr1
    marr2['data'] = arr2
    marr1['mask'] = arr1m
    marr2['mask'] = arr2m

    merged = np.concatenate([marr1,marr2])
    merged.sort(order='data')

    return merged['data'], merged['mask']


def main(argv):

    if len(sys.argv) != 7:
        print("Usage SSF_init.py p e kappa n lmin lmax, where:")
        print("\tp is the semi-latus rectum")
        print("\te is the orbital eccentricity")
        print("\tkappa the acceleration coefficient where dt/dchi = kappa*dt/dchi_geo, i.e., set to 1 for geodesic motion")
        print("\tn is the grid resolution")
        print("\tlmin is the minimum l-mode to compute")
        print("\tlmax is the maximum l-mode to compute")
        exit()

    p       = float(sys.argv[1])
    e       = float(sys.argv[2])
    kappa   = float(sys.argv[3])
    n_res   = float(sys.argv[4])
    lmin    = int(sys.argv[5])
    lmax    = int(sys.argv[6])

    schw = libSchw.Schw(p, e, kappa)

    # Output the orbital parmeters
    schw.print_orbit_parameters()

    # Output the phi value corresponding to the initial chi value.
    # This is needed to start the fortran code correctly.
    print("phi value at chi=pi: ", schw.phi_k[schw.N//2])

    # particle is initially at chi = pi
    chi0 = pi;
    
    # sets the minimum m-mode to compute - only change to non-zero for testing
    mmin = 0
    
    # sets the maximum value that n can reach (in case convergence isn't met). Maybe have a range of values for different eccentricities?
    nmax = 20
    
    # Load the grid data
    input_filename = 'data/input/coords_p%g_e%g_n%g_accel.dat' % (p, e, n_res)
    print("Reading coordinate data from %s" % input_filename)
    try:
        data = np.loadtxt(input_filename)
    except IOError:
        print("Coordinate data not found")
        exit()
    
    # format the various grids
    grid_rho    = data[:,0]
    grid_x      = data[:,1]                     # x is r_*
    grid_r      = data[:,2] + 2
    grid_t      = data[:,3] + schw.t(chi0)
   
    r0_index    = (np.abs(grid_r-schw.r_max)).argmin()          #FIXME this should be made generic to allow other starting values for chi

    # There are always 2 grid points at the particle location as it is always
    # at an element boundary. Make sure to use the index corresponding to the 
    # right most element.
    if grid_rho[r0_index] == grid_rho[r0_index+1]:
        r0_index += 1

    # test = bc.innerBCs(2, 2+1e-5, -22.41213529106035, 2*schw.Omega_phi)
 #    print test
 #    exit()
    
    rIn_index   = (np.abs(grid_r-2.01)).argmin()		    #FIXME this value currently seems to be too small leading to bad data for the radial derivative.
    rOut_index  = (np.abs(grid_r-1000.)).argmin()       #FIXME rOut needs to be in the wavezone for each mode so this splitting of the array will need to take place inside the lmn loop
    
    grid_rhoIn1 = grid_rho[0:rIn_index]
    grid_rhoIn2 = grid_rho[rIn_index:r0_index]
    
    grid_rIn1   = grid_r[0:rIn_index]
    grid_rIn2   = grid_r[rIn_index:r0_index]
    
    grid_tIn1   = grid_t[0:rIn_index]
    grid_tIn2   = grid_t[rIn_index:r0_index]
    
    # r0_index is the index of the first data point that needs the solution
    # exterior to the particle.
    grid_rhoOut = grid_rho[r0_index:rOut_index+1]
    grid_rOut   = grid_r[r0_index:rOut_index+1]
    grid_tOut   = grid_t[r0_index:rOut_index+1]
    
    grid_rOut2  = grid_r[rOut_index+1:-1]
    grid_tOut2  = grid_t[rOut_index+1:-1]
    grid_rhoOut2= grid_rho[rOut_index+1:-1]
        
    rIn         = grid_r[rIn_index]
    r0          = grid_r[r0_index]
    rOut        = grid_r[rOut_index]
    
    rkIn, maskIn   =  merge_arrays(schw.r_k[:schw.N//2+1], grid_rIn2)
    rkOut, maskOut =  merge_arrays(schw.r_k[:schw.N//2+1], grid_rOut)

    grid_rIn2_pos = np.nonzero(maskIn)[0]
    grid_rOut_pos = np.nonzero(maskOut)[0]
    
    lib_rkIn2_pos = np.nonzero(np.invert(maskIn))[0]
    lib_rkOut_pos = np.nonzero(np.invert(maskOut))[0]    
                        
    for l in range1(lmin, lmax):
        for m in range1(mmin, l):
            if (l+m)%2 == 1:
                continue
                
            converged = 0
            # arrays to store the results in
            fieldIn1  = np.zeros(len(grid_rIn1),dtype=complex)
            fieldIn2  = np.zeros(len(grid_rIn2),dtype=complex)
            fieldOut  = np.zeros(len(grid_rOut),dtype=complex)
            fieldOut2 = np.zeros(len(grid_rOut2),dtype=complex)
            
            field_rs_In1 = np.zeros(len(grid_rIn1),dtype=complex)
            field_rs_In2 = np.zeros(len(grid_rIn2),dtype=complex)
            field_rs_Out = np.zeros(len(grid_rOut),dtype=complex)
            field_rs_Out2 = np.zeros(len(grid_rOut2),dtype=complex)
            
            field_t_In1  = np.zeros(len(grid_rIn1),dtype=complex)
            field_t_In2  = np.zeros(len(grid_rIn2),dtype=complex)
            field_t_Out  = np.zeros(len(grid_rOut),dtype=complex)
            field_t_Out2 = np.zeros(len(grid_rOut2),dtype=complex)
            
            n = 0
            while converged != 1:
                print("Calculating (%d, %d, %d)-mode" % (l,m,n))
                omega = m*schw.Omega_phi + n*schw.Omega_r
                
                # compute the homogeneous solutions via numerically integrating from the boundary (for radiative modes) or analytically (for static modes)  
                if(omega !=0):         
                    RkIn  = solIn(l, omega, rIn, r0, rkIn)
                    RkOut = solOut(l, omega, rOut, r0, rkOut)
                    if(e == 0): # fix for e=0 as the ODE solver stops at r0 (coming from rOut). The below line populates the `libration region' in this case
                        RkOut[lib_rkOut_pos[0]:lib_rkOut_pos[-1]+1] = RkOut[grid_rOut_pos[0]]
                else:
                    RkIn  = solInStatic(l, rkIn)
                    RkOut = solOutStatic(l, rkOut)
                        
                Wronskian = schw.f(schw.r_min) * (RkOut[lib_rkOut_pos[0],0]*RkIn[lib_rkIn2_pos[0],1] - RkIn[lib_rkIn2_pos[0],0]*RkOut[lib_rkOut_pos[0],1])    
                                
                # spherical harmonic - note it has an unusual argument order: m, l, phi, theta
                S = sph_harm(m, l, 0.0, pi/2.0)
                                
                # compute the weight coefficients which multiply the homogeneous solutions
                CIn, COut = 0, 0
                for k in range(0,schw.N):
                    CsumFactor = 1.0/(schw.r_k[k]*schw.ut_k[k]) * cos(omega * schw.t_k[k] - m * schw.phi_k[k]) * schw.dtdchi_k2[k]
                    if k <= schw.N/2:
                        k1 = lib_rkOut_pos[k]
                        k2 = lib_rkIn2_pos[k]
                    else:
                        k1 = lib_rkOut_pos[schw.N - k]
                        k2 = lib_rkIn2_pos[schw.N - k]     
                                      
                    CIn  += RkOut[k1,0]*CsumFactor
                    COut += RkIn[k2,0]*CsumFactor
                    #FIXME Surely there is no reason to integrate from chi = 0..2pi and instead can use 0..pi?
                
                CpreFactor     = -4.*pi*S/Wronskian * schw.Omega_r/(2.0*schw.N)                           
                CIn  *= CpreFactor
                COut *= CpreFactor
                
                 
                # construct the field near the horizon
                if omega !=0 :
                    for k in range(0,len(grid_rIn1)):
                        innerBCs        = bc.innerBCs(l, grid_rIn1[k], grid_x[k], omega)
                        r               = grid_rIn1[k]
                        field           = -2.*CIn*innerBCs[0]/r*cmath.exp(-1j*omega*grid_tIn1[k])
                        fieldIn1[k]     += n_fold_func(m,n,field)
                        field_rs_In1[k] += n_fold_func(m,n,-2.*CIn*(innerBCs[2]/r - (1.-2./r)*innerBCs[0]/r**2)*cmath.exp(-1j*omega*grid_tIn1[k]))
                        field_t_In1[k]  += n_fold_func(m,n,-1j*omega*field)                        
                else:
                    inner   = solInStatic(l,grid_rIn1)
                    for k in range(0,len(grid_rIn1)):
                        r               = grid_rIn1[k]
                        field           = -2.*CIn*inner[k,0]/r*cmath.exp(-1j*omega*grid_tIn1[k])
                        fieldIn1[k]     += n_fold_func(m,n,field)
                        field_rs_In1[k] += (1.-2./r)*n_fold_func(m,n,-2.*CIn*(inner[k,1]/r - inner[k,0]/r**2)*cmath.exp(-1j*omega*grid_tIn1[k]))
                        field_t_In1[k]  += n_fold_func(m,n,-1j*omega*field)
                                                
                # construct the field to the left of the particle
                for k in range(0,len(grid_rIn2)):
                    r               = grid_rIn2[k]
                    field           = -2.*CIn*RkIn[grid_rIn2_pos[k],0]/r*cmath.exp(-1j*omega*grid_tIn2[k])
                    fieldIn2[k]     += n_fold_func(m,n,field)
                    field_rs_In2[k] += (1.-2./r)*n_fold_func(m,n,-2.*CIn*(RkIn[grid_rIn2_pos[k],1]/r - RkIn[grid_rIn2_pos[k],0]/r**2 )*cmath.exp(-1j*omega*grid_tIn2[k]))
                    field_t_In2[k]  += n_fold_func(m,n,-1j*omega*field)                    

                # construct the field to the right of the particle
                for k in range(0,len(grid_rOut)):
                    r               = grid_rOut[k]
                    field           = -2.*COut*RkOut[grid_rOut_pos[k],0]/r*cmath.exp(-1j*omega*grid_tOut[k])
                    fieldOut[k]     += n_fold_func(m,n,field)
                    field_rs_Out[k] += (1.-2./r)*n_fold_func(m,n,-2.*COut*(RkOut[grid_rOut_pos[k],1]/r - RkOut[grid_rOut_pos[k],0]/r**2)*cmath.exp(-1j*omega*grid_tOut[k]))
                    field_t_Out[k]  += n_fold_func(m,n,-1j*omega*field)
                
                # construct the field near infinity
                if omega !=0:    
                    for k in range(0,len(grid_rOut2)):
                        r                = grid_rOut2[k]
                        outerBCs         = bc.outerBCs(l,grid_rOut2[k], omega)
                        field            = -2.*COut*outerBCs[0]/r*cmath.exp(-1j*omega*grid_tOut2[k])
                        fieldOut2[k]     += n_fold_func(m,n,field)
                        field_rs_Out2[k] += (1.-2./r)*n_fold_func(m,n,-2.*COut*(outerBCs[1]/r - outerBCs[0]/r**2)*cmath.exp(-1j*omega*grid_tOut2[k]))
                        field_t_Out2[k]  += n_fold_func(m,n,-1j*omega*field)
                else:
                    outer   = solOutStatic(l,grid_rOut2)
                    for k in range(0,len(grid_rOut2)):
                        r                = grid_rOut2[k]
                        field            = -2.*COut*outer[k,0]/r*cmath.exp(-1j*omega*grid_tOut2[k])
                        fieldOut2[k]     += n_fold_func(m,n,field)
                        field_rs_Out2[k] += (1.-2./r)*n_fold_func(m,n,-2.*COut*(outer[k,1]/r - outer[k,0]/r**2)*cmath.exp(-1j*omega*grid_tOut2[k]))
                        field_t_Out2[k]  += n_fold_func(m,n,-1j*omega*field)
                
                # check for convergence                
                if n == nmax:
                    converged = 1
                if abs(1. - fieldIn2[-1]/fieldOut[0]) < 1e-10:
                    converged = 1
                if n > 0:
                    if m == 0:
                        n += 1
                    else:
                        n = -n
                else:
                    n = -n + 1
                if e == 0:
                    converged = 1
                
            
            # plot the results
            # plt.plot(grid_rhoIn1, np.real(fieldIn1))
            # plt.plot(grid_rhoIn2, np.real(fieldIn2))
            # plt.plot(grid_rhoOut, np.real(fieldOut))
            # plt.plot(grid_rhoOut2, np.real(fieldOut2))
            # plt.grid()
            # plt.show()
            
            # output the results to disk    
            re_field = np.concatenate((np.real(fieldIn1), np.real(fieldIn2), np.real(fieldOut), np.real(fieldOut2)))
            im_field = np.concatenate((np.imag(fieldIn1), np.imag(fieldIn2), np.imag(fieldOut), np.imag(fieldOut2)))
            
            re_field_rs = np.concatenate((np.real(field_rs_In1), np.real(field_rs_In2), np.real(field_rs_Out), np.real(field_rs_Out2)))
            im_field_rs = np.concatenate((np.imag(field_rs_In1), np.imag(field_rs_In2), np.imag(field_rs_Out), np.imag(field_rs_Out2)))
            
            re_field_t = np.concatenate((np.real(field_t_In1), np.real(field_t_In2), np.real(field_t_Out), np.real(field_t_Out2)))
            im_field_t = np.concatenate((np.imag(field_t_In1), np.imag(field_t_In2), np.imag(field_t_Out), np.imag(field_t_Out2)))
           
            output_dir = 'data/output/' 
            output_filename = '%sSSF_init_data_p%g_e%g_n%g_l%dm%d.dat' % (output_dir, p, e, n_res, l, m)
            print("Outputting data to %s" % output_filename)
            print("")
            
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                            
            # create the output rho grid. This differs from grid_rho by only having two values of rmax rather than the 3 found in grid_rho
            grid_rho2 = np.concatenate(( grid_rhoIn1, grid_rhoIn2, grid_rhoOut, grid_rhoOut2))
            
            np.savetxt(output_filename, np.transpose( (grid_rho2,  re_field, im_field, re_field_rs, im_field_rs, re_field_t, im_field_t)  ))


def n_fold_func(m,n,field):
    if m == 0 and n != 0:
        return 2*np.real(field)
    else:
        return field            
            
# Analytically computes the static in solution and its radial derivative
def solInStatic(l, rkIn):
    RkIn = np.zeros(shape=(len(rkIn),2), dtype=complex)
    for k in range(0,len(rkIn)):
        r = rkIn[k]
        LegendreP = np.transpose(lpn(l, r - 1.))[-1]
        RkIn[k] = (r*LegendreP[0], r*LegendreP[1] + LegendreP[0]) 
    return RkIn   
    
# Analytically computes the static out solution and its radial deriviative
def solOutStatic(l, rkOut):
    RkOut = np.zeros(shape=(len(rkOut),2), dtype=complex)
    for k in range(0,len(rkOut)):
        r = rkOut[k]
        if l == 0:      #deal with a bug in scipy for l=0
            LegendreQ = np.transpose(lqn(2, r - 1.))[0]
        else:
             LegendreQ = np.transpose(lqn(l, r - 1.))[-1]
             
        RkOut[k] = (r*LegendreQ[0], r*LegendreQ[1] + LegendreQ[0])
    return RkOut                    
               
# Numerically solves for the solution and its radial derivative                    
def solIn(l, omega, rIn, r0, rkIn):
    ode_int = ode(field_eq).set_integrator('zvode', method='bdf', with_jacobian=False, rtol=1e-12)
    y0 = bc.innerBCs(l, rIn, 0, omega)[0:2]
    RkIn = np.zeros(shape=(len(rkIn),2), dtype=complex)
    RkIn[0] = y0    
    ode_int.set_initial_value(y0, rIn).set_f_params(l,omega)
    i = 1    
    while ode_int.successful() and i < len(rkIn):
        ode_int.integrate(rkIn[i])
        RkIn[i] = ode_int.y
        i += 1
    return RkIn                  
         
    
# Numerically solves for the out solution and its radial derivative    
def solOut(l, omega, rOut, r0, rkOut):
    ode_int = ode(field_eq).set_integrator('zvode', method='bdf', nsteps=1000000, rtol=1e-12)
    y0 = bc.outerBCs(l, rOut, omega)
    RkOut = np.zeros(shape=(len(rkOut),2), dtype=complex)
    RkOut[len(rkOut) - 1] = y0
    ode_int.set_initial_value(y0, rOut).set_f_params(l,omega)
    i = len(rkOut) - 2
    while ode_int.successful() and ode_int.t > rkOut[0]:
        ode_int.integrate(rkOut[i])
        RkOut[i] = ode_int.y
        i -= 1
    return RkOut

# The radial scalar field equation
def field_eq(r, y, l, omega):
    f = 1.-2./r
    return [y[1], -2*y[1]/((r-2.)*r) - ((2.-r)*(2.+r*l*(l+1.)) + r**4 * omega**2 )/(r**2*(r-2.)**2)*y[0] ]

    
if __name__ == "__main__": 
    main(sys.argv[1:])
