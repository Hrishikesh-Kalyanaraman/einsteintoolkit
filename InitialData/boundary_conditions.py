import numpy as np
from math import *
import cmath

# computes the tortoise coordinate r_* as a function of the Schwarzschild radial coordinate
def rstar(r):
    return r + 2.*log(r/2. - 1.)

# computes the inner boundary conditions
def innerBCs(l, rIn, rsIn, omega):
    b = [0,0,1,0]
    k = 1
    phi_in = b[2]
    dphi_in = 0
    dphi_drs_in = 0
    fIn = 1. - 2./rIn
    if rIn - 2. != 0:
        while k < 12:               #FIXME add some convergence criteria
            b[3] = -1./(2*k*(k - 4j*omega))*(((2*k - 1)*(k-2) - l*(l+1) - 12j*omega*(k-1))*b[2] + ((k-2)*(k-3)/2. - l*(l+1)/2. - 6j*omega*(k-2))*b[1] - 1j*omega*(k-3)*b[0])
        
            phi_in  += b[3]*(rIn-2.)**k
            dphi_in += b[3]*k*(rIn-2.)**(k-1) - 1j*omega/fIn*b[2]*(rIn - 2.)**(k-1)
            dphi_drs_in += 1./rIn*(b[3]*k*(rIn-2.)**(k) - 1j*omega*b[2]*rIn*(rIn - 2.)**(k-1))
    
            b = np.roll(b,-1)
            k += 1
    else:
        phi_in = 1
        dphi_drs_in = -1j*omega
        
    if rsIn == 0:
        rsIn = rstar(rIn)
                
    phi_in  *= cmath.exp(-1j*omega*rsIn)
    dphi_in *= cmath.exp(-1j*omega*rsIn)
    
    dphi_drs_in *= cmath.exp(-1j*omega*rsIn)
        
    return [phi_in, dphi_in, dphi_drs_in]
    
# computes the outer boundary conditions    
def outerBCs(l, rOut, omega):
    a = [0, 1, 0]
    k=1
    phi_out = a[1]
    dphi_out = 0
    fOut = 1. - 2./rOut
    while k < 12:           #FIXME add some convergence criteria
        a[2] = 1j/(2.*omega*k)*((-k*(k-1.) + l*(l+1.))*a[1] + 2.*(k-1.)**2*a[0])
        
        phi_out  += a[2]*rOut**(-k)    
        dphi_out += -k*rOut**(-k-1)*a[2] + (1j*omega/fOut)*a[1]*rOut**(-k+1) 
        a = np.roll(a,-1)
        k += 1
                
    rsOut = rstar(rOut)    
        
    phi_out  *= cmath.exp(1j*omega*rsOut)
    dphi_out *= cmath.exp(1j*omega*rsOut)    
    
    return [phi_out, dphi_out]
    