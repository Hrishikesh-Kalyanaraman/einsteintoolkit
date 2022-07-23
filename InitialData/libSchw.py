# -*- coding: utf-8 -*-
from math import *
from scipy.fftpack import dct
import numpy as np

class Schw:
    
    def __init__(self,p,e,kappa):
        self.p = p
        self.e = e
        self.kappa = kappa          # kappa is used to for accelerated orbits. Set to 1 for geodesic motion
        self.r_min = p/(1+e)
        self.r_max = p/(1-e)
        
        self.N = 400     #FIXME this should depend upon the eccentricity -- seems ok to make it high with machine precision numbers
        self.chi_k      = pi*np.linspace(0, 1, num=self.N)
        self.chi_k2     = 2.*pi*np.linspace(0, 1, num=self.N+1)
        
        self.r_v        = np.vectorize(self.r)        
        self.r_k        = self.r_v(self.chi_k2)
                        
        dtdchi_v        = np.vectorize(self.dtdchi)
        self.dtdchi_k   = dtdchi_v(self.chi_k)
                
        dphidchi_v      = np.vectorize(self.dphidchi)
        self.dphidchi_k = dphidchi_v(self.chi_k)
        
        self.dtdchi_k2  = dtdchi_v(self.chi_k2)
        self.dphidchi_k2= dphidchi_v(self.chi_k2)
        
        self.Gtn        = self.DCT_I(self.dtdchi_k)
        self.Gphin      = self.DCT_I(self.dphidchi_k)
        
        self.t_v        = np.vectorize(self.t)
        self.t_k        = self.t_v(self.chi_k2)
        
        self.phi_v      = np.vectorize(self.phi)
        self.phi_k      = self.phi_v(self.chi_k2)
        
        self.ut_v       = np.vectorize(self.ut)
        self.ut_k       = self.ut_v(self.chi_k2)
        
        self.T          = self.t(2*pi)
        self.Deltaphi   = self.phi(2*pi)
        
        self.Omega_r    = 2*pi/self.T
        self.Omega_phi  = self.Deltaphi/self.T        
        
    # Output the orbital parameters    
    def print_orbit_parameters(self):
        print("p:    %.15lf" % self.p)
        print("e:    %.15lf" % self.e)
        print("rmin: %.15lf" % self.r_min)
        print("rmax: %.15lf" % self.r_max)
        print("Ωφ:   %.15lf" % self.Omega_phi)
        print("Ωr:   %.15lf" % self.Omega_r)
        print("")
    
    # The radial coordinate along the orbit    
    def r(self, chi):
        return self.p/(1 + self.e*cos(chi))
    
    # The Schwarzschild f(r) function    
    def f(self, r):
        return 1. - 2./r
    
    # The Schwarzschild f(r(chi)) function    
    def fchi(self, chi):
        return 1 - 2/self.r(chi)
    
    # Discrete cosine transform with unitary normalization     
    def DCT_I(self, array):
        n = len(array)
        return dct(array,1)*sqrt(2./(n-1))/2.
    
    # The geodesic dt/dchi    
    def dtdchi_geo(self, chi):
        p=self.p
        e=self.e
        return (p**2*sqrt(-2 - 2*e + p)*sqrt(-2 + 2*e + p))/(sqrt(-6 + p - 2*e*cos(chi))*(-2 + p - 2*e*cos(chi))*(1 + e*cos(chi))**2)
    
    # dt/dchi
    def dtdchi(self, chi):
        return self.kappa*self.dtdchi_geo(chi)
    
    # dt/dtau   
    def dtaudchi(self, chi):
        p=self.p
        e=self.e
        return ((1. - 2./self.r(chi))*self.dtdchi(chi)**2 - ((p*e*sin(chi))/(1.+e*cos(chi))**2)**2/(1.-2./self.r(chi)) - self.r(chi)**2*p*(p-6-2*e*cos(chi))**-1)**0.5
    
    # the contravariant t-component of the four-velocity    
    def ut(self, chi):
        return self.dtdchi(chi)/self.dtaudchi(chi)
        
    # dphi/dchi
    def dphidchi(self, chi):
        return sqrt(self.p)/sqrt(-6 + self.p - 2*self.e*cos(chi))
    
    # coordinate time t along the orbit   
    def t(self, chi):
        N = self.N
        sum = 0;
        for n in range(1, N-1):
            sum += 1./n*self.Gtn[n]*sin(n*chi)
        return sqrt(2./(N -1))*(0.5*self.Gtn[0]*chi + 0.5*self.Gtn[N-1]*sin((N -1.)*chi)/(N-1.) + sum)
    
    # the phi coordinate along the orbit    
    def phi(self, chi):
        N = self.N
        sum = 0;
        for n in range(1, N-1):
            sum += 1./n*self.Gphin[n]*sin(n*chi)
        return sqrt(2./(N -1))*(0.5*self.Gphin[0]*chi + 0.5*self.Gphin[N-1]*sin((N -1.)*chi)/(N-1.) + sum)
