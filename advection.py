# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 07:51:23 2016

@author: thomas

make sure to use Delta Psi_(n-1) and not for n!
"""

import numpy as np

class Box():
    
    def __init__(self,func,I,N,a):
        ''' @param:
           f: initial condition, real function 
           I: computational domain, Interval [xmin,xmax]
           N: Number of grid points
           a: constant velocity
        '''
        self.func = func;
        # Parameter:
        self.N = N;
        self.I = I;
        self.a = a;
        self.sigma = 0.8;   # Pick sigma.
        self.dx = (self.I[1] - self.I[0])/self.N;
        # Variables:
        self.t = 0; # integration time
        self.X = np.zeros(self.N+4);
        # Use two ghost cells on each side.
        self.Psi = np.zeros(self.N+4);
        self.flux = np.zeros(self.N+4); # Allocate array for flux to be calculated by chosen method.
        # Initialize grid and function:
        for n in range(self.N):
            self.X[n+2] = self.I[0] + n*self.dx;
        for n in range(self.N-1):            
            self.Psi[n+2] = self.func(self.X[n+2]+0.5/self.N);
        self.update_boundary_Psi();   
        
        self.flux_functions = {"upwind" : self.flux_upwind,'lax_wendroff' : self.flux_lax_wendroff,
                               'upwind_2nd_loop' : self.flux_upwind_2nd_loop, 'upwind_2nd_arr' : self.flux_upwind_2nd_array};
    
    
    def update_boundary_Psi(self):
        ''' Use periodic boundary conditions. '''
        # Update ghost cells:
        self.Psi[0] = self.Psi[-4];
        self.Psi[1] = self.Psi[-3];
        self.Psi[-2] = self.Psi[2];
        self.Psi[-1] = self.Psi[3];

    def update_boundary_flux(self):
        ''' Use periodic boundary conditions. '''
        # Update ghost cells:
        self.flux[0] = self.flux[-4];
        self.flux[1] = self.flux[-3];
        self.flux[-2] = self.flux[2];
        self.flux[-1] = self.flux[3];
    
    def integrate(self,T,method):
        ''' @param:
           T: time to integrate
           method: method of choice
        '''
        # Pick method.
        calc_flux = self.flux_functions[method];

        # Calculate dt:
        dt = self.sigma/self.a*self.dx;
        
        while self.t<T: # integrate over time T.
            # Calculate flux:            
            self.t += dt;
            calc_flux();
            self.Psi[2:-2] += self.sigma*(self.flux[2:-2]-self.flux[3:-1]);
            self.update_boundary_Psi();


    def flux_upwind(self):
        ''' Upwind method
            Calculate flux divided by velocity a.
            For Upwind Delta_Psi is zero.
        '''
        self.flux[2:-2] = self.Psi[1:-3];     
        self.update_boundary_flux();
            
    def flux_upwind_2nd_loop(self):
        ''' 2nd order Upwind method
            Calculate flux divided by velocity a.
        '''
        kappa = 1-self.sigma;
        for n in range(2,self.N+2):
            k = (self.Psi[n]-self.Psi[n-1])*(self.Psi[n-1]-self.Psi[n-2]);
            if k > 0:
                self.flux[n] = self.Psi[n-1]+kappa*k/(self.Psi[n]-self.Psi[n-2]);
            else:
                self.flux[n] = self.Psi[n-1];
        self.update_boundary_flux();
        
    def flux_upwind_2nd_array(self):
        ''' 2nd order Upwind method
            Calculate flux divided by velocity a.
            Use numpy arrays for performance.
        '''
        self.flux[2:-2] = (1-self.sigma)*(self.Psi[2:-2]-self.Psi[1:-3])*(self.Psi[1:-3]-self.Psi[0:-4]);
        self.flux[2:-2] *= (self.flux[2:-2]>0)*(self.Psi[2:-2]-self.Psi[0:-4]!=0)/(self.Psi[2:-2]-self.Psi[0:-4] + (self.Psi[2:-2]-self.Psi[0:-4]==0));
        self.flux[2:-2] += self.Psi[1:-3];
        self.update_boundary_flux();
    
    def flux_lax_wendroff(self):
        ''' Lax-Wendroff method
            Calculate flux divided by velocity a.
        '''
        self.flux[2:-2] = self.Psi[1:-3]+0.5*(1 - self.sigma)*(self.Psi[2:-2]-self.Psi[1:-3]);
        self.update_boundary_flux();

            
    def get_Psi(self):
        return self.Psi[2:-2];
        
    def get_x(self):
        return self.X[2:-2];
            