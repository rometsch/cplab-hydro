# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:22:23 2016

@author: thomas

Simulation box for a 1D shocktube problem with methods to solve
the Euler equations using operator splitting and the upwind method.
"""

import numpy as np

class Box():
    
    def __init__(self,rho0,u0,eps0,I,N,dt):
        ''' 
            Initialize object with given parameters.
            @param:
           rho0,u0,eps0: initial conditions, real function 
           I: computational domain, Interval [xmin,xmax]
           N: Number of grid points
           dt: constant time step
        '''
        self.rho0 = rho0;
        self.u0 = u0;
        self.eps0 = eps0;
        # Parameter:
        self.N = N;
        self.I = I;
        self.dt = dt;
        self.dx = (self.I[1] - self.I[0])/self.N;
        # Variables:
        self.t = 0; # integration time
        self.X = np.zeros(self.N+4);
        # Use two ghost cells on each side.
        self.rho = np.zeros(self.N+4);
        self.u = np.zeros(self.N+4);
        self.eps = np.zeros(self.N+4);
        self.flux = np.zeros(self.N+4); # Allocate array for flux to be calculated by chosen method.
        # Initialize grid and function:
        for n in range(self.N):
            self.X[n+2] = self.I[0] + n*self.dx;
        for n in range(self.N-1):            
            self.rho[n+2] = self.rho0(self.X[n+2]+0.5/self.N);
            self.u[n+2] = self.u0(self.X[n+2]+0.5/self.N);
            self.eps[n+2] = self.eps0(self.X[n+2]+0.5/self.N);
        self.update_boundaries();   
        
        self.flux_functions = {"upwind" : self.flux_upwind,'lax_wendroff' : self.flux_lax_wendroff,
                               'upwind_2nd_loop' : self.flux_upwind_2nd_loop, 'upwind_2nd_arr' : self.flux_upwind_2nd_array};
    
    
    
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

    def update_boundary_all(self):
        '''Update boundaries of all variables.'''
        self.update_boundary_rho();
        self.update_boundary_u();
        self.update_boundary_eps();

    def update_boundary_rho(self):
        '''Update ghost cells of rho.
            Use periodic boundary conditions.
        '''
        self.update_boundary(self.rho);

    def update_boundary_u(self):
        '''Update ghost cells of u.
            Use periodic boundary conditions.
        '''
        self.update_boundary(self.u);
        
    def update_boundary_eps(self):
        '''Update ghost cells of eps.
            Use periodic boundary conditions.
        '''
        self.update_boundary(self.eps);

    def update_boundary(arr):
        ''' Update ghost cells of variable stored in arr. 
        Use periodic boundary conditions. '''
        # Update ghost cells:
        arr[0] = arr[-4];
        arr[1] = arr[-3];
        arr[-2] = arr[2];
        arr[-1] = arr[3];
            
    def get_rho(self):
        '''Return array with rho values.'''
        return self.rho[2:-2];

    def get_u(self):
        '''Return array with u values.'''
        return self.u[2:-2];
        
    def get_eps(self):
        '''Return array with eps values.'''
        return self.eps[2:-2];
        
    def get_x(self):
        return self.X[2:-2];
            