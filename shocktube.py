# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:22:23 2016

@author: thomas

Simulation box for a 1D shocktube problem with methods to solve
the Euler equations using operator splitting and the upwind method.
"""

import numpy as np

class Box():
    
    def __init__(self,rho0,u0,eps0,gamma,I,N,dt):
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
        self.gamma = gamma;
        self.N = N;
        self.I = I;
        self.dt = dt;
        self.dx = (self.I[1] - self.I[0])/self.N;
        # Variables:
        self.t = 0; # integration time
        self.X = np.zeros(self.N+4);
        # Use two ghost cells on each side.
        # Allocate array for flux to be calculated by chosen method.
        self.rho = np.zeros(self.N+4);
        self.u = np.zeros(self.N+4);
        self.eps = np.zeros(self.N+4);
        self.delta = np.zeros(self.N+4);
        self.flux = np.zeros(self.N+4);
        # Initialize grid and function:
        self.X[2:-2] = self.I[0] + self.dx*np.arange(self.N);
        self.update_boundary(self.X);
#        self.rho[2:-2] = self.rho0(self.X[2:-2]);
#        self.u[2:-2] = self.u0(self.X[2:-2]);
#        self.eps[2:-2] = self.eps0(self.X[2:-2]);
        for f,f0 in ((self.rho,self.rho0),(self.u,self.u0),(self.eps,self.eps0)):
            f[2:-2] = f0(self.X[2:-2]);
        self.update_boundary_all();   
        
  
  
    def integrate(self,T):
        ''' @param:
           T: time to integrate up to
        '''
        
        while self.t<T: # integrate over time T.
            # Calculate flux:            
            self.t += self.dt;
            self.integration_step();
            
        
    def integration_step(self):
        '''Manage advection of density, velocity and energy.'''
        
        # Advection with 2nd order Upwind scheme
        rho_new = np.zeros(self.N+4);
        flux_mass = self.advection_scalar(self.rho)*self.u;  
        rho_new[2:-2] = self.rho[2:-2]-self.dt/self.dx*(flux_mass[3:-1]-flux_mass[2:-2]);
        self.update_boundary(rho_new);

        u_new = np.zeros(self.N+4);
        flux_momentum = np.zeros(self.N+4);
        flux_momentum[1:-2] = 0.5*self.advection_velocity(self.u)[1:-2]*(flux_mass[1:-2]+flux_mass[2:-1]);
        rho_avg = np.zeros(self.N+4);
        rho_avg[2:-2] = 0.5*(rho_new[1:-3]+rho_new[2:-2]);
        if (rho_avg[2:-2]==0).any():
            print("Division by zero encountered in advection of velcoity");
        u_new[2:-2] = (self.u[2:-2]*rho_avg[2:-2]-self.dt/self.dx*(flux_momentum[2:-2]-flux_momentum[1:-3]))/rho_avg[2:-2];
        self.update_boundary(u_new);
        
        eps_new = np.zeros(self.N+4);
        flux_eps = flux_momentum*self.advection_scalar(self.eps);
        if (rho_avg[2:-2]==0).any():
            print("Division by zero encountered in advection of energy");
        eps_new[2:-2] = (self.eps[2:-2]*self.rho[2:-2]-self.dt/self.dx*(flux_eps[3:-1]-flux_eps[2:-2]))/rho_new[2:-2];
        self.update_boundary(eps_new);
        
        self.rho = rho_new;        
        
        # force and pressure terms
        p = (self.gamma-1)*rho_new*eps_new;
        self.u[2:-2] = u_new[2:-2] - self.dt/self.dx/rho_avg[2:-2]*(p[2:-2]-p[1:-3]);
        self.update_boundary(self.u);

        self.eps[2:-2] = eps_new[2:-2] - self.dt/self.dx*p[2:-2]/rho_new[2:-2]*(u_new[3:-1]-u_new[2:-2]);        
        self.update_boundary(self.eps);         

    def advection_velocity(self,Y):
        '''Calculate the flux of a vector variable using the
        2nd order upwind scheme'''
        delta = np.zeros(self.N+4);
        avg = np.zeros(self.N+4);
        adv = np.zeros(self.N+4);
        # calculate geometrical mean:
        for n in range(1,self.N+3):
            nom = (Y[n+1]-Y[n])*(Y[n]-Y[n-1]);
            denom = (Y[n+1]-Y[n-1]);
            if nom>0:
                if (denom == 0):
                    print("Encountered division by 0 in advection step of ",Y.__name__);
                    delta[n] = 0;
                else:
                    delta[n] = 2*nom/denom;
            else:
                delta[n] = 0;
        # calculate average at cell boundary:
        avg[0:-1] = 0.5*(Y[1:]+Y[0:-1]);
        # calculate advected quantities:
        for n in range(1,self.N+2):
            if avg[n] > 0 :
                adv[n] = Y[n]+0.5*(1-avg[n]*self.dt/self.dx)*delta[n];
            else:
                adv[n] = Y[n+1] - 0.5*(1-avg[n]*self.dt/self.dx)*delta[n+1];
        # return resulting flux:
        return adv;
        
    def advection_scalar(self,Y):
        '''Calculate the flux of a scalar variable using the
        2nd order upwind scheme'''
        delta = np.zeros(self.N+4);
        adv = np.zeros(self.N+4);
        # calculate geometrical mean:
        for n in range(1,self.N+3):
            nom = (Y[n+1]-Y[n])*(Y[n]-Y[n-1]);
            denom = (Y[n+1]-Y[n-1]);
            if nom>0:
                if (denom == 0):
                    print("Encountered division by 0 in advection step of ",Y.__name__);
                    delta[n] = 0;
                else:
                    delta[n] = 2*nom/denom;
            else:
                delta[n] = 0;
        # calculate advected quantities:
        for n in range(2,self.N+3):
            if self.u[n] > 0 :
                adv[n] = Y[n-1]+0.5*(1-self.u[n]*self.dt/self.dx)*delta[n-1];
            else:
                adv[n] = Y[n] - 0.5*(1-self.u[n]*self.dt/self.dx)*delta[n];
        # return resulting flux:
        return adv;
                            
        
    def update_boundary_all(self):
        '''Update boundaries of all variables.'''
        self.update_boundary_rho();
        self.update_boundary_u();
        self.update_boundary_eps();

    def update_boundary_rho(self):
        '''Update ghost cells of rho.'''
        self.update_boundary(self.rho);

    def update_boundary_u(self):
        '''Update ghost cells of u.'''
        self.update_boundary(self.u);
        
    def update_boundary_eps(self):
        '''Update ghost cells of eps'''
        self.update_boundary(self.eps);

    def update_boundary(self,arr):
        ''' Update ghost cells of variable stored in arr. 
        Use periodic boundary conditions. '''
        # Update ghost cells:
        arr[0] = arr[-4];
        arr[1] = arr[-3];
        arr[-2] = arr[2];
        arr[-1] = arr[3];
            
    def get_rho(self):
        '''Return array with mass density values.'''
        return self.rho[2:-2];

    def get_u(self):
        '''Return array with velocity values.'''
        return self.u[2:-2];
        
    def get_eps(self):
        '''Return array with energy density values.'''
        return self.eps[2:-2];
        
    def get_p(self):
        '''Return array with pressure values.'''
        return (self.gamma-1)*self.rho[2:-2]*self.rho[2:-2];
        
    def get_T(self):
        '''Return array with pressure values.'''
        return (self.gamma-1)*self.eps[2:-2];
        
    def get_x(self):
        return self.X[2:-2];
            