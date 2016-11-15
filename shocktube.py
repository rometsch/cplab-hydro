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
        self.Xa = np.zeros(self.N+4);
        self.Xb = np.zeros(self.N+4);
        # Use two ghost cells on each side.
        # Allocate array for flux to be calculated by chosen method.
        self.rho = np.zeros(self.N+4);
        self.u = np.zeros(self.N+4);
        self.eps = np.zeros(self.N+4);
        self.delta = np.zeros(self.N+4);
        self.flux = np.zeros(self.N+4);
        # Initialize grid and function:
        self.Xa = self.I[0] + self.dx*(-2+np.arange(self.N+4));
        self.Xb = self.Xa +0.5*self.dx;
        self.rho[2:-2] = self.rho0(self.Xb[2:-2]);
        self.update_boundary_rho(self.rho);
        self.u[2:-2] = self.u0(self.Xa[2:-2]);
        self.update_boundary_velocity(self.u);
        self.eps[2:-2] = self.eps0(self.Xb[2:-2]);
        self.update_boundary_eps(self.eps);

        
  
  
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
        N = self.N;
        # Advection with 2nd order Upwind scheme
        rho_new = np.zeros(N+4);
        flux_mass = np.zeros(N+4);
        flux_mass[2:N+3] = self.advection_scalar(self.rho)[2:N+3]*self.u[2:N+3];  
        rho_new[2:N+2] = self.rho[2:N+2]-self.dt/self.dx*(flux_mass[3:N+3]-flux_mass[2:N+2]);
#        self.update_boundary_rho(rho_new);

        u_new = np.zeros(N+4);
        flux_momentum = np.zeros(N+4);
        flux_momentum[2:N+2] = 0.5*self.advection_velocity(self.u)[2:N+2]*(flux_mass[2:N+2]+flux_mass[3:N+3]);
        rho_avg_new = self.average_density(rho_new);
        rho_avg_old = self.average_density(self.rho);        
        if (rho_avg_new[3:N+2]==0).any():
            print("t = {:.3e} : Division by zero encountered in advection of velcoity".format(self.t));
        u_new[3:N+2] = (self.u[3:N+2]*rho_avg_old[3:N+2]-self.dt/self.dx*(flux_momentum[3:N+2]-flux_momentum[2:N+1]))/rho_avg_new[3:N+2];      
#        u_new[3:-2] = (self.u[3:-2]*0.5*(self.rho[2:-3]+self.rho[3:-2])-self.dt/self.dx*(flux_momentum[3:-2]-flux_momentum[2:-3]))/0.5/(rho_new[2:-3]+rho_new[3:-2]);      
#        self.update_boundary_velocity(u_new);        
        
        eps_new = np.zeros(self.N+4);
        flux_eps = np.zeros(self.N+4);
        flux_eps[2:-1] = flux_mass[2:-1]*self.advection_scalar(self.eps)[2:-1];
        if (rho_new[2:-2]==0).any():
            print("t = {:.3e} : Division by zero encountered in advection of energy".format(self.t));
        eps_new[2:-2] = (self.eps[2:-2]*self.rho[2:-2]-self.dt/self.dx*(flux_eps[3:-1]-flux_eps[2:-2]))/rho_new[2:-2];
#        self.update_boundary_eps(eps_new);
        
        # force and pressure terms
        p = self.calc_pressure(rho_new,eps_new);
        if (rho_avg_new[3:-2]==0).any():
            print("t = {:.3e} : Division by zero encountered force and pressure calculation for velocity".format(self.t));
        self.u[3:-2] = u_new[3:-2] - self.dt/self.dx*(p[3:-2]-p[2:-3])/rho_avg_new[3:-2];
#        self.update_boundary_velocity(self.u);

        self.eps[2:-2] = eps_new[2:-2]*(1 - self.dt/self.dx*(self.gamma-1)*(u_new[3:-1]-u_new[2:-2]));        
#        self.eps[2:-2] = eps_new[2:-2] - self.dt/self.dx*p[2:-2]/rho_new[2:-2]*(u_new[3:-1]-u_new[2:-2]);        
#        self.update_boundary_eps(self.eps);

        self.rho = rho_new;
        self.update_boundary_rho(self.rho);
        self.update_boundary_velocity(self.u);
        self.update_boundary_eps(self.eps);
        
    def calc_pressure(self,rho,eps):
        '''Calculate the pressure at the cell center'''
        p = np.zeros(self.N+4);
        p[:] = (self.gamma-1)*rho[:]*eps[:];
        return p;
        
    def average_density(self,rho):
        '''Calculate average density for cell boundaries'''
        rho_avg = np.zeros(self.N+4);
        rho_avg[1:] = 0.5*(rho[1:]+rho[:-1]);
        return rho_avg;
        
    def average_velocity(self,u):
        '''Calcualte average velocity at cell center'''
        u_avg = np.zeros(self.N+4);
        u_avg[:-1] = 0.5*(u[:-1]+u[1:]);
        return u_avg;

    def advection_velocity(self,Y):
        '''Calculate the flux of a vector variable using the
        2nd order upwind scheme'''
        adv = np.zeros(self.N+4);
        delta = self.undivided_difference(Y);
        # calculate average at cell boundary:
        u_avg = self.average_velocity(self.u);
        # calculate advected quantities:
        for n in range(2,self.N+2):
            if u_avg[n] > 0 :
                adv[n] = Y[n]+0.5*(1-u_avg[n]*self.dt/self.dx)*delta[n];
            else:
                adv[n] = Y[n+1]-0.5*(1+u_avg[n]*self.dt/self.dx)*delta[n+1];
        # return resulting flux:        
        return adv;
        
    def advection_scalar(self,Y):
        '''Calculate the flux of a scalar variable using the
        2nd order upwind scheme'''
        adv = np.zeros(self.N+4);
        delta = self.undivided_difference(Y);
        # calculate advected quantities:
        for n in range(2,self.N+3):
            if self.u[n] > 0 :
                adv[n] = Y[n-1]+0.5*(1-self.u[n]*self.dt/self.dx)*delta[n-1];
            else:
                adv[n] = Y[n]-0.5*(1+self.u[n]*self.dt/self.dx)*delta[n];
        return adv;
                    
    def undivided_difference(self,x):
        '''Calculate the undivided difference of quantity x for use in advection step following van Leer 1972.
        x must be an array of length self.N+4.'''
        delta = np.zeros(self.N+4);
        for n in range(1,self.N+3):
            nom = (x[n+1]-x[n])*(x[n]-x[n-1]);
            denom = (x[n+1]-x[n-1]);
            if nom>0:
                if (denom == 0):
                    print("Encountered division by 0 in calculation of undivided difference of ",x.__name__);
                    delta[n] = 0;
                else:
                    delta[n] = 2*nom/denom;
            else:
                delta[n] = 0;
        return delta;
        
    def update_boundary_all(self):
        '''Update boundaries of all variables.'''
        self.update_boundary_rho(self.rho);
        self.update_boundary_velocity(self.u);
        self.update_boundary_eps(self.eps);

    def update_boundary_rho(self,rho):
        '''Update ghost cells of rho.'''
        rho[1] = rho[2];
        rho[self.N+2] = rho[self.N+1];

    def update_boundary_velocity(self,u):
        '''Update ghost cells of u.'''
        u[1] = -u[3];
        u[2] = 0;
        u[self.N+3] = -u[self.N+1];
        u[self.N+2] = 0;
        
    def update_boundary_eps(self,eps):
        '''Update ghost cells of eps'''
        eps[1] = eps[2];
        eps[self.N+2] = eps[self.N+1];
            
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
        return (self.gamma-1)*self.rho[2:-2]*self.eps[2:-2];
        
    def get_T(self):
        '''Return array with pressure values.'''
        return (self.gamma-1)*self.eps[2:-2];
        
    def get_x(self):
        return self.Xb[2:-2];
            