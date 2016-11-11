# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 19:55:44 2016

@author: thomas
"""

import advection;
import shocktube;
import matplotlib.pyplot as plt;

def rho0(x):
    if (x<-1./3 or x>1./3):
        return 0.0;
    else:
        return 1.0;
        
def u0(x):
    return 1.0;
    
def zero_function(x):
    return 0.0;

T = 0.016;

test_advection_tube = advection.Box(rho0,[-1,1],100,1.0);
test_advection_tube.integrate(T,'upwind_2nd_arr');

dt = 0.8/1*2.0/100

test_shocktube_tube = shocktube.Box(rho0,u0,zero_function,[-1,1],100,dt);
while test_shocktube_tube.t < T:
    test_shocktube_tube.t += dt;
    test_shocktube_tube.advection_rho();
    
print(sum(test_advection_tube.get_Psi()-test_shocktube_tube.get_rho()));

plt.plot(test_shocktube_tube.get_x(),test_shocktube_tube.u[2:-2]);
plt.draw();