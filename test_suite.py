# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 19:55:44 2016

@author: thomas
"""

import advection;
import shocktube;
import numpy as np;
import matplotlib.pyplot as plt;

def rho0(x):
    return np.abs(x)<1.0/3;
        
def u0(x):
    return 1.0;
    
def zero_function(x):
    return 0.0;

dt = 0.8/1*2.0/100

T = 0.15;

test_advection_tube = advection.Box(rho0,[-1,1],100,1.0);
test_advection_tube.integrate(T,'upwind_2nd_arr');


test_shocktube_tube = shocktube.Box(rho0,u0,zero_function,[-1,1],100,dt);
test_shocktube_tube.integrate(T);

diff=sum(np.abs(test_advection_tube.get_Psi()-test_shocktube_tube.get_rho()));

print("Sum of difference of the modulus of density arrays",diff);

plt.plot(test_shocktube_tube.get_x(),test_shocktube_tube.get_rho());
plt.plot(test_advection_tube.get_x(),test_advection_tube.get_Psi());
plt.legend(['new','old']);
plt.draw();

#plt.plot(test_shocktube_tube.get_x(),test_shocktube_tube.delta_rho[2:-2]);
#plt.draw()