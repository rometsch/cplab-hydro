# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 15:01:47 2016

@author: thomas
"""

import shocktube;
#import numpy as np;
import matplotlib.pyplot as plt;

def rho0(x):
    return 1.0*(x<=0.5) + 0.125*(x>0.5);
        
def u0(x):
    return 0.0;
    
def eps0(x):
    return 2.5*(x<=0.5) + 2.0*(x>0.5);

dt = 0.001;
I=[0.0,1.0];
t_end = 0.228;
N_cells = 100;
gamma = 1.4;

box = shocktube.Box(rho0,u0,eps0,gamma,I,N_cells,dt);
box.integrate(t_end);

x = box.get_x();
rho = box.get_rho();
u = box.get_u();
T = box.get_T();
p = box.get_p();
eps = box.get_eps();

data_string = "# Numerical solution of shocktube problem after time t = {0:.3e}\n".format(t_end);
data_string += "# I    x(i)      u(i)       rho(i)   Temp(i)   Pgas(i)\n"
for n in range(N_cells):
    data_string += "{0:03d}\t{1:.3e}\t{2:.3e}\t{3:.3e}\t{4:.3e}\t{5:.3e}\n".format(n,x[n],u[n],rho[n],T[n],p[n]);

with open("shocktube_sol_num.txt", "w") as text_file:
    text_file.write(data_string);


plt.plot(x,rho);
plt.plot(x,u);
plt.plot(x,eps);
plt.legend(["rho","u","eps"]);
plt.draw();
