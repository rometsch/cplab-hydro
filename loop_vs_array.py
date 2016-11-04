# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 14:57:23 2016

@author: thomas
"""

import matplotlib.pyplot as plt
import advection

def initial_function(x):
    if (x<-1./3 or x>1./3):
        return 0.0;
    else:
        return 1.0;
        
#   Make sure upwind 2nd order with loop calculation 
#   yields same result as with array calculation
def compare_loop_vs_array():
    plt.figure();
    plt.clf();
    labels = ["analytic"];
    style = ['.','x'];
    methods = ['upwind_2nd_loop','upwind_2nd_arr'];
    T=40;
    # Plot initial function:
    initial = advection.Box(initial_function,[-1,1],400,1.0);
    plt.plot(initial.get_x(),initial.get_Psi(),'-')
    for n in range(2):
        box = advection.Box(initial_function,[-1,1],400,1.0);
        box.integrate(T,methods[n]);
        labels.append(methods[n]);
        plt.plot(box.get_x(),box.get_Psi(),style[n])
    plt.legend(labels,loc=(0.75,0.5));
    plt.xlabel("x");
    plt.ylabel("u");
    plt.title("T = {}".format(T))
    plt.draw();
    plt.savefig("loop vs array.pdf".format(T));
    
compare_loop_vs_array();

#%% Study time cost of loop vs array for 2nd order Upwind:

import timeit
N = 2;
# Measure time to integrate untill T=4, average over N runs:
T_upwind = timeit.timeit('box.integrate(4,\'upwind\')',setup='import advection; from __main__ import initial_function; box=advection.Box(initial_function,[-1,1],4000,1.0);',number=N)/N
T_lax_wendroff = timeit.timeit('box.integrate(4,\'lax_wendroff\')',setup='import advection; from __main__ import initial_function; box=advection.Box(initial_function,[-1,1],4000,1.0);',number=N)/N
T_upwind_2nd_loop = timeit.timeit('box.integrate(4,\'upwind_2nd_loop\')',setup='import advection; from __main__ import initial_function; box=advection.Box(initial_function,[-1,1],4000,1.0);',number=N)/N
T_upwind_2nd_arr = timeit.timeit('box.integrate(4,\'upwind_2nd_arr\')',setup='import advection; from __main__ import initial_function; box=advection.Box(initial_function,[-1,1],4000,1.0);',number=N)/N

print("Excecution times of integration from T=0 to T=4, (N=400,sigma=0.8)");
print("Upwind:\t\t\t {:.9f}".format(T_upwind));
print("Lax Werndroff:\t\t {:.9f}".format(T_lax_wendroff));
print("Upwind 2nd ord loop:\t {:.9f}".format(T_upwind_2nd_loop));
print("Upwind 2nd ord array:\t {:.9f}".format(T_upwind_2nd_arr));
print("\nRatio loop to array =\t {}".format(T_upwind_2nd_loop/T_upwind_2nd_arr))
