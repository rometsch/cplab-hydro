# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 17:33:06 2016

@author: thomas
"""

import matplotlib.pyplot as plt
import advection

def initial_function(x):
    if (x<-1./3 or x>1./3):
        return 0.0;
    else:
        return 1.0;


def plot_after_time(T,methods):
    # Plot integrated function after time T for all methods in method:
    plt.figure();
    plt.clf();
    labels = ["analytic"];
    # Plot initial function:
    initial = advection.Box(initial_function,[-1,1],400,1.0);
    plt.plot(initial.get_x(),initial.get_Psi(),'-')
    for method in methods:
        box = advection.Box(initial_function,[-1,1],400,1.0);
        box.integrate(T,method);
        labels.append(method);
        plt.plot(box.get_x(),box.get_Psi(),'.')
    plt.legend(labels,loc=(0.75,0.5));
    plt.xlabel("x");
    plt.ylabel("u");
    plt.title("T = {}".format(T))
    plt.draw();
    plt.savefig("T{}.pdf".format(T));
    

plot_after_time(0.4,['upwind','upwind_2nd_arr','lax_wendroff'])    
plot_after_time(40,['upwind','upwind_2nd_arr','lax_wendroff'])   
