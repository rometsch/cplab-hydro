# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 16:46:54 2016

@author: thomas

This script is to test how numpy arrays behave when altered inside a function.
Since python objects are passed by reference, the expectation is that
numpy arrays behave just the same way, but its unclear how operations
like += behave.
"""

import numpy as np;

def add_to(arr):
    '''
    Add something to the array using the += operator.
    '''
    arr += np.array([1,1,2,2]);

def mult_by(arr):
    arr *= np.array([2,2,3,3]);
    
def set_to(arr):
    ''' use range operator to set new values '''
    arr[:] = np.array([1,2,3,4]);

outer_array = np.array([2.,3.,4.,5.]);
print("test +=");
print('before',outer_array);
add_to(outer_array);
print("after",outer_array)

print("");
print("test *=");
print("before",outer_array);
mult_by(outer_array);
print("after",outer_array);

print("");
print("test setting new values");
print("before",outer_array);
set_to(outer_array);
print("after",outer_array);