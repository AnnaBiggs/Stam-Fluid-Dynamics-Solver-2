# boundary settings

import numpy as np


# pushes same boundary values into walls cells

def set_bnd0 (N_loc, array):
    for i in range (0,N_loc+2):
        array[0,i] = array[1,i]
        array[N_loc+1,i] = array[N_loc,i]
        array[i,0] = array[i,1]
        array[i,N_loc+1] = array[i,N_loc]

# negates left and right wall values

def set_bnd1 (N_loc, array):
    for i in range (0,N_loc+2):
        array[0,i] = array[1,i]
        array[N_loc+1,i] = array[N_loc,i]
        array[i,0] = -array[i,1]
        array[i,N_loc+1] = -array[i,N_loc]

# negates top and bottom wall values

def set_bnd2 (N_loc, array):
    for i in range (0,N_loc+2):
        array[0,i] = -array[1,i]
        array[N_loc+1,i] = -array[N_loc,i]
        array[i,0] = array[i,1]
        array[i,N_loc+1] = array[i,N_loc]
        
        
# testing boundary conditions
        
a = np.arange(100).reshape(10,10)

print (a)

set_bnd0 (8,a)

print (a)

set_bnd1(8,a)

print (a)

set_bnd2(8,a)

print (a)