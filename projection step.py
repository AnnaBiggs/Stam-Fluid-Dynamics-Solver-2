#projection step

import numpy as np
import matplotlib.pyplot as plt



def project(N_loc, h_vel_loc, v_vel_loc, p, div):
    h = 1.0/N_loc
    k=0
    for i in range (1,N_loc):
        for j in range (1,N_loc):
            div[i,j] = -0.5*h*(h_vel_loc[i+1,j] - h_vel_loc[i-1,j] + v_vel_loc[i,j+1] - v_vel_loc[i,j-1])
            p[i,j] = 0
    set_bnd(N, 0, div)
    set_bnd(N, 0, p)
    
    while k<20:
        for i in range (1,N_loc):
            for j in range (1,N_loc):
                p[i,j] = (div[i,j]+p[i-1,j] + p[i+1,j] + p[i,j-1] + p[i,j+1])/4
        set_bnd(N, 0, p)
        k = k+1
    
    for i in range (1, N_loc):
        for j in range (1,N_loc):
            h_vel_loc[i,j] -= 0.5*(p[i+1,j] - p[i-1,j])/h
            v_vel_loc[i,j] -= 0.5*(p[i,j+1] - p[i,j-1])/h
    set_bnd(N, 1, h_vel_loc)
    set_bnd(N, 2, v_vel_loc)