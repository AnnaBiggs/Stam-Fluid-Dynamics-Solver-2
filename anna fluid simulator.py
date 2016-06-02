#Anna density solver

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
#import cProfile
#import re
import advection as ad


#VARIABLE ASSIGNMENT

a0 = np.zeros((100,100))
a = np.zeros((100,100))


N = a0.shape[0]-2   #sets N as 2 less than the number of grid cells per plot side

dt = 1   #defining dt

visc = .1   #defining viscosity constant

bouy_con = 0.5

h0 = np.zeros((100,100))
v0 = np.zeros((100,100))

v = np.zeros((100,100))
h = np.zeros((100,100))


#FUNCTION DEFINITIONS

def draw_array(a1):   #flips i and j in the plot to correspond with desired conventions
    plt.gcf().clear()
    plt.imshow(np.transpose(a1), origin = "lower", cmap="spectral", interpolation = "none")
    plt.colorbar()
    plt.draw()
    plt.show(block=False)

    
def draw_vector_field(h, v, scale):    #h and v are the arrays containing horizontal and vertical components of the vector field
    plt.quiver(h.transpose(),v.transpose(), units = "x", scale = scale )
    plt.show()

# pushes same boundary values into walls cells

def set_bnd0 (N_loc, array):
    array[0,0:(N_loc + 2)] = array[1,0:(N_loc + 2)]
    array[N_loc+1,0:(N_loc + 2)] = array[N_loc,0:(N_loc + 2)]
    array[0:(N_loc + 2),0] = array[0:(N_loc + 2),1]
    array[0:(N_loc + 2),N_loc+1] = array[0:(N_loc + 2),N_loc]
    
#    for i in range (0,N_loc+2):
#        array[0,i] = array[1,i]
#        array[N_loc+1,i] = array[N_loc,i]
#        array[i,0] = array[i,1]
#        array[i,N_loc+1] = array[i,N_loc]

# negates left and right wall values

def set_bnd1 (N_loc, array):
    array[0,0:(N_loc + 2)] = -array[1,0:(N_loc + 2)]
    array[N_loc+1,0:(N_loc + 2)] = -array[N_loc,0:(N_loc + 2)]
    array[0:(N_loc + 2),0] = array[0:(N_loc + 2),1]
    array[0:(N_loc + 2),N_loc+1] = array[0:(N_loc + 2),N_loc]

#    for i in range (0,N_loc+2):
#        array[0,i] = -array[1,i]
#        array[N_loc+1,i] = -array[N_loc,i]
#        array[i,0] = array[i,1]
#        array[i,N_loc+1] = array[i,N_loc]

# negates top and bottom wall values

def set_bnd2 (N_loc, array):
    array[0,0:(N_loc + 2)] = array[1,0:(N_loc + 2)]
    array[N_loc+1,0:(N_loc + 2)] = array[N_loc,0:(N_loc + 2)]
    array[0:(N_loc + 2),0] = -array[0:(N_loc + 2),1]
    array[0:(N_loc + 2),N_loc+1] = -array[0:(N_loc + 2),N_loc]
    
#    for i in range (0,N_loc+2):
#        array[0,i] = array[1,i]
#        array[N_loc+1,i] = array[N_loc,i]
#        array[i,0] = -array[i,1]
#        array[i,N_loc+1] = -array[i,N_loc]
        

#density solver: find the densitites which when diffused backward in time yield the densities we started with

def diffuse(array, array0, diff, N_loc, dt):
    k=0
    while k < 100:
        array[1:(N_loc+1), 1:(N_loc+1)] = (array0[1:(N_loc + 1), 1:(N_loc + 1)] + diff*(array[2:(N_loc + 2), 1:(N_loc + 1)] + array[0:N_loc, 1:(N_loc + 1)] + array[1:(N_loc + 1), 2:(N_loc + 2)] + array[1:(N_loc + 1), 0:N_loc]))/(1+4*diff)
#        for i in range (1, N_loc + 1):
#            for j in range (1, N_loc + 1):
#                array[i,j] = (array0[i,j] + diff*(array[i+1, j] + array[i-1,j] + array[i,j+1] + array[i, j-1]))/(1+4*diff)
        k = k + 1
        set_bnd0(N_loc, array)

    
# interpolate and advect functions compiled in CPython extension module, imported above
        
#def interpolate(array, x, y):
#    i0 = int(x)     #the i & j indices of the cell containing the initial point
#    j0 = int(y)      #remember: int rounds down
#            
#    i1 = i0 + 1
#    j1 = j0 + 1
#        
#    #linear interpolation in the x direction
#    Lx1 = (i1 - x)*array[i0,j0] + (x-i0)*array[i1,j0]
#    Lx2 = (i1 - x)*array[i0,j1] + (x - i0)*array[i1,j1]
#            
#    #Linear interpolation in the y direction
#    interp_val = (j1 - y)*Lx1 + (y-j0)*Lx2
#    return interp_val
    
#takes some decimal indices x and y of an array, returns estimated value of that point based on values of surrounding cells


#def advect(h_vel, v_vel, array, array0, N_loc, dt):   #h_vel is the horizontal velocity field, v_vel is the vertical velocity field, dt is the time step, N_loc is the number of grid cells, array0 is the initial field, and array is the current field
#    for i in range (1, N_loc + 1 ):
#        for j in range (1, N_loc + 1):
#            x_mid = i - 0.5*h_vel[i,j]*dt      #implementing second-order Runge-Kutta method
#            y_mid = j - 0.5*v_vel[i,j]*dt      
#            if x_mid<0.5:
#                x_mid = 0.5
#            if x_mid > N_loc + 0.5:
#                x_mid = N_loc + 0.5
#            if y_mid < 0.5:
#                y_mid = 0.5
#            if y_mid > N_loc + 0.5:
#                y_mid = N_loc + 0.5
#                
#            x_mid = int(x_mid)
#            y_mid = int(y_mid)
#            
#            x = i - h_vel[x_mid, y_mid]*dt    #x and y are the decimal indices of the iniital particle position
#            y = j - v_vel[x_mid, y_mid]*dt
#            
#            if x < 0.5:
#                x = 0.5
#            if x > N_loc + 0.5:
#                x = N_loc + 0.5
#            if y < 0.5:
#                y = 0.5
#            if y > N_loc + 0.5:
#                y = N_loc + 0.5
#            
#            array[i,j] = interpolate(array0, x, y)
            
            
def project(h_vel, v_vel, N_loc):    
    
    #calculate the divergence of the velocity vector field using centered differences
    div = np.zeros((N_loc + 2, N_loc + 2))
    div[1:(N_loc+1),1:(N_loc+1)] = (N_loc/2)*(h_vel[2:(N_loc+2), 1:(N_loc+1)] - h_vel[0:N_loc, 1:(N_loc+1)] + v_vel[1:(N_loc+1), 2:(N_loc + 2)] - v_vel[1:(N_loc+1), 0:N_loc])

#    for i in range(1, N_loc + 1):
#        for j in range (1, N_loc + 1):
#            div[i,j] = (N_loc/2)*(h_vel[i+1,j] - h_vel[i-1,j] + v_vel[i, j+1] - v_vel[i, j-1])
    #divergence left as 0 at boundaries
        
    #solve for pressure
    p = np.zeros((N_loc + 2, N_loc + 2))
    k = 0
    while k < 100:
        p[1:(N_loc+1),1:(N_loc+1)] = 0.25*(p[0:N_loc,1:(N_loc+1)] + p[2:(N_loc + 2),1:(N_loc+1)] + p[1:(N_loc+1),0:N_loc] + p[1:(N_loc+1), 2:(N_loc + 2)] - div[1:(N_loc+1),1:(N_loc+1)]/(N_loc**2))

#        for i in range (1, N_loc + 1):
#            for j in range (1, N_loc + 1):
#                p[i,j] = 0.25*(p[i-1,j] + p[i+1,j] + p[i,j-1] + p[i,j+1] - div[i,j]/(N_loc**2))
        
        k = k + 1
        set_bnd0(N_loc, p)
    
    h_vel[1:(N_loc+1),1:(N_loc+1)] -= (p[2:(N_loc + 2), 1:(N_loc+1)] - p[0:N_loc, 1:(N_loc+1)])*(N_loc/2)
    v_vel[1:(N_loc+1),1:(N_loc+1)] -= (p[1:(N_loc+1), 2:(N_loc + 2)] - p[1:(N_loc+1),0:N_loc])*(N_loc/2)

    #calculate the horizontal and vertical pressure gradients using centered differences, then use it to update the horizontal and vertical velocity fields
#    for i in range (1, N_loc + 1):
#        for j in range (1, N_loc + 1):
#            h_vel[i,j] -= (p[i+1,j] - p[i-1,j])*(N_loc/2)
#            v_vel[i,j] -= (p[i,j+1] - p[i,j-1])*(N_loc/2)
                
    set_bnd1(N_loc, h_vel)
    set_bnd2(N_loc, v_vel)
     
       
    
def set_source(array, N_loc):
    midpnt = int(N_loc/4)
    half_width = int(N_loc/15)
    buffer = int(N/20)
    array[(midpnt-half_width):(midpnt + half_width), buffer:(buffer + 2*half_width)] = 1
        
#    for i in range (midpnt - half_width, midpnt + half_width):
#        for j in range (1, 1 + 2*half_width):
            
def set_vel(h_vel, v_vel, N_loc):
    midpnt = int(N_loc/4)
    half_width = int(N_loc/10)
    h_vel[(midpnt-half_width):(midpnt + half_width), 1:(1 + 2*half_width)] = 2
    v_vel[(midpnt-half_width):(midpnt + half_width), 1:(1 + 2*half_width)] = 1.5
    
    
def bouyancy(array, v_vel, b_num, N_loc):
    v_vel[1:(N_loc + 1),1:(N_loc + 1)] +=b_num*array[1:(N_loc + 1),1:(N_loc + 1)]
#    for i in range (1, N_loc + 1):
#        for j in range (1, N_loc + 1):
#            v_vel[i,j] += b_num*array[i,j]
    set_bnd2(N_loc, v_vel)
    

def animate(h_vel0, v_vel0, h_vel, v_vel, dens_array0, dens_array, N_loc, visc, bouy_con, dt):
    project(h_vel0,v_vel0, N_loc)   #updates h0 and v0
    ad.advect(h_vel0, v_vel0, h_vel, h_vel0, N_loc, dt)   #advects velocity fields along updated h0 v0 velocity field, updating h and v
    set_bnd1(N_loc, h_vel)
    ad.advect(h_vel0, v_vel0, v_vel, v_vel0, N_loc, dt)
    set_bnd2(N_loc, v_vel)
    project(h_vel, v_vel, N_loc)    #makes updated h, v velocity field incompressible again
    diffuse(dens_array, dens_array0, visc, N_loc, dt)     #updates a from a0
    dens_array0, dens_array = dens_array, dens_array0     #a0 now the updated density array
    bouyancy(dens_array0, v_vel, bouy_con, N_loc)    #updates v proportional to the density quantity in the cell
    ad.advect(h_vel, v_vel, dens_array, dens_array0, N_loc, dt) #updates a from a0 by advecting along updated h, v velocity field
    set_bnd0(N_loc, dens_array)
    dens_array0, dens_array = dens_array, dens_array0
    h_vel0, h_vel = h_vel, h_vel0
    v_vel0, v_vel = v_vel, v_vel0
    draw_array(dens_array0)
    set_source(dens_array0, N_loc)
    print ("hello")



#STATEMENTS TO EXECUTE


set_source(a0, N)
set_vel(h0, v0, N)
#draw_array(a0)
#draw_vector_field(h0, v0, 0.5)
#

#while True:
#    project(h0,v0,N)   #updates h0 and v0
#    ad.advect(h0, v0, h, h0, N, dt)   #advects velocity fields along updated h0 v0 velocity field, updating h and v
#    set_bnd1(N, h)
#    ad.advect(h0, v0, v, v0, N, dt)
#    set_bnd2(N, v)
#    project(h, v, N)    #makes updated h, v velocity field incompressible again
#    diffuse(a, a0, visc, N, dt)     #updates a from a0
#    a0, a = a, a0     #a0 now the updated density array
#    bouyancy(a0, v, bouy_con, N)    #updates v proportional to the density quantity in the cell
#    ad.advect(h, v, a, a0, N, dt) #updates a from a0 by advecting along updated h, v velocity field
#    set_bnd0(N, a)
#    a0, a = a, a0
#    h0, h = h, h0
#    v0, v = v, v0
#    draw_array(a0)
#    set_source(a0, N)
#    set_vel(h0, v0, N)




#ANIMATION ATTEMPT


# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, N+2), ylim=(0, N+2))


#draw_array(a0)

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate,fargs = (h0, v0, h, v, a0, a, N, visc, bouy_con, dt), interval = 1)

plt.show()


#PYTHON PROFILER

#
#cProfile.run("animate(h0, v0, h, v, a0, a, N, visc, bouy_con, dt)", filename = "profile_animate", sort = ("cumulative"))
#
##built-in method builtins.exec had much longer cumulative time took way longer than everything else (8.343 seconds, while no other function went above 1). Only called once
#
#import pstats
#
#p = pstats.Stats("profile_animate")
#p.strip_dirs().sort_stats("cumulative").print_stats(10)
