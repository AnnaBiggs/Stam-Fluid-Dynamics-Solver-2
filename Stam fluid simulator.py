#compiled smoke simulation code

import numpy as np
import matplotlib.pyplot as plt

# Creates initial density and velocity arrays

a0 = np.zeros((20,20))    #create iniital density array
a = np.zeros((20,20))        #create next density array

for i in range (8,12):
    for j in range (8,12):
        a0[i,j]=1 


N = a0.shape[0]-2   #sets N as 2 less than the number of grid cells per plot side

dt = .04   #defining dt

visc =  .1   #defining viscosity constant


#h0 = np.zeros((100,100)) #creates initial horizontal and vertical velocity array
#v0 = np.zeros((100,100))
#
#for i in range (45,55):
#    for j in range (45,55):
#        v0[i,j] = 1
#        h0[i,j] = 1
#
#v = np.zeros((100,100))
#h = np.zeros((100,100))

        
# Smaller velocity array for visualization
        
h0 = np.zeros((20,20)) #creates initial horizontal and vertical velocity array
v0 = np.zeros((20,20))

for i in range (8,12):
    for j in range (8,12):
        v0[i,j] = 1
        h0[i,j] = 1

v = np.zeros((20,20))
h = np.zeros((20,20))

#N = 18

#visualizes initial velocity array

#plt.quiver(v0, h0, units = 'x', scale = .5)
#plt.show()
        
#N = 18   #sets N as 2 less than the number of grid cells per plot side


# using "loc" to indicate the local variable version of N

def draw_array(a1):   #flips i and j in the plot to correspond with desired conventions
    plt.imshow(np.transpose(a1), origin = "lower", cmap="spectral", interpolation = "none")
    plt.colorbar()
    plt.show(block=False)
    plt.draw()
    
def draw_vector_field(h, v, scale):    #h and v are the arrays containing horizontal and vertical components of the vector field
    plt.quiver(h.transpose(),v.transpose(), units = "x", scale = scale )
    plt.show()

def set_bnd (N_loc, b, array):
    for i in range (1,N_loc+1):
        if b==1:
            array[0,i]= -array[1,i]
            array[N_loc+1, i] = -array[N_loc, i]
        else:
            array[0,i] = array[1,i]
            array[N_loc+1 ,i] = array[N_loc, i]
            
        if b==2:
            array[i,0] = -array[i,1]
            array[i,N_loc + 1] = -array[i, N_loc]
        else:
            array[i,0] = array[i,1]
            array[i, N_loc+1] = array[i, N_loc]
            
    array[0,0] = 0.5*(array[1,0]+ array[0,1])
    array[0, N_loc + 1] = 0.5*(array[1,N_loc + 1] + array[0, N_loc])
    array[N_loc + 1, 0] = 0.5*(array[N_loc,0] + array[N_loc + 1,1])
    array[N_loc + 1, N_loc + 1] = 0.5*(array[N_loc, N_loc + 1] + array[N_loc + 1, N_loc])



def diffuse(N_loc, dens_array, dens_array0, diff, dt, b):
    c = dt*diff*N_loc*N_loc
    k=0
    while (k<20):
        for i in range (1,N_loc+1):
            for j in range (1,N_loc+1):
                dens_array[i,j]=(dens_array0[i,j]+c*(dens_array[i-1,j] + dens_array[i+1,j] + dens_array[i,j-1] + dens_array[i, j+1]))/(1+4*c)
        k = k + 1
        
    set_bnd(N_loc,b,dens_array)
    
    
        
def advect(N_loc, b, dens_array, dens_array0, h_vel, v_vel, dt):
    dt0 = dt*N_loc
    for i in range (1,N_loc+1):
        for j in range (1, N_loc+1):
            x = i - dt0*v_vel[i,j]
            y = j - dt0*h_vel[i,j]

            
            if (x<0.5):
                x=0.5
                
            if (x > N_loc + 0.5):
                x = N_loc + 0.5
                
            i0 = int(x)
            i1 = i0 + 1
            
            if (y<0.5):
                y=0.5
            if (y > N_loc + 0.5):
                y = N_loc + 0.5
            
            j0 = int(y)
            j1 = j0 +1
            
            s1 = x-i0
            s0 = 1-s1
            t1 = y-j0
            t0 = 1-t1
            
            dens_array[i,j] = s0*(t0*dens_array0[i0,j0] + t1*dens_array0[i0,j1]) + s1*(t0*dens_array0[i1,j0] + t1*dens_array0[i1,j1])
            
            
    set_bnd(N_loc, b, dens_array)
    

def project(N_loc, h_vel, v_vel, p, div):
    h = 1.0/N_loc
    k=0
    for i in range (1,N_loc+1):
        for j in range (1,N_loc+1):
            div[i,j] = -0.5*h*(h_vel[i+1,j] - h_vel[i-1,j] + v_vel[i,j+1] - v_vel[i,j-1])
            p[i,j] = 0
    set_bnd(N_loc, 0, div)
    set_bnd(N_loc, 0, p)
    
    draw_array(div)
    
    while k<20:
        for i in range (1,N_loc+1):
            for j in range (1,N_loc+1):
                p[i,j] = (div[i,j]+p[i-1,j] + p[i+1,j] + p[i,j-1] + p[i,j+1])/4
        set_bnd(N_loc, 0, p)
        k = k+1
        
    draw_array(p)
    
    for i in range (1, N_loc+1):
        for j in range (1, N_loc+1):
            h_vel[i,j] -= 0.5*(p[i+1,j] - p[i-1,j])/h
            v_vel[i,j] -= 0.5*(p[i,j+1] - p[i,j-1])/h
    set_bnd(N_loc, 1, h_vel)
    set_bnd(N_loc, 2, v_vel)



def vel_step(N_loc, h_vel, v_vel, h_vel0, v_vel0, visc, dt):
    
    diffuse(N_loc, h_vel, h_vel0, visc, dt, 1)
    draw_vector_field(h_vel,v_vel,0.5)

    diffuse(N_loc, v_vel, v_vel0, visc, dt, 2)
    draw_vector_field(h_vel,v_vel,0.5)

    project(N_loc, h_vel, v_vel, h_vel0, v_vel0)
    draw_vector_field(h_vel,v_vel,0.5)
    
    h_vel0, h_vel = h_vel, h_vel0
    v_vel0, v_vel = v_vel, v_vel0
    
    advect (N_loc, 1, h_vel, h_vel0, h_vel0, v_vel0, dt)
    draw_vector_field(h_vel,v_vel,0.5)

    advect (N_loc, 2, v_vel, v_vel0, h_vel0, v_vel0, dt)
    draw_vector_field(h_vel,v_vel,0.5)

    project (N_loc, h_vel, v_vel, h_vel0, v_vel0)
    draw_vector_field(h_vel,v_vel,0.5)


#Smoke simulation loop

#while True:
#    vel_step(N, h, v, h0, v0, visc, dt)
#    diffuse (N, a, a0, visc, dt, 0)
#    a0, a = a, a0
#    advect(N, 0, a, a0, h0, v0, dt) #Stam's code has h and v rather than h0 and v0 passed in here, although for whatever reason this seems to make the advection more effective. Why?
#    a0, a = a, a0
#    draw_array(a0)
#    draw_vector_field(h0,v0,0.5)

    
#Visualizing velocity fields at each stage through "vel_step" function

#vel_step (N, h, v, h0, v0, visc, dt)
#draw_vector_field(h0,v0,0.5)
#
#project(N, h0, v0, h, v)
#draw_vector_field(h0,v0,0.5)

while True:
    diffuse(N, a, a0, visc, dt, 0)
    a, a0 = a0, a
    draw_array(a0)