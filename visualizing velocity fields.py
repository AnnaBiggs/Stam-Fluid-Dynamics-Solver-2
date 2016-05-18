#Visualizing velocity fields

import numpy as np
import matplotlib.pyplot as plt


h0 = np.ones((20,20)) #creates initial horizontal and vertical velocity array
v0 = np.ones((20,20))


v = np.zeros((20,20))
h = np.zeros((20,20))

#X,Y = meshgrid(arange(0,10),arange(0,10))

for i in range (8,12):
    for j in range (8,12):
        v[i,j] = 1
        h[i,j] = 1
        
N = 18   #sets N as 2 less than the number of grid cells per plot side

dt = 1   #defining dt

visc = 1   #defining viscosity constant
        


def set_bnd (N_loc, b, array):
    for i in range (1,N_loc):
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
        for i in range (1,N_loc):
            for j in range (1,N_loc):
                dens_array[i,j]=(dens_array0[i,j]+c*(dens_array[i-1,j] + dens_array[i+1,j] + dens_array[i,j-1] + dens_array[i, j+1]))/(1+4*c)
        k = k + 1
        
    set_bnd(N_loc,b,dens_array)
    
    
        
def advect(N_loc, b, dens_array, dens_array0, h_vel, v_vel, dt):
    dt0 = dt*N_loc
    for i in range (1,N_loc):
        for j in range (1, N_loc):
            x = i - dt0*v_vel[i,j]  #could potentially have v_vel and h_vel switched around here...
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
    for i in range (1,N_loc):
        for j in range (1,N_loc):
            div[i,j] = -0.5*h*(h_vel[i+1,j] - h_vel[i-1,j] + v_vel[i,j+1] - v_vel[i,j-1])
            p[i,j] = 0
    set_bnd(N_loc, 0, div)
    set_bnd(N_loc, 0, p)
    
    plt.imshow(div, cmap="gray") # draw the new plot
    plt.show(block=False)
    plt.draw()
    
    while k<20:
        for i in range (1,N_loc):
            for j in range (1,N_loc):
                p[i,j] = (div[i,j] + p[i-1,j] + p[i+1,j] + p[i,j-1] + p[i,j+1])/4
        
        set_bnd(N_loc, 0, p)
        k = k+1
        
    plt.imshow(p, cmap="gray") # draw the new plot
    plt.show(block=False)
    plt.draw()
    
    for i in range (1, N_loc):
        for j in range (1, N_loc):
            h_vel[i,j] -= 0.5*(p[i+1,j] - p[i-1,j])/h
            v_vel[i,j] -= 0.5*(p[i,j+1] - p[i,j-1])/h
    set_bnd(N_loc, 1, h_vel)
    set_bnd(N_loc, 2, v_vel)



def vel_step(N_loc, h_vel, v_vel, h_vel0, v_vel0, visc, dt):
    diffuse(N_loc, h_vel, h_vel0, visc, dt, 1)
    diffuse(N_loc, v_vel, v_vel0, visc, dt, 2)
    project(N_loc, h_vel, v_vel, h_vel0, v_vel0)
    h_vel0 = h_vel
    v_vel0 = v_vel
    advect (N_loc, 1, h_vel, h_vel0, h_vel0, v_vel0, dt)
    advect (N_loc, 2, v_vel, v_vel0, h_vel0, v_vel0, dt)
    project (N_loc, h_vel, v_vel, h_vel0, v_vel0)
    
    
#while True:
#    vel_step(N, h, v, h0, v0, visc, dt)
#    plt.show(block=False)
    
plt.quiver(v, h, units = 'x', scale = .5)
plt.show()

diffuse(18,h,h0,visc,dt,1)
plt.quiver(v, h, units = 'x', scale = .5)
plt.show()

diffuse(18,v,v0,visc,dt,2)
plt.quiver(v, h, units = 'x', scale = .5)
plt.show()


project(18, h, v, h0, v0)
plt.quiver(v, h, units = 'x', scale = .5)
plt.show()
