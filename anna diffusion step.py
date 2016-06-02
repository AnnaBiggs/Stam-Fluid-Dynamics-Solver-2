#Anna diffusion step

#find the densitites which when diffused backward in time yield the densities we started with

#How: given some cell with some density, I wrote an equation for the cell's previous density based on 
#backward exchange of densities with its four direct neighbors.
#Then, I solved for the current cell density, getting an expression in terms of the cell's previous density
#and the current density of its four neighbors


import numpy as np
import matplotlib.pyplot as plt

def diffuse(array, array0, diff, N_loc, dt):
    k=0
    
    while k < 20:
        for i in range (1, N_loc + 1):
            for j in range (1, N_loc + 1):
                array[i,j] = 0.25*(array0[i,j]/(diff*dt) + array[i+1,j] + array[i-1,j] + array[i, j+1] + array[i, j-1])
        k = k + 1
                
#Testing diffusion step:

def draw_array(a1):   #flips i and j in the plot to correspond with desired conventions
    plt.imshow(np.transpose(a1), origin = "lower", cmap="spectral", interpolation = "none")
    plt.colorbar()
    plt.show(block=False)
    plt.draw()
    

a0 = np.zeros((100,100))    #create iniital density array
a = np.zeros((100,100))        #create next density array

for i in range (45,55):
    for j in range (45,55):
        a0[i,j]=1 


N = a0.shape[0]-2   #sets N as 2 less than the number of grid cells per plot side

dt = .1   #defining dt

visc =  .1   #defining viscosity constant

draw_array(a0)

while True:
    diffuse(a,a0,visc,N,dt)   #update a with a0
    a, a0 = a0, a              #call the updated array a0
    draw_array(a0)             #draw a0


#Notes: 
# 1) likely a bad idea to divide by "diff" in case the viscosity contant passed in
#is zero...but how to change that? can I just avoid making the viscosity zero?
# 2) The blob is diffusing into the top rightmost corner...why?