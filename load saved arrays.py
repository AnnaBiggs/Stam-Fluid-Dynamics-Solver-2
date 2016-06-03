
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import animation

def draw_array(a1):   #flips i and j in the plot to correspond with desired conventions
    plt.gcf().clear()
    plt.imshow(np.transpose(a1), origin = "lower", cmap="spectral", interpolation = "none")
    plt.colorbar()
    plt.draw()
    plt.show(block=False)

    
def draw_vector_field(h, v, scale):    #h and v are the arrays containing horizontal and vertical components of the vector field
    plt.quiver(h.transpose(),v.transpose(), units = "x", scale = scale )
    plt.show()

def load_dens_array(f, dt, N):
    directory = 'C:\\Users\\16AnnaBB\\Desktop\\SENIOR BACKUP\\Senior Project\\Stam Fluid Dynamics Solver\\Saved Arrays\\dt = ' + str(dt) + ', N = ' + str(N) + '\\DensityArray'
    os.chdir(directory)
    filename = "densarray" + "_" + str(f) + ".npy"
    return np.load(str(filename))


def load_hvel_array(n):
    os.chdir('C:\\Users\\16AnnaBB\\Desktop\\SENIOR BACKUP\\Senior Project\\Stam Fluid Dynamics Solver\\Saved Arrays\\VelocityArray\\Horizontal')
    filename = "H_velarray_" + str(n) + ".npy"
    return np.load(str(filename))

def load_vvel_array(n):
    os.chdir('C:\\Users\\16AnnaBB\\Desktop\\SENIOR BACKUP\\Senior Project\\Stam Fluid Dynamics Solver\\Saved Arrays\\VelocityArray\\Vertical')
    filename = "V_velarray_" + str(n) + ".npy"
    return np.load(str(filename))
    
    
def animate(i):
    #multipy i by a constant to only animate every other, every third, etc. frame
    draw_array(load_dens_array(i, 0.2, 98))
    
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
#ax = plt.axes(xlim=(0, N+2), ylim=(0, N+2))


#draw_array(a0)

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, interval = 0)

plt.show()
    
#fr_num = 0
#fr_total = 73
#
#while fr_num <= fr_total:
#    cur_plot = load_dens_array(fr_num)
#    draw_array(cur_plot)
#    fr_num +=1
#
#fr_num = 0
#
#while fr_num <= fr_total:
#    draw_vector_field(load_hvel_array(fr_num), load_vvel_array(fr_num), 0.5)
#    fr_num +=1