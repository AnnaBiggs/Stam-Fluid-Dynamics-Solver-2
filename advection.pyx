#advection routine in Python for Cython compilation

def interpolate(array, x, y):
    i0 = int(x)     #the i & j indices of the cell containing the initial point
    j0 = int(y)      #remember: int rounds down
            
    i1 = i0 + 1
    j1 = j0 + 1
        
    #linear interpolation in the x direction
    Lx1 = (i1 - x)*array[i0,j0] + (x-i0)*array[i1,j0]
    Lx2 = (i1 - x)*array[i0,j1] + (x - i0)*array[i1,j1]
            
    #Linear interpolation in the y direction
    interp_val = (j1 - y)*Lx1 + (y-j0)*Lx2
    return interp_val
    

def advect(h_vel, v_vel, array, array0, N_loc, dt):   #h_vel is the horizontal velocity field, v_vel is the vertical velocity field, dt is the time step, N_loc is the number of grid cells, array0 is the initial field, and array is the current field
    for i in range (1, N_loc + 1 ):
        for j in range (1, N_loc + 1):
            x_mid = i - 0.5*h_vel[i,j]*dt      #implementing second-order Runge-Kutta method
            y_mid = j - 0.5*v_vel[i,j]*dt      
            if x_mid<0.5:
                x_mid = 0.5
            if x_mid > N_loc + 0.5:
                x_mid = N_loc + 0.5
            if y_mid < 0.5:
                y_mid = 0.5
            if y_mid > N_loc + 0.5:
                y_mid = N_loc + 0.5
                
            x_mid = int(x_mid)
            y_mid = int(y_mid)
            
            x = i - h_vel[x_mid, y_mid]*dt    #x and y are the decimal indices of the iniital particle position
            y = j - v_vel[x_mid, y_mid]*dt
            
            if x < 0.5:
                x = 0.5
            if x > N_loc + 0.5:
                x = N_loc + 0.5
            if y < 0.5:
                y = 0.5
            if y > N_loc + 0.5:
                y = N_loc + 0.5
            
            array[i,j] = interpolate(array0, x, y)