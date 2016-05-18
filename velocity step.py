#velocity step

import numpy as np
import matplotlib.pyplot as plt

def vel_step(N_loc, h_vel_loc, v_vel_loc, h_vel_loc0, v_vel_loc0, visc, dt):
    diffuse(N_loc, h_vel_loc, h_loc0, visc, dt, 1)
    diffuse(N_loc, v_vel_loc, v_vel_loc0, visk, dt, 1)
    project(N_loc, h_vel_loc, v_vel_loc, h_vel_loc0, v_vel_loc0)
    SWAP (N_loc, h_vel_loc0, h_vel_loc)
    SWAP (N_loc, v_vel_loc0, v_vel_loc)
    advect (N_loc, h_vel_loc0, v_vel_loc0, h_vel_loc, h_vel_loc0, dt, 1)
    advect (N_loc, h_vel_loc0, v_vel_loc0, v_vel_loc, v_vel_loc0, dt, 2)
    project (N_loc, h_vel_loc, v_vel_loc, h_vel_loc0, v_vel_loc0)