import sys
import os

# Get the directory of the script file
script_dir = os.path.dirname(os.path.realpath(__file__))

# Add the parent directory to sys.path
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from read_simulation_data.get_foreground_variable import get_foreground_variable
from read_simulation_data.get_background_variable import get_background_variable
from read_simulation_data.get_info import get_info

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

def plot_foreground_variable_1D(folder, ax, snap_number, key, **kwargs):
        # Get the foreground variable and time
        foreground_variable, unit, t = get_foreground_variable(folder, snap_number, key)
        radius, _ = get_background_variable(folder, "r")
        physical_constants, global_parameters, grid_info = get_info(folder)

        K_B = physical_constants['K_B'][0]
        MU = physical_constants['MU'][0]
        R_SUN = physical_constants['R_SUN'][0]
        M_U = physical_constants['M_U'][0]
        GAMMA = global_parameters['GAMMA'][0]

        z0 = global_parameters['R_START'][0]*R_SUN
        z1 = global_parameters['R_END'][0]*R_SUN

        dz = grid_info['dz'][0]

        # Normalize the foreground variable
        if kwargs.get("normalize") == True:
            if key == 'T1':
                background_variable, _ = get_background_variable(folder, "T0")
                foreground_variable = foreground_variable/background_variable
                ax.set_ylabel("$T_1/T_0$")
            elif key == 'rho1':
                background_variable, _ = get_background_variable(folder, "rho0")
                foreground_variable = foreground_variable/background_variable
                ax.set_ylabel(r"$\rho_1/\rho_0$")
            elif key == 'p1':
                background_variable, _ = get_background_variable(folder, "p0")
                foreground_variable = foreground_variable/background_variable
                ax.set_ylabel("$p_1/p_0$")
            elif key == 's1':
                c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA)
                foreground_variable = foreground_variable/c_p
                ax.set_ylabel("$s_1/c_p$")
            elif key == 'vy':
                ax.set_ylabel(f"$v_x$ [{unit}]")
            elif key == 'vz':
                ax.set_ylabel(f"$v_z$ [{unit}]")
        else:
            if key == 'T1':
                ax.set_ylabel(f"$T_1$ [{unit}]")
            elif key == 'rho1':
                ax.set_ylabel(f"$\\rho_1$ [{unit}]")
            elif key == 'p1':
                ax.set_ylabel(f"$p_1$ [{unit}]")
            elif key == 's1':
                ax.set_ylabel(f"$s_1$ [{unit}]")
            elif key == 'vy':
                ax.set_ylabel(f"$v_x$ [{unit}]")
            elif key == 'vz':
                ax.set_ylabel(f"$v_z$ [{unit}]")

        # Handle kwargs
        if kwargs.get("font_size") == None:
            font_size = 13
        else:
            font_size = kwargs.get("font_size")
        if kwargs.get('grid_on') == None:
            ax.grid(True)
        else:
            ax.grid(kwargs.get('grid_on'))

        t = t*u.s
        if kwargs.get("t_scale") == None:
            t = t.to(u.hour)
        else:
            t = t.to(kwargs.get("t_scale"))        
        
        ax.set_xlabel("z [Solar radii]", fontsize=font_size)

        im = plt.plot(radius/R_SUN, foreground_variable)

        return im, t