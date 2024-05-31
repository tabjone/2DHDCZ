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

def plot_foreground_variable_cross_section_3D(folder, ax, snap_number, key, direction='vertical_yz', cross_section_placement=0.5, **kwargs):
    
    physical_constants, global_parameters, grid_info = get_info(folder)

    K_B = physical_constants['K_B'][0]
    MU = physical_constants['MU'][0]
    R_SUN = physical_constants['R_SUN'][0]
    M_U = physical_constants['M_U'][0]
    GAMMA = global_parameters['GAMMA'][0]

    x0 = 0.0
    x1 = global_parameters['X_SIZE'][0]*R_SUN
    y0 = 0.0
    y1 = global_parameters['Y_SIZE'][0]*R_SUN
    z0 = global_parameters['R_START'][0]*R_SUN
    z1 = global_parameters['R_END'][0]*R_SUN

    nx = global_parameters['NX'][0]
    ny = global_parameters['NY'][0]
    nz = global_parameters['NZ'][0]

    dz = grid_info['dz'][0]
    dy = grid_info['dy'][0]
    dx = grid_info['dx'][0]

    foreground_variable, unit, t = get_foreground_variable(folder, snap_number, key)

    if kwargs.get('cross_section_index') != None:
        cross_section_index = kwargs.get('cross_section_index')
        if direction == 'vertical_xz':
            cross_section_placement_actual = (cross_section_index*dy)/R_SUN
        elif direction == 'vertical_yz':
            cross_section_placement_actual = (cross_section_index*dx)/R_SUN
        elif direction == 'horizontal_xy':
            cross_section_placement_actual = (cross_section_index*dz + z0)/R_SUN
    else:
        if direction == 'vertical_xz':
            cross_section_index = int((ny-1)*cross_section_placement)
            cross_section_placement_actual = (cross_section_index*dy)/R_SUN
        elif direction == 'vertical_yz':
            cross_section_index = int((nx-1)*cross_section_placement)
            cross_section_placement_actual = (cross_section_index*dx)/R_SUN
        elif direction == 'horizontal_xy':
            cross_section_index = int((nz-1)*cross_section_placement)
            cross_section_placement_actual = (cross_section_index*dz + z0)/R_SUN
        else:
            raise ValueError("Invalid direction. Choose either 'vertical_xz', 'vertical_yz' or 'horizontal_xy'")
    
    if direction == 'vertical_xz':
        foreground_variable = foreground_variable[:,cross_section_index,:]
    if direction == 'vertical_yz':
        foreground_variable = foreground_variable[:,:,cross_section_index]
    if direction == 'horizontal_xy':
        foreground_variable = foreground_variable[cross_section_index,:,:]

    # Normalize the foreground variable
    if kwargs.get("normalize") == True:
        if key == 'T1':
            background_variable, _ = get_background_variable(folder, "T0")
            foreground_variable = foreground_variable/background_variable[:,np.newaxis]
            ax.set_title("$T_1/T_0$")
        elif key == 'rho1':
            background_variable, _ = get_background_variable(folder, "rho0")
            foreground_variable = foreground_variable/background_variable[:,np.newaxis]
            ax.set_title(r"$\rho_1/\rho_0$")
        elif key == 'p1':
            background_variable, _ = get_background_variable(folder, "p0")
            foreground_variable = foreground_variable/background_variable[:,np.newaxis]
            ax.set_title("$p_1/p_0$")
        elif key == 's1':
            c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA)
            foreground_variable = foreground_variable/c_p
            ax.set_title("$s_1/c_p$")
        elif key == 'vy':
            ax.set_title(f"$v_y$ [{unit}]")
            kwargs['quiver_on'] = False
        elif key == 'vz':
            ax.set_title(f"$v_z$ [{unit}]")
            kwargs['quiver_on'] = False
    else:
        if key == 'T1':
            ax.set_title(f"$T_1$ [{unit}]")
        elif key == 'rho1':
            ax.set_title(f"$\\rho_1$ [{unit}]")
        elif key == 'p1':
            ax.set_title(f"$p_1$ [{unit}]")
        elif key == 's1':
            ax.set_title(f"$s_1$ [{unit}]")
        elif key == 'vy':
            ax.set_title(f"$v_y$ [{unit}]")
            kwargs['quiver_on'] = False
        elif key == 'vz':
            ax.set_title(f"$v_z$ [{unit}]")
            kwargs['quiver_on'] = False

        # Handle kwargs
        if kwargs.get("vmin") == None:
            vmin = -np.max(foreground_variable)
        else:
            vmin = kwargs.get("vmin")

        if kwargs.get("vmax") == None:
            vmax = np.max(foreground_variable)
        else:
            vmax = kwargs.get("vmax")

        if kwargs.get("cmap") == None:
            cmap = "RdBu"
        else:
            cmap = kwargs.get("cmap")

        if kwargs.get("norm") == None:
            norm = None
        else:
            norm = kwargs.get("norm")

        if kwargs.get("aspect") == None:
            aspect = (y1-y0)/(z1-z0) * 1.2
        else:
            aspect = kwargs.get("aspect")

        if kwargs.get("font_size") == None:
            font_size = 13
        else:
            font_size = kwargs.get("font_size")

        if kwargs.get("extent") == None:
            extent = np.array([y0,y1,z0,z1])/R_SUN
        else:
            extent = kwargs.get("extent")

        if kwargs.get("quiver_on") == True:
            vy, _, _ = get_foreground_variable(folder, snap_number, "vy")
            vz, _, _ = get_foreground_variable(folder, snap_number, "vz")
            if direction == 'vertical_xz':
                vy = vy[:,cross_section_index,:]
                vz = vz[:,cross_section_index,:]
            if direction == 'vertical_yz':
                vy = vy[:,:,cross_section_index]
                vz = vz[:,:,cross_section_index]
            if direction == 'horizontal_xy':
                vy = vy[cross_section_index,:,:]
                vz = vz[cross_section_index,:,:]

            z_shape, y_shape = vy.shape

            if kwargs.get("num_quivers") == None:
                num_quivers = 10
            else:
                num_quivers = kwargs.get("num_quivers")

            quiver_stride_y = int(y_shape//num_quivers)
            quiver_stride_z = int(z_shape//num_quivers)

            Z, Y = np.mgrid[0:z_shape:quiver_stride_z, 0:y_shape:quiver_stride_y]
            # Convert pixel indices to data coordinates
            Y_data = (Y * dy)/R_SUN
            Z_data = (Z * dz)/R_SUN + z0/R_SUN
            ax.quiver(Y_data, Z_data, vy[::quiver_stride_z, ::quiver_stride_y], vz[::quiver_stride_z, ::quiver_stride_y])

        t = t*u.s
        if kwargs.get("t_scale") == None:
            t = t.to(u.hour)
        else:
            t = t.to(kwargs.get("t_scale"))        
        
        ax.set_xlabel("x [Solar radii]", fontsize=font_size)
        ax.set_ylabel("z [Solar radii]", fontsize=font_size)

        if vmin < 0 and vmax > 0 and norm == TwoSlopeNorm:
            norm = norm(vcenter=0, vmin=vmin, vmax=vmax)
            im =ax.imshow(foreground_variable, origin="lower", extent=extent, aspect=aspect,norm=norm, cmap=cmap)
        else:
            im =ax.imshow(foreground_variable, origin="lower", extent=extent, aspect=aspect,vmin=vmin, vmax=vmax, cmap=cmap)

        return im, t, cross_section_placement_actual, cross_section_index