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
from matplotlib import colors

def plot_foreground_variable_2D(folder, ax, snap_number, key, **kwargs):
        # Get the foreground variable and time
        foreground_variable, unit, t = get_foreground_variable(folder, snap_number, key)
        physical_constants, global_parameters, grid_info = get_info(folder)
        unit = u.Unit(unit)

        K_B = physical_constants['K_B'][0]
        MU = physical_constants['MU'][0]
        R_SUN = physical_constants['R_SUN'][0]
        M_U = physical_constants['M_U'][0]
        GAMMA = global_parameters['GAMMA'][0]

        if kwargs.get('plot_range') != None:
            plot_range = kwargs.get('plot_range')
            r, unit = get_background_variable(folder, "r")
            r_mask = np.where((r > plot_range[0]*R_SUN) & (r < plot_range[1]*R_SUN))
            z0 = plot_range[0]*R_SUN
            z1 = plot_range[1]*R_SUN
        else:
            r_mask = slice(None)
            z0 = global_parameters['R_START'][0]*R_SUN
            z1 = global_parameters['R_END'][0]*R_SUN

        y0 = 0.0
        y1 = global_parameters['Y_SIZE'][0]*R_SUN
        
        dz = grid_info['dz'][0]
        dy = grid_info['dy'][0]

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
                unit = u.km/u.s
                foreground_variable = foreground_variable/1.0e5
                ax.set_title(f"$v_x$ [{unit}]")
                kwargs['quiver_on'] = False
            elif key == 'vz':
                unit = u.km/u.s
                foreground_variable = foreground_variable/1.0e5
                ax.set_title(f"$v_z$ [{unit}]")
                kwargs['quiver_on'] = False
        else:
            if key == 'T1':
                ax.set_title(f"$T_1$ [{unit.to_string(format='latex_inline')}]")
            elif key == 'rho1':
                ax.set_title(f"$\\rho_1$ [{unit.to_string(format='latex_inline')}]")
            elif key == 'p1':
                ax.set_title(f"$p_1$ [{unit.to_string(format='latex_inline')}]")
            elif key == 's1':
                ax.set_title(f"$s_1$ [{unit.to_string(format='latex_inline')}]")
            elif key == 'vy':
                unit = u.km/u.s
                foreground_variable = foreground_variable/1.0e5
                ax.set_title(f"$v_x$ [{unit.to_string(format='latex_inline')}]")
                kwargs['quiver_on'] = False
            elif key == 'vz':
                unit = u.km/u.s
                foreground_variable = foreground_variable/1.0e5
                ax.set_title(f"$v_z$ [{unit.to_string(format='latex_inline')}]")
                kwargs['quiver_on'] = False

        # Handle kwargs
        if kwargs.get("vmin") == None:
            vmin = -np.max(np.abs(foreground_variable[r_mask]))
        else:
            vmin = kwargs.get("vmin")

        if kwargs.get("vmax") == None:
            vmax = np.max(np.abs(foreground_variable[r_mask]))
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
            if kwargs.get('quiver_width') == None:
                quiver_width = 0.5
            else:
                quiver_width = kwargs.get('quiver_width')
            if kwargs.get('quiver_headwidth') == None:
                quiver_headwidth = 5
            else:
                quiver_headwidth = kwargs.get('quiver_headwidth')
            if kwargs.get('quiver_scale') == None:
                quiver_scale = 1.0
            else:
                quiver_scale = kwargs.get('quiver_scale')
            
            vy, _, _ = get_foreground_variable(folder, snap_number, "vy")
            vz, _, _ = get_foreground_variable(folder, snap_number, "vz")
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
            ax.quiver(Y_data, Z_data, vy[::quiver_stride_z, ::quiver_stride_y], vz[::quiver_stride_z, ::quiver_stride_y], linewidths=quiver_width, headwidth=quiver_headwidth, scale=quiver_scale)

        t = t*u.s
        if kwargs.get("t_scale") == None:
            t = t.to(u.hour)
        else:
            t = t.to(kwargs.get("t_scale"))        
        
        ax.set_xlabel("x [Solar radii]", fontsize=font_size)
        ax.set_ylabel("z [Solar radii]", fontsize=font_size)

        if kwargs.get("plot_range") != None:
            ax.set_ylim(plot_range[0], plot_range[1])
            ax.set_xlim(plot_range[2], plot_range[3])
        
        foreground_variable = foreground_variable[r_mask]

        if vmin < 0 and vmax > 0 and norm == TwoSlopeNorm:
            norm = norm(vcenter=0, vmin=vmin, vmax=vmax)
            im =ax.imshow(foreground_variable, origin="lower", extent=extent, aspect=aspect,norm=norm, cmap=cmap)
        elif norm == "log":
            im = ax.imshow(foreground_variable,
                        origin="lower",
                        extent=extent, 
                        aspect=aspect, 
                        norm=colors.SymLogNorm(
                            linthresh=0.03, linscale=0.03,
                            vmin=-1.0, vmax=1.0, base=10),
                        cmap=cmap)
        elif norm == "power":
            if kwargs.get("gamma") == None:
                gamma = 0.5
            else:
                gamma = kwargs.get("gamma")
            im = ax.imshow(foreground_variable,
                        origin="lower",
                        extent=extent, 
                        aspect=aspect, 
                        norm=colors.PowerNorm(gamma=gamma),
                        cmap=cmap)

        else:
            im =ax.imshow(foreground_variable, origin="lower", extent=extent, aspect=aspect,vmin=vmin, vmax=vmax, cmap=cmap)

        return im, t