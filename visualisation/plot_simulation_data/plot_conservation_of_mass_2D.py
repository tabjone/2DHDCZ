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
from read_simulation_data.get_num_snaps import get_num_snaps

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


def plot_conservation_of_mass_2D(folder, ax, snap_number, key, **kwargs):
    rho1_initial, unit, _ = get_foreground_variable(folder, 0, 'rho1')
    rho1_initial = rho1_initial * u.Unit(unit)

    num_snaps = get_num_snaps(folder)

    t = np.zeros(num_snaps) * u.s
    t[0] = 0.0 * u.s
    mass_difference = np.zeros(num_snaps) * rho1_initial.unit
    mass_difference[0] = 0.0 * rho1_initial.unit

    for i in range(1, num_snaps):
        rho1, unit, t_current = get_foreground_variable(folder, i, 'rho1')
        rho1 = rho1 * u.Unit(unit)

        mass_difference[i] = np.sum(rho1 - rho1_initial)

        t[i] = t_current * u.s

    t_return = t[snap_number].value
    
    t = t.to(u.h)

    ax.set_title(r"$\int_{V}(\rho_1(t)-\rho_1(0)) dV$")
    ax.set_xlabel("Time [{}]".format(t.unit))
    ax.set_ylabel("Mass difference [{}]".format(mass_difference.unit))

    im = ax.plot(t, mass_difference)

    return im, t_return