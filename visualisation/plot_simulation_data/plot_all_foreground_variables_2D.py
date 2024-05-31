import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter

import astropy.units as u

from .plot_foreground_variable_2D import plot_foreground_variable_2D

def plot_all_foreground_variables_2D(folder, fig, snap_number, **kwargs):
    # Handle kwargs
    if kwargs.get("title_size") == None:
        title_size = 15
    else:
        title_size = kwargs.get("title_size")
    if kwargs.get("thesis") == None:
        gs = gridspec.GridSpec(2, 6, width_ratios=[1, 0.05, 1, 0.05, 1, 0.05], wspace=0.6, hspace=0.3)

        # Order of keys for plotting
        plot_order = ["T1", "rho1", "vy", "p1", "s1", "vz"]

        for idx, key in enumerate(plot_order):
            i, j = divmod(idx, 3)  # Convert 1D index to 2D indices
            # Create the subplot using GridSpec indexing
            ax = fig.add_subplot(gs[i, 2*j])

            # Plot the foreground variable
            im, t = plot_foreground_variable_2D(folder, ax, snap_number, key, **kwargs)
            
            # Create a new axes for the colorbar next to the current subplot
            pos = ax.get_position()
            cax = fig.add_axes([pos.x1 + 0.01, pos.y0, 0.01, pos.height])
            cbar = plt.colorbar(im, cax=cax)

            # Enforce scientific notation
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            cbar.ax.yaxis.set_major_formatter(formatter)

        # Title for the entire plot with the time
        fig.suptitle(f"t={t.value:.1f} [{t.unit.to_string(format='latex_inline')}]", fontsize=title_size, y=0.96)
    else:
        gs = gridspec.GridSpec(3, 4, width_ratios=[1, 0.05, 1, 0.05], wspace=0.1, hspace=0.3)

        # Order of keys for plotting
        plot_order = ["s1", "p1", "rho1", "T1", "vy", "vz"]

        for idx, key in enumerate(plot_order):
            i, j = divmod(idx, 2)  # Convert 1D index to 2D indices
            # Create the subplot using GridSpec indexing
            ax = fig.add_subplot(gs[i, 2*j])

            # Plot the foreground variable
            im, t = plot_foreground_variable_2D(folder, ax, snap_number, key, **kwargs)
            
            # Create a new axes for the colorbar next to the current subplot
            pos = ax.get_position()
            cax = fig.add_axes([pos.x1 + 0.01, pos.y0, 0.01, pos.height])
            cbar = plt.colorbar(im, cax=cax)

            # Enforce scientific notation
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            cbar.ax.yaxis.set_major_formatter(formatter)

        # Title for the entire plot with the time
        fig.suptitle(f"t={t.value:.1f} [{t.unit.to_string(format='latex_inline')}]", fontsize=title_size, y=0.96)