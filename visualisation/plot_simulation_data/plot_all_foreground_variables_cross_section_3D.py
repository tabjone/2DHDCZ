import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter

import astropy.units as u

from .plot_foreground_variable_cross_section_3D import plot_foreground_variable_cross_section_3D

def plot_all_foreground_variables_cross_section_3D(folder, fig, snap_number, direction='vertical_yz', cross_section_placement=0.5, **kwargs):
    # Handle kwargs
    if kwargs.get("title_size") == None:
        title_size = 15
    else:
        title_size = kwargs.get("title_size")

    gs = gridspec.GridSpec(2, 6, width_ratios=[1, 0.05, 1, 0.05, 1, 0.05], wspace=0.6, hspace=0.3)

    if kwargs.get("plot_order") == None:
        # Order of keys for plotting
        plot_order = ["T1", "rho1", "vy", "p1", "s1", "vz"]
    else:
        plot_order = kwargs.get("plot_order")

    for idx, key in enumerate(plot_order):
        i, j = divmod(idx, 3)  # Convert 1D index to 2D indices
        # Create the subplot using GridSpec indexing
        ax = fig.add_subplot(gs[i, 2*j])

        # Plot the foreground variable
        im, t, cross_section_placement_actual, cross_section_index = plot_foreground_variable_cross_section_3D(folder, ax, snap_number, key, direction=direction, cross_section_placement=cross_section_placement, **kwargs)
        
        # Create a new axes for the colorbar next to the current subplot
        pos = ax.get_position()
        cax = fig.add_axes([pos.x1 + 0.01, pos.y0, 0.01, pos.height])
        cbar = plt.colorbar(im, cax=cax)

        # Enforce scientific notation
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        cbar.ax.yaxis.set_major_formatter(formatter)

    # Find values for the cross section placement
    if direction == 'vertical_yz':
        cross_section_title = f'x={cross_section_placement_actual:.2f}R$_\odot$'
    elif direction == 'vertical_xz':
        cross_section_title = f'y={cross_section_placement_actual:.2f}R$_\odot$'
    elif direction == 'horizontal_xy':
        cross_section_title = f'z={cross_section_placement_actual:.2f}R$_\odot$'

    # Title for the entire plot with the time
    fig.suptitle(f"t={t.value:.1f} [{t.unit.to_string(format='latex_inline')}]," + cross_section_title + f' Index={cross_section_index}', fontsize=title_size, y=0.96)