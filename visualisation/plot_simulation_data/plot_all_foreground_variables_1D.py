import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import astropy.units as u

from .plot_foreground_variable_1D import plot_foreground_variable_1D

def plot_all_foreground_variables_1D(folder, fig, snap_number, **kwargs):
    # Handle kwargs
    if kwargs.get("title_size") == None:
        title_size = 15
    else:
        title_size = kwargs.get("title_size")

    gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 1], wspace=0.6, hspace=0.3)

    standard_order = ["T1", "rho1", "p1", "s1", "vz"]
    # Order of keys for plotting
    if kwargs.get("plot_order") == None:
        plot_order = standard_order
    else:
        plot_order = kwargs.get("plot_order")

    for idx, key in enumerate(plot_order):
        i, j = divmod(idx, 3)  # Get the row and column index
        # Create the subplot using GridSpec indexing
        ax = fig.add_subplot(gs[i, j])

        # Plot the foreground variable
        if key in standard_order:
            im, t = plot_foreground_variable_1D(folder, ax, snap_number, key, **kwargs)
        elif key in kwargs.get('custom_plots', {}):
            # It's a custom plot, call the corresponding function
            custom_plot_func = kwargs['custom_plots'][key]
            im, t = custom_plot_func(folder, ax, snap_number, key, **kwargs)

    # Title for the entire plot with the time
    fig.suptitle(f"t={t.value:.1f} [{t.unit.to_string(format='latex_inline')}]", fontsize=title_size, y=0.96)
