from .plot_all_foreground_variables_2D import plot_all_foreground_variables_2D

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

import sys
import os

# Get the directory of the script file
script_dir = os.path.dirname(os.path.realpath(__file__))

# Add the parent directory to sys.path
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from read_simulation_data.get_num_snaps import get_num_snaps

def animate_all_foreground_variables_2D(folder, save=False, save_name=None, fps=4, save_interval=10, anim_interval=250, **kwargs):

    fig = plt.figure(figsize=(16, 9))  # Create a figure
    num_snaps = get_num_snaps(folder)
    if kwargs.get('snap_range') is not None:
        snap_range = kwargs.get('snap_range')
        init_snap = snap_range[0]
    else:
        snap_range = range(0, num_snaps, save_interval)
        init_snap = 0

    def init_animation():
        plot_all_foreground_variables_2D(folder, fig, init_snap, **kwargs)

    def update_animation(snap_number):
        """Update the plots for each frame."""
        for ax in fig.get_axes():
            ax.clear()  # Clear previous data
        fig.clear()
        plot_all_foreground_variables_2D(folder, fig, snap_number, **kwargs)
    
    anim = FuncAnimation(fig, update_animation, interval=anim_interval, frames=snap_range, init_func=init_animation)

    if save:
        anim.save(save_name, writer='ffmpeg', fps=fps, extra_args=['-vcodec', 'libx264'])
    else:
        plt.show()