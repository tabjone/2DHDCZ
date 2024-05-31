from .plot_all_foreground_variables_cross_section_3D import plot_all_foreground_variables_cross_section_3D

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

from read_simulation_data.get_foreground_variable import get_foreground_variable

def animate_all_foreground_variables_across_axis_3D(folder, snap_number, save=False, save_name=None, fps=4, save_interval=1, anim_interval=250, axis='x', **kwargs):
    
    foreground_variable, _, _ = get_foreground_variable(folder, 0, 'p1')
    
    if axis == 'x':
        direction = 'vertical_yz'
        axis_size = foreground_variable.shape[2]
    elif axis == 'y':
        direction = 'vertical_xz'
        axis_size = foreground_variable.shape[1]
    elif axis == 'z':
        direction = 'horizontal_xy'
        axis_size = foreground_variable.shape[0]
    else:
        raise ValueError('axis must be x, y or z')

    fig = plt.figure(figsize=(16, 9))  # Create a figure

    def init_animation():
        kwargs['cross_section_index'] = 0
        plot_all_foreground_variables_cross_section_3D(folder, fig, snap_number, direction=direction, **kwargs)

    def update_animation(index_nr):
        """Update the plots for each frame."""
        for ax in fig.get_axes():
            ax.clear()  # Clear previous data
        fig.clear()
        kwargs['cross_section_index'] = index_nr
        plot_all_foreground_variables_cross_section_3D(folder, fig, snap_number, direction=direction, **kwargs)

    anim = FuncAnimation(fig, update_animation, interval=anim_interval, frames=range(0,axis_size,save_interval), init_func=init_animation)

    if save:
        anim.save(save_name, writer='ffmpeg', fps=fps, extra_args=['-vcodec', 'libx264'])
    else:
        plt.show()