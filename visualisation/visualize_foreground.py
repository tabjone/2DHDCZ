import os
import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from IPython.display import HTML

from astropy import units as u

# Format colorbar ticks in scientific notation
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)

R_sun = 6.957e10

def read_fg(file_path):
    with h5py.File(file_path, 'r') as f:
        T1 = np.array(f['/T1'])
        rho1 = np.array(f['/rho1'])
        p1 = np.array(f['/p1'])
        s1 = np.array(f['/s1'])
        vx = np.array(f['/vx'])
        vz = np.array(f['/vz'])
        info = f['info'][0].decode("utf-8")

    return T1, rho1, p1, s1, vx, vz, info

def read_bg(file_path):
    with h5py.File(file_path, 'r') as f:
        r = np.array(f['/r'])
        T0 = np.array(f['/T0'])
        rho0 = np.array(f['/rho0'])
        p0 = np.array(f['/p0'])
        g = np.array(f['/g'])
        grad_s0 = np.array(f['/grad_s0'])
        
    return r, T0, rho0, p0, g, grad_s0

def get_info(info):
    t = eval(info.split("\n")[0].split(":")[1])
    nx = eval(info.split("\n")[1].split(":")[1])
    nz = eval(info.split("\n")[2].split(":")[1])
    nz_ghost = eval(info.split("\n")[3].split(":")[1])
    dx = eval(info.split("\n")[4].split(":")[1])
    dz = eval(info.split("\n")[5].split(":")[1])
    return t, nx, nz, nz_ghost, dx, dz

def format_title(key):
    # Dictionary for special replacements
    special = {"rho": r"\rho", "p": "p", "T": "T", "v": "v"}

    # Using the last character as the subscript
    prefix = key[:-1]
    suffix = key[-1]

    # Replace with special character if available
    prefix = special.get(prefix, prefix)

    return "${}_{}$".format(prefix, suffix)

class Visualize_Foreground:
    def __init__(self, folder):
        self.folder = folder + "snap{}.h5"
        self.num_snaps = len(os.listdir(folder))-1
        
        self.set_plot_params()
        
    def set_plot_params(self):
        self.quiver_stride = 10
        self.font_size = 13
        self.title_size = 15
        
        max_values = [0] * 6
        min_values = [0] * 6
        for i in range(self.num_snaps):
            T1, rho1, p1, s1, vx, vz, info = read_fg(self.folder.format(i))
            
            max_values = np.maximum(max_values, [T1.max(), rho1.max(), p1.max(), s1.max(), vx.max(), vz.max()])
            min_values = np.minimum(min_values, [T1.min(), rho1.min(), p1.min(), s1.min(), vx.min(), vz.min()])
        
        T1, rho1, p1, s1, vx, vz, info = read_fg(self.folder.format(self.num_snaps-1))
        t, nx, nz, nz_ghost, dx, dz = get_info(info)
        
        """
        THIS NEEDS TO BE IN THE INFO FILE
        """
        self.x_0 = 0
        self.x_1 = dx*nx/R_sun
        self.z_0 = 0.6
        self.z_1 = 0.6 + nz*dz/R_sun
        
        self.aspect = (self.x_1-self.x_0)/(self.z_1-self.z_0) * 1.2
        
        t_end = t*u.s
        if t>1e4:
            self.t_end = t_end.to("h")
        if t>1e5:
            self.t_end = t_end.to("day")
        
        self.quiver_scale = 4*max(max(max_values[4], max_values[5]), max(np.abs(min_values[4]), np.abs(min_values[5])))
        
        # holding touple(var_name, quiver true/false, vmin, vmax)
        self.plot_params = {
            "T1": (True, min_values[0], max_values[0], "Temperature [K]"),
            "rho1": (True, min_values[1], max_values[1], r"Density [g/cm$^3$]"),
            "p1": (True, min_values[2], max_values[2], r"Pressure [dyn/cm$^2$]"),
            "s1": (True, min_values[3], max_values[3], "Entropy [erg/K]"),
            "vx": (False, min_values[4], max_values[4], "x-velocity [cm/s]"),
            "vz": (False, min_values[5], max_values[5], "z-velocity [cm/s]")
        }     

    def plot(self, fig, ax, snap_nr, key):
        T1, rho1, p1, s1, vx, vz, info = read_fg(self.folder.format(snap_nr))
        d = {"T1": T1, "rho1": rho1, "p1": p1, "s1": s1, "vx":vx, "vz":vz}
        t, nx, nz, nz_ghost, dx, dz = get_info(info)
    
        d = d[key][nz_ghost:-1-nz_ghost+1,:]
        if self.plot_params[key][0]:
            vx = vx[nz_ghost:-1-nz_ghost+1,:]
            vz = vz[nz_ghost:-1-nz_ghost+1,:]

            Y, X = np.mgrid[0:d.shape[0]:self.quiver_stride, 0:d.shape[1]:self.quiver_stride]
            # Convert pixel indices to data coordinates
            X_data = (X * dx)/R_sun
            Y_data = (Y * dz)/R_sun+self.z_0
            ax.quiver(X_data, Y_data, vx[::self.quiver_stride, ::self.quiver_stride], vz[::self.quiver_stride, ::self.quiver_stride], scale=self.quiver_scale)
        
        t = t*u.s
        t = t.to(self.t_end.unit)
        
        ax.set_title(format_title(key))
        ax.set_xlabel("z [Solar radii]", fontsize=self.font_size)
        ax.set_ylabel("x [Solar radii]", fontsize=self.font_size)
        
        vmin = self.plot_params[key][1]
        vmax = self.plot_params[key][2]
        im =ax.imshow(d, origin="lower", extent=[self.x_0,self.x_1,self.z_0,self.z_1], aspect=self.aspect,vmin=vmin, vmax=vmax)
        return im, t

    def plot_all(self, fig, snap_nr):
        gs = gridspec.GridSpec(2, 6, width_ratios=[1, 0.05, 1, 0.05, 1, 0.05], wspace=0.6, hspace=0.3)

        # Order of keys for plotting
        plot_order = ["T1", "rho1", "vx", "p1", "s1", "vz"]

        for idx, key in enumerate(plot_order):
            i, j = divmod(idx, 3)  # Convert 1D index to 2D indices
            # Create the subplot using GridSpec indexing
            ax = fig.add_subplot(gs[i, 2*j])
            
            im, t = self.plot(fig, ax, snap_nr, key)
            
            # Create a new axes for the colorbar next to the current subplot
            pos = ax.get_position()
            cax = fig.add_axes([pos.x1 + 0.01, pos.y0, 0.01, pos.height])
            cbar = plt.colorbar(im, cax=cax)

            # Set label for the colorbar
            cbar.set_label(self.plot_params[key][3], fontsize=self.font_size)

            # Enforce scientific notation
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            cbar.ax.yaxis.set_major_formatter(formatter)

        # Title for the entire plot with the time
        fig.suptitle(f"t={t.value:.1f} [{t.unit.to_string(format='latex_inline')}]", fontsize=self.title_size, y=0.96)

    def save_image_all(self, snap_nr, filename):
        fig = plt.figure(figsize=(16, 9))
        self.plot_all(fig, snap_nr)
        fig.savefig(filename, dpi=80, bbox_inches='tight')

        
    def init_animation(self):
        self.plot_all(self.fig, 0)
        
    def update_animation(self, snap_nr):
        """Update the plots for each frame."""
        for ax in self.fig.get_axes():
            ax.clear()  # Clear previous data
        self.fig.clear()
        self.plot_all(self.fig, snap_nr)  # Plot the new snapshot
            
    def animate(self, save=False, save_name=None, fps=4):
        self.fig = plt.figure(figsize=(16, 9))  # Create a figure
        #self.plot_all(self.fig,self.num_snaps-1)
        anim = FuncAnimation(self.fig, self.update_animation, interval=250, frames=self.num_snaps, init_func=self.init_animation)
        if save:
            anim.save(save_name, writer='ffmpeg', fps=fps, extra_args=['-vcodec', 'libx264'])
        else:
            plt.show()

vf = Visualize_Foreground("../data/velocity_right/")
vf.animate(save=True, save_name="velocity_right_1000.mp4")