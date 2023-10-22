import os
import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from IPython.display import HTML
from matplotlib.colors import TwoSlopeNorm, NoNorm

from astropy import units as u

# Format colorbar ticks in scientific notation
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)

R_sun = 6.957e10

def read_fg(file_path):
    # NOW THIS IS A LITTLE BIT HACKY TO MAKE NEW DATA WORK WITH OLD CODE
    with h5py.File(file_path, 'r') as f:
        T1 = np.array(f['variables/T1'])
        rho1 = np.array(f['variables/rho1'])
        p1 = np.array(f['variables/p1'])
        s1 = np.array(f['variables/s1'])
        vy = np.array(f['variables/vy'])
        vz = np.array(f['variables/vz'])
        
        t = np.array(f['grid_info/t'])
        nx = np.array(f['grid_info/ny'])
        nz = np.array(f['grid_info/nz'])
        nz_full = np.array(f['grid_info/nz_full'])
        nz_ghost = np.array(f['grid_info/nz_ghost'])
        dx = np.array(f['grid_info/dy'])
        dz = np.array(f['grid_info/dz'])
        z0 = np.array(f['grid_info/z0'])
        z1 = np.array(f['grid_info/z1'])
        x0 = np.array(f['grid_info/y0'])
        x1 = np.array(f['grid_info/y1'])
        
        variables = {"T1": T1, "rho1": rho1, "p1": p1, "s1": s1, "vx": vy, "vz": vz}
        info = {"t": t, "nx": nx, "nz": nz, "nz_full": nz_full, "nz_ghost": nz_ghost, "dx": dx, "dz": dz, "z0": z0, "z1": z1, "x0": x0, "x1": x1}
    return variables, info

def read_bg(file_path):
    with h5py.File(file_path, 'r') as f:
        r = np.array(f['variables/r'])
        T0 = np.array(f['variables/T0'])
        rho0 = np.array(f['variables/rho0'])
        p0 = np.array(f['variables/p0'])
        g = np.array(f['variables/g'])
        grad_s0 = np.array(f['variables/grad_s0'])

        nz = np.array(f['grid_info/nz'])
        nz_ghost = np.array(f['grid_info/nz_ghost'])
        nz_full = np.array(f['grid_info/nz_full'])
        dz = np.array(f['grid_info/dz'])
        z0 = np.array(f['grid_info/z0'])
        z1 = np.array(f['grid_info/z1'])

        variables = {"r": r, "T0": T0, "rho0": rho0, "p0": p0, "g": g, "grad_s0": grad_s0}
        info = {"nz": nz, "nz_ghost": nz_ghost, "nz_full": nz_full, "dz": dz, "z0": z0, "z1": z1}
        
    return variables, info


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
        self.cmap = "RdBu"
        self.norm = TwoSlopeNorm
        self.set_plot_params()
        
    def set_plot_params(self):
        self.num_quivers = 10
        variables, info = read_fg(self.folder.format(0))
        #print(variables['vx'])
        z_shape, x_shape = variables['vx'].shape
        self.quiver_stride_x = int(x_shape//self.num_quivers)
        self.quiver_stride_z = int(z_shape//self.num_quivers)


        self.quiver_stride = 100
        self.font_size = 13
        self.title_size = 15
        
        max_values = [0] * 6
        min_values = [0] * 6
        for i in range(self.num_snaps):
            variables, info = read_fg(self.folder.format(i))
            T1, rho1, p1, s1, vx, vz = variables["T1"], variables["rho1"], variables["p1"], variables["s1"], variables["vx"], variables["vz"]
            
            max_values = np.maximum(max_values, [T1.max(), rho1.max(), p1.max(), s1.max(), vx.max(), vz.max()])
            min_values = np.minimum(min_values, [T1.min(), rho1.min(), p1.min(), s1.min(), vx.min(), vz.min()])
        
        #max_values = [max(abs(min_val), abs(max_val)) for min_val, max_val in zip(min_values, max_values)]
        #min_values = [-val for val in max_values]



        variables, info = read_fg(self.folder.format(self.num_snaps-1))
        T1, rho1, p1, s1, vx, vz = variables["T1"], variables["rho1"], variables["p1"], variables["s1"], variables["vx"], variables["vz"]

        t, nx, nz, nz_full, nz_ghost, dx, dz, z0, z1, x0, x1 = info['t'], info['nx'], info['nz'], info['nz_full'], info['nz_ghost'], info['dx'], info['dz'], info['z0'], info['z1'], info['x0'], info['x1']
        
        self.x_0 = x0/R_sun
        self.x_1 = x1/R_sun
        self.z_0 = z0/R_sun
        self.z_1 = z1/R_sun
        
        self.aspect = (self.x_1-self.x_0)/(self.z_1-self.z_0) * 1.2
        
        self.t_end = t*u.s
        if t>1e4:
            self.t_end = self.t_end.to("h")
        if t>1e5:
            self.t_end = self.t_end.to("day")
        
        # multiply by big number to decrease quiver scale, multiply by small number to increase quiver scale
        self.quiver_scale = 0.1*max(max(max_values[4], max_values[5]), max(np.abs(min_values[4]), np.abs(min_values[5])))
        
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
        variables, info = read_fg(self.folder.format(snap_nr))
        t, nx, nz, nz_ghost, dx, dz = info['t'], info['nx'], info['nz'], info['nz_ghost'], info['dx'], info['dz']   

        d = variables[key][nz_ghost:-1-nz_ghost+1,:]
        if self.plot_params[key][0]:
            vx = variables["vx"]
            vz = variables["vz"]
            vx = vx[nz_ghost:-1-nz_ghost+1,:]
            vz = vz[nz_ghost:-1-nz_ghost+1,:]

            Y, X = np.mgrid[0:d.shape[0]:self.quiver_stride_z, 0:d.shape[1]:self.quiver_stride_x]
            # Convert pixel indices to data coordinates
            X_data = (X * dx)/R_sun
            Y_data = (Y * dz)/R_sun+self.z_0
            self.quiver_scale = None
            ax.quiver(X_data, Y_data, vx[::self.quiver_stride_z, ::self.quiver_stride_x], vz[::self.quiver_stride_z, ::self.quiver_stride_x], scale=self.quiver_scale)
        
        t = t*u.s
        t = t.to(self.t_end.unit)
        
        ax.set_title(format_title(key))
        ax.set_xlabel("x [Solar radii]", fontsize=self.font_size)
        ax.set_ylabel("z [Solar radii]", fontsize=self.font_size)
        
        vmin = self.plot_params[key][1]
        vmax = self.plot_params[key][2]
        if self.norm == TwoSlopeNorm:
            norm = self.norm(vcenter=0, vmin=vmin, vmax=vmax)
            im =ax.imshow(d, origin="lower", extent=[self.x_0,self.x_1,self.z_0,self.z_1], aspect=self.aspect,norm=norm, cmap=self.cmap)
        else:
            im =ax.imshow(d, origin="lower", extent=[self.x_0,self.x_1,self.z_0,self.z_1], aspect=self.aspect,vmin=vmin, vmax=vmax, cmap=self.cmap)
        
        
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

        
    def animate_one(self, snap_nr, key, save=False, save_name=None, fps=4):
        pass
        
        def init_animation():
            im, t = self.plot(self.fig, self.ax, 0, key)
            cbar = plt.colorbar(im)
            # add colorbar

        def update_animation(snap_nr):
            """Update the plots for each frame."""
            for ax in self.fig.get_axes():
                ax.clear()
                #clear colorbar ect.

        self.fig, self.ax = plt.subplots(figsize=(6,6))

    def animate_all(self, save=False, save_name=None, fps=4, save_interval=10):
        def init_animation():
            self.plot_all(self.fig, 0)
            

        def update_animation(snap_nr):
            """Update the plots for each frame."""
            for ax in self.fig.get_axes():
                ax.clear()  # Clear previous data
            self.fig.clear()
            self.plot_all(self.fig, snap_nr)  # Plot the new snapshot

        self.fig = plt.figure(figsize=(16, 9))  # Create a figure
        #self.plot_all(self.fig,self.num_snaps-1)
        anim = FuncAnimation(self.fig, update_animation, interval=250, frames=range(0,self.num_snaps, save_interval), init_func=init_animation)
        if save:
            anim.save(save_name, writer='ffmpeg', fps=fps, extra_args=['-vcodec', 'libx264'])
        else:
            plt.show()

    def animate_all_no_vmin_vmax(self, save=False, save_name=None, fps=4, save_interval=10):
        for key in self.plot_params.keys():
            self.plot_params[key] = (self.plot_params[key][0], None, None, self.plot_params[key][3])
        self.animate_all(save, save_name, fps, save_interval)

    def plot_v_of_t(self):
        t = np.zeros(self.num_snaps)
        vx = np.zeros(self.num_snaps)
        vz = np.zeros(self.num_snaps)
        for i in range(self.num_snaps):
            variables, info = read_fg(self.folder.format(i))
            t[i] = info["t"]
            vx[i] = np.max(np.abs(variables["vx"]))
            vz[i] = np.max(np.abs(variables["vz"]))

        t = t*u.s
        t = t.to(self.t_end.unit)
        plt.plot(t, vx, label="vx")
        plt.plot(t, vz, label="vz")
        plt.xlabel("t [%s]" % self.t_end.unit.to_string(format='latex_inline'))
        #plt.xlabel("t [s]")
        plt.ylabel("Max Velocity [cm/s]")
        plt.legend()
        plt.show()

    def plot_p_of_t(self):
        t = np.zeros(self.num_snaps)
        p = np.zeros(self.num_snaps)
        for i in range(self.num_snaps):
            variables, info = read_fg(self.folder.format(i))
            t[i] = info["t"]
            p[i] = np.max(np.abs(variables["p1"]))
        plt.plot(t, p)
        plt.xlabel("t [s]")
        plt.ylabel("Max Pressure [dyn/cm$^2$]")
        plt.show()

    def plot_rho_of_t(self):
        t = np.zeros(self.num_snaps)
        rho = np.zeros(self.num_snaps)
        for i in range(self.num_snaps):
            variables, info = read_fg(self.folder.format(i))
            t[i] = info["t"]
            rho[i] = np.max(np.abs(variables["rho1"]))
        plt.plot(t, rho)
        plt.xlabel("t [s]")
        plt.ylabel(r"Max Density [g/cm$^3$]")
        plt.show()

    def plot_T_of_t(self):
        t = np.zeros(self.num_snaps)
        T = np.zeros(self.num_snaps)
        for i in range(self.num_snaps):
            variables, info = read_fg(self.folder.format(i))
            t[i] = info["t"]
            T[i] = np.max(np.abs(variables["T1"]))
        plt.plot(t, T)
        plt.xlabel("t [s]")
        plt.ylabel("Max Temperature [K]")
        plt.show()

    def plot_s_of_t(self):
        t = np.zeros(self.num_snaps)
        s = np.zeros(self.num_snaps)
        for i in range(self.num_snaps):
            variables, info = read_fg(self.folder.format(i))
            t[i] = info["t"]
            s[i] = np.max(np.abs(variables["s1"]))
            
        plt.plot(t, s)
        plt.xlabel("t [s]")
        plt.ylabel("Max Entropy [erg/K]")
        plt.show()

if __name__ == "__main__":
    directory = "../data/soft_wall_rk2_upw2_bigger/"
    save_name = directory.split("/")[-2] + ".mp4"

    vf = Visualize_Foreground(directory)
    #vf.plot_all(plt.figure(figsize=(16,9)), 1)
    #plt.show()
    if True:
        vf.plot_params['T1'] = (True, None, None, vf.plot_params['T1'][3])
        vf.plot_params['rho1'] = (True, None, None, vf.plot_params['rho1'][3])
        vf.plot_params['s1'] = (True, None, None, vf.plot_params['s1'][3])
        vf.plot_params['p1'] = (True, None, None, vf.plot_params['p1'][3])
        vf.plot_params['vx'] = (False, None, None, vf.plot_params['vx'][3])
        vf.plot_params['vz'] = (False, None, None, vf.plot_params['vz'][3])

        #vf.norm = NoNorm
        vf.animate_all(save=True, save_name=save_name, fps=6, save_interval=1)
    if False:
        vf.animate_all_no_vmin_vmax(save=True, save_name=save_name, fps=4, save_interval=1)