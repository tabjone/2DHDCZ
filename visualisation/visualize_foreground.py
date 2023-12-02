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
from scipy.ndimage import zoom

from astropy import units as u

# Format colorbar ticks in scientific notation
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)

R_sun = 6.957e10

def resize_array(arr, target_size=300):
    if arr.shape[0] > target_size:
        # Calculate the zoom factor
        zoom_factor = target_size / arr.shape[0]
        # Apply zoom only on the first axis
        resized_arr = zoom(arr, zoom_factor)
        return resized_arr
    else:
        return arr

def read_mpi_info(file_path):
    with h5py.File(file_path, 'r') as f:
        n_procs = f['/total_processes'][()]
    return n_procs

def read_fg_mpi(snap, n_procs, folder):
    variables_list = {key: [] for key in ['T1', 'rho1', 'p1', 's1', 'vz', 'vy']}
    infos = []

    for i in range(n_procs):
        variable, info = read_fg("{}snap{}_{}.h5".format(folder, snap, i))
        infos.append(info)

        for key in variables_list.keys():
            if i == 0:
                variables_list[key].append(variable[key][:-info['nz_ghost']])
            elif i == n_procs - 1:
                variables_list[key].append(variable[key][info['nz_ghost']:])
            else:
                variables_list[key].append(variable[key][info['nz_ghost']:-info['nz_ghost']])

    # Concatenating the data for each variable
    total_variables = {key: np.concatenate(variables_list[key], axis=0) for key in variables_list.keys()}
    
    # Aggregating the information
    info = infos[0]
    info['z0'] = infos[0]['z0']
    info['z1'] = infos[-1]['z1']

    return total_variables, info

def read_bg_mpi(n_procs, folder):
    variables_list = {key: [] for key in ['r', 'T0', 'rho0', 'p0', 'g', 'grad_s0']}
    infos = []

    for i in range(n_procs):
        variable, info = read_bg("{}background_{}.h5".format(folder, i))
        infos.append(info)

        for key in variables_list.keys():
            if i == 0:
                variables_list[key].append(variable[key][:-info['nz_ghost']])
            elif i == n_procs - 1:
                variables_list[key].append(variable[key][info['nz_ghost']:])
            else:
                variables_list[key].append(variable[key][info['nz_ghost']:-info['nz_ghost']])

    # Concatenating the data for each variable
    total_variables = {key: np.concatenate(variables_list[key], axis=0) for key in variables_list.keys()}
    
    # Aggregating the information
    info = infos[0]
    info['z0'] = infos[0]['z0']
    info['z1'] = infos[-1]['z1']

    return total_variables, info

def read_fg(file_path):
    with h5py.File(file_path, 'r') as f:
        T1 = np.array(f['variables/T1'])
        rho1 = np.array(f['variables/rho1'])
        p1 = np.array(f['variables/p1'])
        s1 = np.array(f['variables/s1'])
        vy = np.array(f['variables/vy'])
        vz = np.array(f['variables/vz'])
        
        t = np.array(f['grid_info/t'])
        ny = np.array(f['grid_info/ny'])
        nz = np.array(f['grid_info/nz'])
        nz_full = np.array(f['grid_info/nz_full'])
        nz_ghost = np.array(f['grid_info/nz_ghost'])
        dy = np.array(f['grid_info/dy'])
        dz = np.array(f['grid_info/dz'])
        z0 = np.array(f['grid_info/z0'])
        z1 = np.array(f['grid_info/z1'])
        y0 = np.array(f['grid_info/y0'])
        y1 = np.array(f['grid_info/y1'])
        
        variables = {"T1": T1, "rho1": rho1, "p1": p1, "s1": s1, "vy": vy, "vz": vz}
        info = {"t": t, "ny": ny, "nz": nz, "nz_full": nz_full, "nz_ghost": nz_ghost, "dy": dy, "dz": dz, "z0": z0, "z1": z1, "y0": y0, "y1": y1}
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
    def __init__(self, folder, delta=False):
        self.folder = folder
        self.delta = delta

        snap_id = 0
        snaps = []

        self.n_procs = read_mpi_info(folder+"mpi_info.h5")
        
        if self.n_procs > 1:
            exit_flag = False  # Flag to indicate when to exit the while loop
            while True:
                for i in range(self.n_procs):
                    snap_filepath = folder + f"snap{snap_id}_{i}.h5"
                    if os.path.exists(snap_filepath):
                        snaps.append(snap_filepath)
                    else:
                        exit_flag = True  # Set the flag to True when the condition is met
                        break  # Break out of the for loop
                snap_id += 1
                if exit_flag:
                    break  # Break out of the while loop if the flag is set
        else:
            while True:
                snap_filepath = folder + f"snap{snap_id}.h5"
                if os.path.exists(snap_filepath):
                    snaps.append(snap_filepath)
                    snap_id += 1
                else:
                    break        
        self.snaps = snaps
        self.num_snaps = int(len(self.snaps)/self.n_procs)
        self.cmap = "RdBu"
        self.norm = TwoSlopeNorm
        self.set_plot_params()
        
    def set_plot_params(self):
        self.num_quivers = 10
        if self.n_procs > 1:
            variables, info = read_fg_mpi(self.num_snaps-1, self.n_procs, self.folder)
        else:
            variables, info = read_fg(self.folder+f"snap{self.num_snaps-1}.h5")
        
        z_shape, y_shape = variables['vy'].shape
        self.quiver_stride_y = int(y_shape//self.num_quivers)
        self.quiver_stride_z = int(z_shape//self.num_quivers)

        self.font_size = 13
        self.title_size = 15

        t, y0, y1, z0, z1 = info['t'], info['y0'], info['y1'], info['z0'], info['z1']

        self.y0 = y0/R_sun
        self.y1 = y1/R_sun
        self.z0 = z0/R_sun
        self.z1 = z1/R_sun
        
        self.aspect = (self.y1-self.y0)/(self.z1-self.z0) * 1.2
        
        self.t_end = t*u.s
        if t>1e3:
            self.t_end = self.t_end.to("h")
        if t>1e4:
            self.t_end = self.t_end.to("day")
        
        # holding touple(var_name, quiver true/false, vmin, vmax)
        if self.delta:
            self.plot_params = {
            "T1": (False, None, None, "Temperature [K]"),
            "rho1": (False, None, None, r"Density [g/cm$^3$]"),
            "p1": (False, None, None, r"Pressure [dyn/cm$^2$]"),
            "s1": (False, None, None, "Entropy [erg/K]"),
            "vy": (False, None, None, "y-velocity [cm/s]"),
            "vz": (False, None, None, "z-velocity [cm/s]")
        }  
        else:
            self.plot_params = {
                "T1": (True, None, None, "Temperature [K]"),
                "rho1": (True, None, None, r"Density [g/cm$^3$]"),
                "p1": (True, None, None, r"Pressure [dyn/cm$^2$]"),
                "s1": (True, None, None, "Entropy [erg/K]"),
                "vy": (False, None, None, "y-velocity [cm/s]"),
                "vz": (False, None, None, "z-velocity [cm/s]")
            }     

    def plot(self, fig, ax, snap_nr, key):
        if self.n_procs > 1:
            if self.delta:
                variables_new, info = variables, info = read_fg_mpi(snap_nr, self.n_procs, self.folder)
                variables_old, info = variables, info = read_fg_mpi(snap_nr-1, self.n_procs, self.folder)
                variables = {key: variables_new[key]-variables_old[key] for key in variables_new.keys()}
                background, _ = read_bg_mpi(self.n_procs, self.folder)
            else:
                variables, info = read_fg_mpi(snap_nr, self.n_procs, self.folder)
                background, _ = read_bg_mpi(self.n_procs, self.folder)
        else:
            if self.delta:
                variables_new, info = read_fg(self.folder+f"snap{snap_nr}.h5")
                variables_old, info = read_fg(self.folder+f"snap{snap_nr-1}.h5")
                variables = {key: variables_new[key]-variables_old[key] for key in variables_new.keys()}
                background, _ = read_bg(self.folder+"background.h5")
            else:
                variables, info = read_fg(self.folder+f"snap{snap_nr}.h5")
                background, _ = read_bg(self.folder+"background.h5")
        t, nz_ghost, dy, dz = info['t'], info['nz_ghost'], info['dy'], info['dz']  

        d = variables[key][nz_ghost:-1-nz_ghost+1,:]

        if self.delta:
            ax.set_title(r"$\Delta$"+format_title(key))
                
        else:
            if (key=="T1"):
                d = d/(background["T0"][nz_ghost:-1-nz_ghost+1])[:, np.newaxis]
                ax.set_title(format_title(key)+"$/T_0$")
            elif (key=="rho1"):
                d = d/(background["rho0"][nz_ghost:-1-nz_ghost+1])[:, np.newaxis]
                ax.set_title(format_title(key)+r"$/\rho_0$")
            elif (key=="p1"):
                d = d/(background["p0"][nz_ghost:-1-nz_ghost+1])[:, np.newaxis]
                ax.set_title(format_title(key)+"$/p_0$")
            else:
                ax.set_title(format_title(key))

        if self.plot_params[key][0]:
            vy = variables["vy"][nz_ghost:-1-nz_ghost+1,:]
            vz = variables["vz"][nz_ghost:-1-nz_ghost+1,:]

            Y, X = np.mgrid[0:d.shape[0]:self.quiver_stride_z, 0:d.shape[1]:self.quiver_stride_y]
            # Convert pixel indices to data coordinates
            X_data = (X * dy)/R_sun
            Y_data = (Y * dz)/R_sun+self.z0
            self.quiver_scale = None
            ax.quiver(X_data, Y_data, vy[::self.quiver_stride_z, ::self.quiver_stride_y], vz[::self.quiver_stride_z, ::self.quiver_stride_y], scale=self.quiver_scale)
        
        t = t*u.s
        t = t.to(self.t_end.unit)
        
        ax.set_xlabel("x [Solar radii]", fontsize=self.font_size)
        ax.set_ylabel("z [Solar radii]", fontsize=self.font_size)

        if self.plot_params[key][2] == None:
            vmax = np.max(np.abs(d))
            vmin = -vmax        

        extent = [self.y0,self.y1,self.z0,self.z1]


        if vmin < 0 and vmax > 0:
            if self.norm == TwoSlopeNorm:
                norm = self.norm(vcenter=0, vmin=vmin, vmax=vmax)
                im =ax.imshow(d, origin="lower", extent=extent, aspect=self.aspect,norm=norm, cmap=self.cmap)
        else:
            im =ax.imshow(d, origin="lower", extent=extent, aspect=self.aspect,vmin=vmin, vmax=vmax, cmap=self.cmap)
        return im, t

    def plot_all(self, fig, snap_nr):
        gs = gridspec.GridSpec(2, 6, width_ratios=[1, 0.05, 1, 0.05, 1, 0.05], wspace=0.6, hspace=0.3)

        # Order of keys for plotting
        plot_order = ["T1", "rho1", "vy", "p1", "s1", "vz"]

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
            #cbar.ax.yaxis.set_major_formatter(formatter)

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
            if self.delta:
                self.plot_all(self.fig, 1)
            else:
                self.plot_all(self.fig, 0)

        def update_animation(snap_nr):
            """Update the plots for each frame."""
            for ax in self.fig.get_axes():
                ax.clear()  # Clear previous data
            self.fig.clear()
            self.plot_all(self.fig, snap_nr)  # Plot the new snapshot

        self.fig = plt.figure(figsize=(16, 9))  # Create a figure
        
        if self.delta:
            anim = FuncAnimation(self.fig, update_animation, interval=250, frames=range(1,self.num_snaps, save_interval), init_func=init_animation)
        else:
            anim = FuncAnimation(self.fig, update_animation, interval=250, frames=range(0,self.num_snaps, save_interval), init_func=init_animation)
        if save:
            anim.save(save_name, writer='ffmpeg', fps=fps, extra_args=['-vcodec', 'libx264'])
        else:
            plt.show()

    def plot_rho_of_t(self):
        t = np.zeros(self.num_snaps)
        rho = np.zeros(self.num_snaps)
        for i in range(self.num_snaps):
            if self.n_procs > 1:
                variables, info = read_fg_mpi(i, self.n_procs, self.folder)
            else:
                variables, info = read_fg(self.folder + "{i}.h5")
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
            if self.n_procs > 1:
                variables, info = read_fg_mpi(i, self.n_procs, self.folder)
            else:
                variables, info = read_fg(self.folder + "{i}.h5")
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
            if self.n_procs > 1:
                variables, info = read_fg_mpi(i, self.n_procs, self.folder)
            else:
                variables, info = read_fg(self.folder + "{i}.h5")
            t[i] = info["t"]
            s[i] = np.max(np.abs(variables["s1"]))
            
        plt.plot(t, s)
        plt.xlabel("t [s]")
        plt.ylabel("Max Entropy [erg/K]")
        plt.show()

if __name__ == "__main__":
    DATA_FOLDER = "/mn/stornext/d10/data/tabjone/data/"
    RUN_NAME = "bigger_box_high_res_little_higher_k/"

    vf = Visualize_Foreground(DATA_FOLDER+RUN_NAME, delta=False)
    vf.norm = TwoSlopeNorm
    fig, ax = plt.subplots(figsize=(6,6))

    vf.plot(fig, ax, 1, "p1")
    plt.savefig("vy.png")