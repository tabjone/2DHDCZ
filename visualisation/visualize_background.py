import matplotlib.pyplot as plt
import h5py
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
#import seaborn as sns
from visualize_foreground import *

R_sun = 6.957e10
G = 6.6743e-8

def read_solar_S_hdf5(file_path):
    with h5py.File(file_path, 'r') as f:
        r_over_R = np.array(f['/r_over_R'])
        c_s = np.array(f['/c_s'])
        Gamma_1 = np.array(f['/Gamma_1'])
        T0 = np.array(f['/T'])
        rho0 = np.array(f['/rho'])
        p0 = np.array(f['/p'])
        
        r_over_R = np.flip(r_over_R)
        c_s = np.flip(c_s)
        Gamma_1 = np.flip(Gamma_1)
        T0 = np.flip(T0)
        rho0 = np.flip(rho0)
        p0 = np.flip(p0)

        variables = {'r_over_R': r_over_R, 'c_s': c_s, 'Gamma_1': Gamma_1, 'T0': T0, 'rho0': rho0, 'p0': p0}
    return variables

solar_S = read_solar_S_hdf5("../additional_files/solar_s.h5")

def plot_background(SAVE_DIR, RUN_NAME, SAVE_NAME):
    n_procs = read_mpi_info(SAVE_DIR+RUN_NAME+"mpi_info.h5")
    var, info = read_bg_mpi(n_procs, SAVE_DIR + RUN_NAME)

    solar_S['M'] = 4*np.pi*R_sun**3 * cumtrapz(solar_S['rho0']*(solar_S['r_over_R'])**2, solar_S['r_over_R'], initial=0)
    solar_S['g'] = -G*solar_S['M']/((solar_S['r_over_R']*R_sun)**2)

    var['H'] = var['p0']/(var['rho0']*var['g'])
    solar_S['H'] = -solar_S['p0']/(solar_S['rho0']*solar_S['g'])

    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(16, 9))

    ax[0][0].plot(solar_S['r_over_R'], solar_S['T0'], color="black", linestyle=":", linewidth=2, label="Solar S")
    ax[0][0].plot(var['r']/R_sun, var['T0'], color="black", linestyle="--", linewidth=2, label="Our model")

    ax[1][0].semilogy(solar_S['r_over_R'], solar_S['p0'], color="black", linestyle=":", linewidth=2, label="Solar S")
    ax[1][0].semilogy(var['r']/R_sun, var['p0'], color="black", linestyle="--", linewidth=2, label="Our model")

    ax[0][1].semilogy(solar_S['r_over_R'], solar_S['rho0'], color="black", linestyle=":", linewidth=2, label="Solar S")
    ax[0][1].semilogy(var['r']/R_sun, var['rho0'], color="black", linestyle="--", linewidth=2, label="Our model")
    ax[0][1].set_xlim(0.6, 0.97)
    ax[0][1].set_ylim(1e-3, 1e0)
    #ax.yaxis.set_major_formatter(formatter)
    ax[0][1].set_ylabel("Density [g/cm$^3$]", fontsize=16)

    ax[1][1].plot(solar_S['r_over_R'], np.abs(solar_S['g']), color="black", linestyle=":", linewidth=2, label="Solar S")
    ax[1][1].plot(var['r']/R_sun, var['g'], color="black", linestyle="--", linewidth=2, label="Our model")
    ax[1][1].set_xlim(0.6, 0.97)
    ax[1][1].set_ylim(0.2e5, 0.8e5)
    ax[1][1].set_xlabel("Radius [Solar radii]", fontsize=16)
    ax[1][1].set_ylabel("Gravitational acceleration [cm s$^{-2}$]", fontsize=16)

    ax[0][2].semilogy(solar_S['r_over_R'], solar_S['H'], color="black", linestyle=":", linewidth=2, label="Solar S")
    ax[0][2].semilogy(var['r']/R_sun, var['H'], color="black", linestyle="--", linewidth=2, label="Our model")
    ax[0][2].set_ylabel("Pressure scale height [cm]", fontsize=16)

    ax[1][2].plot(var['r']/R_sun, var['grad_s0'], color="black", linestyle="--", linewidth=2, label="Our model")
    ax[1][2].set_ylabel("Entropy gradient [erg/(cm K)]", fontsize=16)

    for ax_ in ax:
        for axis in ax_:
            axis.tick_params(axis='both', labelsize=14)
            axis.legend(fontsize=14)
            axis.set_xlim(0.6, 1.0)
            #axis.yaxis.set_major_formatter(formatter)
            #axis.grid(True, which="both", ls="-", color='0.65')

    for i in range(3):
        ax[1][i].set_xlabel("Normalized radius R/R$_*$", fontsize=16)

    ax[0][0].set_ylabel("Temperature [K]", fontsize=16)
    ax[1][0].set_ylabel("Pressure [dyn/cm$^2$]", fontsize=16)
    fig.tight_layout()

    plt.savefig(SAVE_NAME)