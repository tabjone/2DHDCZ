import numpy as np
import h5py

# Load the data skipping the first 5 lines
data = np.genfromtxt('cptrho.l5bi.d.15c', skip_header=5)

# Separate data into different arrays
r_over_R = data[:, 0]
c_s = data[:, 1]
rho = data[:, 2]
p = data[:, 3]
Gamma_1 = data[:, 4]
T = data[:, 5]

# Save the data into an HDF5 file
with h5py.File('solar_s.h5', 'w') as f:
    f.create_dataset('r_over_R', data=r_over_R)
    f.create_dataset('c_s', data=c_s)
    f.create_dataset('rho', data=rho)
    f.create_dataset('p', data=p)
    f.create_dataset('Gamma_1', data=Gamma_1)
    f.create_dataset('T', data=T)

    # Add metadata
    f.attrs['Description'] = 'Sound speed, etc. for Model S (Christensen-Dalsgaard et al. 1996; Science 272, 1286)'
    f.attrs['Download date, location'] = 'Tuesday, 11 July 2023 at 15:04 https://users-phys.au.dk/~jcd/solar_models/cptrho.l5bi.d.15c'
    f['r_over_R'].attrs['units'] = 'R/R'
    f['c_s'].attrs['units'] = 'cm/sec'
    f['rho'].attrs['units'] = 'g/cm^3'
    f['p'].attrs['units'] = 'dyn/cm^2'
    f['Gamma_1'].attrs['units'] = 'Unitless'
    f['T'].attrs['units'] = 'K'