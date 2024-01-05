import h5py
import numpy as np
from .get_mpi_info import get_mpi_info

def get_background_variable(folder, variable_key):
    # Get the number of processes
    num_processes = get_mpi_info(folder+'mpi_info.h5')

    # Variable from all processes
    variable = []

    for i in range(num_processes):
        with h5py.File(f'{folder}background_{i}.h5', 'r') as f:
            dataset = f[f'variables/{variable_key}']
            process_variable = np.array(dataset)
            unit_attr = dataset.attrs['unit']
            unit = unit_attr[0].decode('utf-8') if isinstance(unit_attr[0], bytes) else unit_attr[0]
            nz_ghost = np.array(f['grid_info/nz_ghost'])

            if num_processes ==1:
                variable.append(process_variable)
            else:
                if i == 0:
                    variable.append(process_variable[:-nz_ghost])
                elif i == num_processes - 1:
                    variable.append(process_variable[nz_ghost:])
                else:
                    variable.append(process_variable[nz_ghost:-nz_ghost])
    
    # Concatenating the data for each variable
    variable = np.concatenate(variable, axis=0)

    return variable, unit