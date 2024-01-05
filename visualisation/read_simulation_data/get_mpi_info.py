import h5py

def get_mpi_info(file_path):
    with h5py.File(file_path, 'r') as f:
        num_processes = f['/total_processes'][()]
    return num_processes