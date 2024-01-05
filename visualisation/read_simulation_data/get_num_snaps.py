from .get_mpi_info import get_mpi_info
import os

def get_num_snaps(folder):
    """
    Get the number of snapshots in the folder.
    """
    snaps = []
    snap_id = 0
    num_processes = get_mpi_info(folder+"mpi_info.h5")
    num_snaps = 0
    
    exit_flag = False  # Flag to indicate when to exit the while loop
    while True:
        for i in range(num_processes):
            snap_filepath = f'{folder}snap{snap_id}_{i}.h5'
            if os.path.exists(snap_filepath):
                snaps.append(snap_filepath)
            else:
                exit_flag = True  # Set the flag to True when the condition is met
                break  # Break out of the for loop
        snap_id += 1
        if exit_flag:
            break  # Break out of the while loop if the flag is set

    num_snaps = int(len(snaps)/num_processes)
    return num_snaps-1