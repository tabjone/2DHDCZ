import subprocess
import time
import numpy as np

def update_ny_nz(file_path, new_n:int):
    # Read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Modify the lines containing NY and NZ
    for i, line in enumerate(lines):
        if line.strip().startswith('#define NY'):
            lines[i] = f'#define NY {new_n} // Number of grid points in y-direction\n'
        elif line.strip().startswith('#define NZ'):
            lines[i] = f'#define NZ {new_n} // Number of grid points in z-direction\n'

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)


def compile_program():
    """
    Compiles the MPI program by calling the Makefile's `compile` target.
    """
    # Call the make command to compile the MPI program
    subprocess.run("make clean", shell=True, check=True)
    subprocess.run("make compile", shell=True, check=True)

def run_program(num_processes):
    """
    Runs the MPI program with the specified number of processes by setting the NP environment variable
    and then calling the Makefile's `execute` target.
    :param num_processes: Number of processes to run the MPI program with.
    :param make_command: The command to compile and execute the MPI program.
    """
    make_command = f"mpiexec -n {num_processes} ./main_mpi.o"

    # Call the make command to execute the MPI program
    start_time = time.time()
    subprocess.run(make_command, shell=True, check=True)
    end_time = time.time()

    return end_time - start_time

def calculate_grid_points(num_processes):
    return int(np.sqrt(num_processes*50*50))

def main():
    # calculate grid points
    # update grid points
    #compile
    #execute
    #repeat 10 times per process

    # Number of processes to test
    processes = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    
    # Number of runs per number of processes
    runs_per_process = 10

    times = []
    for num in processes:
        times_process = []
        new_n = calculate_grid_points(num)
        update_ny_nz("src/simulation_parameters/global_parameters.h", new_n)
        compile_program()
        for run in range(runs_per_process):
            print(f"Running with {num} processes, iteration {run+1}...")
            exec_time = run_program(num)
            times_process.append(exec_time)
            print(f"Execution time: {exec_time:.4f} seconds")


        times.append(times_process)
    
    # Save results to a file
    with open('mpi_runtimes_weak.txt', 'w+') as file:
        file.write("Processes, Runtimes (s)\n")
        for i, num in enumerate(processes):
            file.write(f"{num}, {', '.join([str(time) for time in times[i]])}\n")

if __name__ == "__main__":
    main()
    #compile_program()