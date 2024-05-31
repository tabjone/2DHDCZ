import subprocess
import time
import numpy as np
import os

def run_mpi_program(num_processes):
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

def main():
    # Number of processes to test
    processes = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    
    # Number of runs per number of processes
    runs_per_process = 10
    
    times = []
    for num in processes:
        times_process = []
        for run in range(runs_per_process):
            print(f"Running with {num} processes, iteration {run+1}...")
            exec_time = run_mpi_program(num)
            times_process.append(exec_time)
            print(f"Execution time: {exec_time:.4f} seconds")
        
        times.append(times_process)
    
    # Save results to a file
    with open('mpi_runtimes_2.txt', 'w+') as file:
        file.write("Processes, Runtimes (s)\n")
        for i, num in enumerate(processes):
            file.write(f"{num}, {', '.join([str(time) for time in times[i]])}\n")

    print("\nResults saved to mpi_runtimes_2.txt")

if __name__ == "__main__":
    main()