#ifndef GLOBAL_PARAMETERS_H__
#define GLOBAL_PARAMETERS_H__

// Run name, order, time ect.
#define RUN_NAME "save_test"

#define T 1e5 // Simulation time in seconds
#define MAX_DT 1 // Maximum time step in seconds
#define SAVE_INTERVAL 10 // Save interval in seconds
#define SAVE_ALL 0 // Save all time steps

#define CFL_CUT 0.1

#define UPWIND_ORDER 1
#define CENTRAL_ORDER 2
#define TIME_ORDER 2

// Dimensions and MHD
#define DIMENSIONS 2
#define MHD 0

// Grid size
#define R_START 0.6 // In units of solar radii
#define R_END 0.97 // In units of solar radii
#define X_SIZE 0.05 // In units of solar radii
#define NX 10 // Number of grid points in x-direction
#define NZ 10 // Number of grid points in z-direction

// Gauss-Seidel tolerance, max iterations
#define GS_TOL 1e-5
#define GS_MAX_ITER 1e6

// Physical parameters
#define DEL_AD 0.4 // Adiabatic temperature gradient
#define GAMMA 5.0/3 // Superadiabatic parameter
#define NABLA_AD 0.4

// Background parameters
#define K 0.01 // superadiabacicity parameter in CZ
#define p_step 0.1 // Number for determining step-size in space

// Extra parameters
// 0 for constant, others not implemented yet
#define EXTRAPOLATE_GHOST_CELLS 0

#define DEBUG 1

#endif