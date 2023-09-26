#ifndef GLOBAL_PARAMETERS_H__
#define GLOBAL_PARAMETERS_H__

#define RUN_NAME "vzvx_ten"
#define LOAD 0 // 0 for not loading, 1 for loading
#define LOAD_SNAP_NUMBER 0 // Snap number to load

#define FLOAT_P double // Choose precision of floating point numbers

#define T 1e5 // Simulation time in seconds
#define MAX_DT 1 // Maximum time step in seconds
#define SAVE_INTERVAL 200 // Save interval in seconds
#define SAVE_ALL 0 // 0 for saving on interval above, 1 for saving all time steps

#define CFL_CUT 0.98 // CFL cut

// Order of the scheme (only upwind 1,2, central 2,4, rk1,2 implemented so far)
#define UPWIND_ORDER 1 // 1 for first order, 2 for second order
#define CENTRAL_ORDER 2 // 2 for second order, 4 for fourth order
#define TIME_ORDER 2 // 1,2,3,4 for RK1, RK2, RK3, RK4

// Dimensions and MHD (only 2D HD implemented so far)
#define DIMENSIONS 2 // 1 for 1D, 2 for 2D, 3 for 3D
#define MHD 0 // 0 for hydro, 1 for MHD

// Grid size
#define R_START 0.6 // In units of solar radii
#define R_END 0.97 // In units of solar radii
#define X_SIZE 0.05 // In units of solar radii
#define NX 100 // Number of grid points in x-direction
#define NZ 1000 // Number of grid points in z-direction

// Gauss-Seidel tolerance, max iterations
#define GS_TOL 1e-5 // Gauss-Seidel tolerance
#define GS_MAX_ITER 1e6 // Gauss-Seidel max iterations

// Physical parameters
#define DEL_AD 0.4 // Adiabatic temperature gradient
#define GAMMA 5.0/3 // Superadiabatic parameter
#define NABLA_AD 0.4 // Adiabatic temperature gradient

// Background parameters
#define K 0.01 // superadiabacicity parameter in CZ
#define p_step 0.1 // Number for determining step-size in space

// Extra parameters
// Extrapolation of ghost cells
// 0 for constant, others not implemented yet
#define EXTRAPOLATE_GHOST_CELLS 0

#define DEBUG 1

#endif