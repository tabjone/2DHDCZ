#ifndef GLOBAL_PARAMETERS_H__
#define GLOBAL_PARAMETERS_H__

#include "hdf5.h"

#define RUN_NAME "test_pert_rk1_upw1"
#define LOAD 0 // 0 for not loading, 1 for loading
#define LOAD_SNAP_NUMBER 0 // Snap number to load

#define T 1e4 // Simulation time in seconds
#define MAX_DT 1 // Maximum time step in seconds
#define SAVE_INTERVAL 10 // Save interval in seconds
#define SAVE_ALL 0 // 0 for saving on interval above, 1 for saving all time steps

#define CFL_CUT 0.96 // CFL cut

// Order of the scheme (only upwind 1,2, central 2,4, rk1,2 implemented so far)
#define UPWIND_ORDER 1 // 1 for first order, 2 for second order
#define CENTRAL_ORDER 2 // 2 for second order, 4 for fourth order
#define TIME_ORDER 1 // 1,2,3,4 for RK1, RK2, RK3, RK4

#define FLOAT_PRECISION 1 // 0 for float, 1 for double, 2 for long double

#if FLOAT_PRECISION == 0
    #define FLOAT_P float
    #define H5_FLOAT_P H5T_NATIVE_FLOAT
#elif FLOAT_PRECISION == 1
    #define FLOAT_P double
    #define H5_FLOAT_P H5T_NATIVE_DOUBLE
#elif FLOAT_PRECISION == 2
    #define FLOAT_P long double
    #define H5_FLOAT_P H5T_NATIVE_LDOUBLE
#else
    #error "Invalid FLOAT_PRECISION type specified."
#endif

// Dimensions and MHD (only 2D HD implemented so far)
#define DIMENSIONS 2 // 1 for 1D, 2 for 2D, 3 for 3D
#define MHD 0 // 0 for hydro, 1 for MHD

// Grid size
#define R_START 0.6 // In units of solar radii
#define R_END 0.97 // In units of solar radii
#define X_SIZE 0.05 // In units of solar radii
#define NX 60 // Number of grid points in x-direction
#define NZ 500 // Number of grid points in z-direction

// Gauss-Seidel tolerance, max iterations
#define GS_TOL 1e-6 // Gauss-Seidel tolerance
#define GS_MAX_ITER 2e6 // Gauss-Seidel max iterations

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