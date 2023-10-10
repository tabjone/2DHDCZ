#ifndef GLOBAL_PARAMETERS_H__
#define GLOBAL_PARAMETERS_H__

#include "hdf5.h"
#include <mpi.h>

#define INITIALIZATION_TYPE 2 // 0 for zeros, 1 for velocity right, 2 for density pertubation, 3 for random

#define MPI_ON 0 // 0 for MPI off, 1 for MPI on

#define RUN_NAME "rk1_test_test"
#define LOAD 0 // 0 for not loading, 1 for loading
#define LOAD_SNAP_NUMBER 0 // Snap number to load

#define T 1e5 // Simulation time in seconds
#define MAX_DT 1e2 // Maximum time step in seconds
#define SAVE_INTERVAL 1e2 // Save interval in seconds
#define SAVE_ALL 0 // 0 for saving on interval above, 1 for saving all time steps

#define CFL_CUT 0.8 // CFL cut

// Order of the scheme (only upwind 1,2, central 2,4, rk1,2 implemented so far)
#define UPWIND_ORDER 1 // 1 for first order, 2 for second order
#define CENTRAL_ORDER 2 // 2 for second order, 4 for fourth order
#define TIME_ORDER 1 // 1,2,3 for RK1, RK2, RK3

#define FLOAT_PRECISION 1 // 0 for float, 1 for double, 2 for long double

// Warning, the units might not be correct in the equations of the system for SI
#define UNITS 0 // 0 for cgs, 1 for SI

// Dimensions and MHD (only 2D HD implemented so far)
#define DIMENSIONS 2 // 1 for 1D, 2 for 2D, 3 for 3D
// THIS WILL BE CHANGED TO MHD INTEGRATED IN FUNCTIONS BUT BFIELD_ON WILL CHOSE TO INCLUDE OR NOT INCLUDE B-FIELD
#define MHD 0 // 0 for hydro, 1 for MHD

#define VERTICAL_BOUNDARY_TYPE 1 // 0 for Hard-wall, 1 for spunge

// Grid size
#define CZ_START 0.7 // In units of solar radii
#define R_START 0.7 // In units of solar radii
#define R_END 0.90 // In units of solar radii
#define X_SIZE 0.05 // In units of solar radii
#define Y_SIZE 0.05 // In units of solar radii
#define NX 40 // Number of grid points in x-direction
#define NY 40 // Number of grid points in y-direction
#define NZ 160 // Number of grid points in z-direction

// Gauss-Seidel tolerance, max iterations
#define GS_TOL 1e-6 // Gauss-Seidel tolerance
#define GS_MAX_ITER 5e6 // Gauss-Seidel max iterations

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

// Debugging

#define GRAVITY_ON 1 // 0 for gravity off, 1 for gravity on
#define ADVECTION_ON 1 // 0 for advection off, 1 for advection on
#define GAS_PRESSURE_ON 1 // 0 for gas pressure off, 1 for gas pressure on
#define THERMAL_DIFFUSIVITY_ON 1 // 0 for thermal diffusivity off, 1 for thermal diffusivity on
#define BFIELD_ON 0 // 0 for B-field off, 1 for B-field on

#define DEBUG 1




// Do not change anything below this line
// --------------------------------------

#if FLOAT_PRECISION == 0
    #define FLOAT_P float
    #define H5_FLOAT_P H5T_NATIVE_FLOAT
    #define MPI_FLOAT_P MPI_FLOAT
#elif FLOAT_PRECISION == 1
    #define FLOAT_P double
    #define H5_FLOAT_P H5T_NATIVE_DOUBLE
    #define MPI_FLOAT_P MPI_DOUBLE
#elif FLOAT_PRECISION == 2
    #define FLOAT_P long double
    #define H5_FLOAT_P H5T_NATIVE_LDOUBLE
    #define MPI_FLOAT_P MPI_LONG_DOUBLE
#else
    #error "Invalid FLOAT_PRECISION type specified."
#endif // FLOAT_PRECISION

#if UPWIND_ORDER > 2
    #error "Upwind order > 2 not implemented."
#endif // UPWIND_ORDER > 2
#if CENTRAL_ORDER > 4
    #error "Central order > 4 not implemented."
#endif // CENTRAL_ORDER > 4
#if DIMENSIONS != 2
    #error "Only 2D implemented."
#endif // DIMENSIONS != 2
#if MHD != 0
    #error "Only hydrodynamics implemented."
#endif // MHD != 0
#if BFIELD_ON != 0
    #error "Only hydrodynamics implemented."
#endif // BFIELD_ON != 0
#if EXTRAPOLATE_GHOST_CELLS != 0
    #error "Only constant extrapolation implemented."
#endif // EXTRAPOLATE_GHOST_CELLS != 0
#endif // GLOBAL_PARAMETERS_H__