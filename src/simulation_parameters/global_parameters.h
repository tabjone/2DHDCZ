#ifndef GLOBAL_PARAMETERS_H__
#define GLOBAL_PARAMETERS_H__

#define MPI_ON 1 // 0 for MPI off, 1 for MPI on

#define RUN_NAME "test_hello_4" // Name of the run
#define SAVE_DIR "/mn/stornext/d10/data/tabjone/data/"
#define LOAD 0 // 0 for not loading, 1 for loading
#define LOAD_SNAP_NUMBER 0 // Snap number to load

#define BACKGROUND_INTEGRATION_STARTING_POINT 0.985 // Starting point for integration of background
#define BACKGROUND_INTEGRATION_TOP 1.10 // Top point for integration of background
#define BACKGROUND_INTEGRATION_BOTTOM 0.65 // Bottom point for integration of background

#define T 1e6 // Simulation time in seconds
#define MAX_DT 1e3 // Maximum time step in seconds
#define SAVE_INTERVAL 1e1 // Save interval in seconds
#define SAVE_ALL 0 // 0 for saving on interval above, 1 for saving all time steps

#define CFL_CUT 0.1 // CFL cut

// Order of the scheme (only upwind 1,2, central 2,4, rk1,2 implemented so far)
#define UPWIND_ORDER 2 // 1 for first order, 2 for second order
#define CENTRAL_ORDER 2 // 2 for second order, 4 for fourth order
#define TIME_ORDER 3 // 1,2,3 for RK1, RK2, RK3

#define FLOAT_PRECISION 1 // 0 for float, 1 for double, 2 for long double

// Warning, the units might not be correct in the equations of the system for SI
#define UNITS 0 // 0 for cgs, 1 for SI

// Dimensions
#define DIMENSIONS 2 // 1 for 1D, 2 for 2D, 3 for 3D

// Grid size
#define CZ_START 0.71 // In units of solar radii
// R_START, R_END must be between 0.3 and 1.5
#define R_START 0.78 // In units of solar radii
#define R_END 0.98 // In units of solar radii
#define X_SIZE 0.1 // In units of solar radii
#define Y_SIZE 0.1 // In units of solar radii
#define NX 400 // Number of grid points in x-direction
#define NY 150 // Number of grid points in y-direction
#define NZ 150 // Number of grid points in z-direction

// Iterative solver parameters
#define ITERATIVE_SOLVER_TYPE 1 // 0 for Jacobi, 1 for Gauss-Seidel, 2 for SOR (not implemented)
#define ITERATIVE_SOLVER_TOLERANCE 1e-5 // Iterative solver tolerance
#define ITERATIVE_SOLVER_MAX_ITERATIONS (int)5e8 // Iterative solver max iterations

// Physical parameters
#define DEL_AD 0.4 // Adiabatic temperature gradient
#define GAMMA 5.0/3 // Superadiabatic parameter
#define NABLA_AD 0.4 // Adiabatic temperature gradient

// Background parameters
#define SUPERAD_PARAM 0.01 // superadiabacicity parameter in CZ
#define p_step 0.001 // Number for determining step-size in space


#define GRAVITY_ON 1 // 0 for gravity off, 1 for gravity on
#define ADVECTION_ON 1 // 0 for advection off, 1 for advection on
#define GAS_PRESSURE_ON 1 // 0 for gas pressure off, 1 for gas pressure on
#define VISCOSITY_ON 1 // 0 for viscosity off, 1 for viscosity on
#define THERMAL_DIFFUSIVITY_ON 1 // 0 for thermal diffusivity off, 1 for thermal diffusivity on
#define BFIELD_ON 0 // 0 for B-field off, 1 for B-field on

#define DEBUG 1
#define SAVE_RHS 1

// Do not change anything below this line
// --------------------------------------

#if UPWIND_ORDER > 2
    #error "Upwind order > 2 not implemented."
#endif // UPWIND_ORDER > 2
#if CENTRAL_ORDER > 4
    #error "Central order > 4 not implemented."
#endif // CENTRAL_ORDER > 4
#if BFIELD_ON != 0
    #error "Only hydrodynamics implemented."
#endif // BFIELD_ON != 0
#if EXTRAPOLATE_GHOST_CELLS != 0
    #error "Only constant extrapolation implemented."
#endif // EXTRAPOLATE_GHOST_CELLS != 0
#endif // GLOBAL_PARAMETERS_H__
