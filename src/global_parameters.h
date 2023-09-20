#ifndef GLOBAL_PARAMETERS_H__
#define GLOBAL_PARAMETERS_H__

#define RUN_NAME "HD_2D_1"

#define DIMENSIONS 2
#define MHD 0

#define UPWIND_ORDER 2
#define CENTRAL_ORDER 2
#define TIME_ORDER 1

#define R_START 0.6 // In units of solar radii
#define R_END 0.97 // In units of solar radii
#define DEL_AD 0.4 // Adiabatic temperature gradient
#define GAMMA 5.0/3 // Superadiabatic parameter

#define K 0.01 // superadiabacicity parameter in CZ
#define p_step 0.1 // Number for determining step-size in space
#define NABLA_AD 0.4
#define GAMMA 5.0/3 // Superadiabatic parameter

#define X_SIZE 0.05 // In units of solar radii

// This should be split for x and z directions
#define RESOLUTION 100 // ...


// 0 for constant, others not implemented yet
#define EXTRAPOLATE_GHOST_CELLS 0

#define T 500 // Simulation time in seconds
#define SAVE_INTERVAL 1.0 // Save interval in seconds

#define DEBUG 1

#endif