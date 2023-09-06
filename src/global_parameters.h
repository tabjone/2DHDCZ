#ifndef GLOBAL_PARAMETERS_H__
#define GLOBAL_PARAMETERS_H__

#define RUN_NAME "HD_2D_1"

#define DIMENSIONS 2
#define MHD 0

#define UPWIND_ORDER 2
#define CENTRAL_ORDER 2

#define R_START 0.6 // In units of solar radii
#define R_END 0.97 // In units of solar radii
#define DEL_AD 0.4 // Adiabatic temperature gradient
#define GAMMA 5.0/3 // Superadiabatic parameter

#define K 0.01 // superadiabacicity parameter in CZ
#define p 0.1 // Number for determining step-size in space
#define NABLA_AD 0.4

#define x_size 0.01 // In units of solar radii
#define dx 0.0001 // Step-size in x-direction in units of solar radii
#define dy 0.0001 // Step-size in y-direction in units of solar radii
#define dz 0.0001 // Step-size in z-direction in units of solar radii

#endif