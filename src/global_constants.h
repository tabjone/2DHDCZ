#ifndef GLOBAL_CONSTANTS_H__
#define GLOBAL_CONSTANTS_H__

#include "global_parameters.h"

#define ETA 0.1 // magnetic diffusivity [cm2 s-1]

#define MU 0.6 // mean molecular weight

#if UNITS == 0
    // CGS units
    #define K_B 1.3807e-16 // boltzmann constant [cm2 g s-2 K-1]
    #define M_U 1.66053907e-24 // atomic mass constant [g]
    #define R_SUN 6.957e10 // solar radius [cm]
    #define M_SUN 1.989e33 // solar mass [g]
    

    #if GRAVITY_ON == 0
        #define G 0.0 // gravitational constant [cm3 g-1 s-2]
    
    #else
        #define G 6.6743e-8 // gravitational constant [cm3 g-1 s-2]
    
    #endif // GRAVITY_ON == 0

#else
    // SI units
    #define K_B 1.3807e-23 // boltzmann constant [m2 kg s-2 K-1]
    #define M_U 1.66053907e-27 // atomic mass constant [kg]
    #define R_SUN 6.957e8 // solar radius [m]
    #define M_SUN 1.989e30 // solar mass [kg]

    #if GRAVITY_ON == 0
        #define G 0.0 // gravitational constant [m3 kg-1 s-2]

    #else
        #define G 6.6743e-11 // gravitational constant [m3 kg-1 s-2]
    
    #endif // GRAVITY_ON == 0
#endif // UNITS == 0
#endif // GLOBAL_CONSTANTS_H__