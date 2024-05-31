#ifndef GLOBAL_CONSTANTS_H__
#define GLOBAL_CONSTANTS_H__

#include "math.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#define ETA 0.1 // magnetic diffusivity [cm2 s-1]

#define MU 0.594 // mean molecular weight

#define OMEGA_EQ 2.87e-6 // equatorial angular velocity [rad/s]
#define OMEGA_CORE 2.9e-6 // core angular velocity [rad/s]

#if UNITS == 0
    // CGS units
    #define K_B 1.380649e-16 // boltzmann constant [erg/K]
    #define M_U 1.66053907e-24 // atomic mass constant [g]
    #define R_SUN 6.957e10 // solar radius [cm]
    #define M_SUN 1.989e33 // solar mass [g]

    #define VISCOSITY_COEFF (1.0e7/5.37) // viscosity coefficient [P] (Poise) (originally 1.0e7/5.37)
    #define THERMAL_DIFFUSIVITY_COEFF 1.866e6 // thermal diffusivity [...] (orginally 1.866e6)

    #if GRAVITY_ON == 0
        #define G 0.0 // gravitational constant [cm3 g-1 s-2]
    
    #else
        #define G 6.6743e-8 // gravitational constant [cm3 g-1 s-2]
    
    #endif // GRAVITY_ON == 0

#else
    // SI units
    #define K_B 1.380649e-23 // boltzmann constant [J/K]
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