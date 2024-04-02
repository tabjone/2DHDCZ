#ifndef GLOBAL_INITIALIZATION_H__
#define GLOBAL_INITIALIZATION_H__

#define IC_SOD_SHOCK 0 // 0 for off, 1 for on
#define IC_SOD_SHOCK_DIRECTION 1 // 0 for horizontal, 1 for vertical

#define IC_RANDOM_OSCILLATIONS 0 // 0 for off, 1 for on

#define IC_ENTROPY_PERTURBATION 0 // 0 for off, 1 for on
#define IC_N_ENTROPY_PERTURBATION 3 // the number of gaussians for the entropy perturbation
// Placement of gaussian
// Percentage of box in x-y. Heigh in solar radii in z.
#define IC_ENTROPY_CENTRE_Z {0.80, 0.90, 0.85}
#define IC_ENTROPY_CENTRE_Y {0.40, 0.60, 0.45}
#define IC_ENTROPY_CENTRE_X {0.0, 0.0}
#define IC_ENTROPY_AMPLITUDE_DIRECTION {1.0,-1.0, -1.0}
// Sigma of gaussian in solar radii
#define IC_ENTROPY_SIGMA_Z 0.01
#define IC_ENTROPY_SIGMA_Y 0.01
#define IC_ENTROPY_SIGMA_X 0.01

#define IC_DENSITY_PERTURBATION 0 // 0 for off, 1 for on
#define IC_N_DENSITY_PERTURBATION 2 // the number of gaussians for the entropy perturbation
// Placement of gaussian between 0 and 1 in each direction for inside the grid
#define IC_DENSITY_CENTRE_Z {0.3, 0.6}
#define IC_DENSITY_CENTRE_Y {0.5, 0.6}
#define IC_DENSITY_CENTRE_X {0.5, 0.6}

#define IC_OSCILLATION_MODES 1 // 0 for off, 1 for on
#define IC_OSCILLATION_MODES_N_NUM 1 // the number of oscillation modes in the z-direction
#define IC_OSCILLATION_MODES_M_NUM 1 // the number of oscillation modes in the y-direction
#define IC_OSCILLATION_MODES_L_NUM 1 // the number of oscillation modes in the x-direction
#define IC_OSCILLATION_MODES_N {1} // the oscillation modes in the z-direction
#define IC_OSCILLATION_MODES_M {1} // the oscillation modes in the y-direction
#define IC_OSCILLATION_MODES_L {1} // the oscillation modes in the x-direction
#define IC_OSCILLATION_MODES_PHI_M {0.0} // the phase shift of the oscillation modes in the y-direction
#define IC_OSCILLATION_MODES_PHI_L {0.0} // the phase shift of the oscillation modes in the x-direction

#define IC_ZEROS 0 // 0 for off, 1 for on

#if (IC_SOD_SHOCK + IC_ENTROPY_PERTURBATION + IC_DENSITY_PERTURBATION + IC_OSCILLATION_MODES + IC_ZEROS) > 1
    #error "You can only have one initialization type"
#endif 

// Background
#define CONSTANT_BACKGROUND 0 // 0 for non-constant background, 1 for constant background
#define CONSTANT_BACKGROUND_DENSITY 1.0e-1 // Constant background density
#define CONSTANT_BACKGROUND_PRESSURE 1.0e13 // Constant background pressure
#define CONSTANT_BACKGROUND_TEMPERATURE 1.0e6 // Constant background temperature

#endif // GLOBAL_INITIALIZATION_H__