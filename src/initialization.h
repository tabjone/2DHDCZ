#ifndef INITIALIZATION_H__
#define INITIALIZATION_H__

#define IC_SOD_SHOCK 0 // 0 for off, 1 for on
#define IC_SOD_SHOCK_DIRECTION 0 // 0 for horizontal, 1 for vertical

#define IC_ENTROPY_PERTUBATION 1 // 0 for off, 1 for on
#define IC_N_ENTROPY_PERTUBATION 2 // the number of gaussians for the entropy pertubation
// Placement of gaussian between 0 and 1 in each direction for inside the grid
#define IC_ENTROPY_CENTRE_Z {0.2, 0.6}
#define IC_ENTROPY_CENTRE_Y {0.2, 0.6}

#define IC_DENSITY_PERTUBATION 0 // 0 for off, 1 for on
#define IC_N_DENSITY_PERTUBATION 2 // the number of gaussians for the entropy pertubation
// Placement of gaussian between 0 and 1 in each direction for inside the grid
#define IC_DENSITY_CENTRE_Z {0.2, 0.6}
#define IC_DENSITY_CENTRE_Y {0.2, 0.6}
#define IC_ENTROPY_CENTRE_X {0.2, 0.6}


#if (IC_SOD_SHOCK + IC_ENTROPY_PERTUBATION + IC_DENSITY_PERTUBATION) > 1
    #error "You can only have one initialization type"
#endif

#endif // INITIALIZATION_H__