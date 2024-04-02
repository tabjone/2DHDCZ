#ifndef GLOBAL_BOUNDARY_H__
#define GLOBAL_BOUNDARY_H__

#define HARD_WALL_VERTICAL 0

#define SOFT_WALL_VERTICAL 1
#define ALPHA_VERTICAL 0.3 // Soft wall parameter, 0 is no damping, 1 is linear damping, everything in between is exponential damping
#define SOFT_WALL_HEIGHT_PERCENTAGE_VERTICAL 0.02 // Percentage of the domain that will be damped at the top and bottom
#define SOFT_WALL_WIDTH 0.003 // Width of the damping region ( width of tanh function )


#define PERIODIC_BOUNDARY_VERTICAL 0
#define NO_BOUNDARY_VERTICAL 0

#define GHOST_CELLS_EXTRAPOLATION_VERTICAL 0 // 0 for constant, 1 for anti-symmetric

// For now this is used for sod shock tube
#define UPPER_PRESSURE_BOUNDARY 1.5e8
#define LOWER_PRESSURE_BOUNDARY -1.0e9
#define UPPER_ENTROPY_BOUNDARY -8.0e-3
#define LOWER_ENTROPY_BOUNDARY 1.0e-5

#if (HARD_WALL_VERTICAL + SOFT_WALL_VERTICAL + PERIODIC_BOUNDARY_VERTICAL + NO_BOUNDARY_VERTICAL) > 1
    #error "You can only have one vertical boundary type"
#endif

#endif // GLOBAL_BOUNDARY_H__