#ifndef FUNCTIONS_H__
#define FUNCTIONS_H__

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>

#include "global_parameters.h"
#include "global_constants.h"
#include "shared_files.h"

#include "./precalculated_variables/precalculated_variables.h"

#include "./solve_elliptic_equation/solve_elliptic_equation.h"
#include "./rhs_functions/rhs_functions.h"
#include "./one_time_step/one_time_step.h"
#include "./equations/equations.h"
#include "./boundary/boundary.h"
#include "./initialization_2D/initialization_2D.h"
#include "io_functions/io_functions.h"

void calculate_grid_info(struct GridInfo2D **grid_info, struct MpiInfo *mpi_info);
void calculate_grid_info_mpi(struct MpiInfo *mpi_info, struct GridInfo2D **grid_info);

#endif // FUNCTIONS_H__