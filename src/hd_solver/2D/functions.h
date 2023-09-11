#ifndef FUNCTIONS_H__
#define FUNCTIONS_H__

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>

#include "./solve_elliptic_equation/solve_elliptic_equation.h"
#include "./solve_diff_eqs/solve_diff_eqs.h"

#include "./initialization/initialization.h"
#include "shared_files.h"
#include "global_parameters.h"
#include "global_constants.h"


#include "io_functions/io_functions.h"


void one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg);

#endif // FUNCTIONS_H__