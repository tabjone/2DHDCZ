#include "data_structures/background_data/background_data.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_data_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_1D.h"
#include "data_structures/precalculated_data/precalculated_data_1D/precalculated_data_1D.h"
#include "solver/rhs_functions/rhs_functions_1D/rhs_functions_1D.h"
#include "global_float_precision.h"
#include "global_parameters.h"
#include <stdio.h>
#include <stdlib.h>

int test_rhs_ds1_dt_1D()
{
    // Test the rhs_ds1_dt_1D function

    // Declare the background variables, foreground variables and grid info
    struct BackgroundVariables *bg = (struct BackgroundVariables *)malloc(sizeof(struct BackgroundVariables));
    struct ForegroundVariables1D *fg = (struct ForegroundVariables1D *)malloc(sizeof(struct ForegroundVariables1D));
    struct GridInfo1D *grid_info = (struct GridInfo1D *)malloc(sizeof(struct GridInfo1D));
    struct PrecalculatedVariables1D *precalc = (struct PrecalculatedVariables1D *)malloc(sizeof(struct PrecalculatedVariables1D));

    int nz = 20;
    FLOAT_P dz = 1.0/(nz-1);

    FLOAT_P *vz, *s1, *grad_s0, *grad_T0;
    vz = malloc(nz*sizeof(FLOAT_P));
    s1 = malloc(nz*sizeof(FLOAT_P));
    grad_s0 = malloc(nz*sizeof(FLOAT_P));
    grad_T0 = malloc(nz*sizeof(FLOAT_P));
    fg->vz = vz;
    fg->s1 = s1;
    bg->grad_s0 = grad_s0;
    precalc->grad_T0 = grad_T0;

    for (int i = 0; i < nz; i++)
    {
        fg->vz[i] = i;
        fg->s1[i] = i*dz; // s1 = z => ds1/dz = 1
        bg->grad_s0[i] = 2.0;
    }

    precalc->one_over_dz = 1.0/dz;
    precalc->one_over_2dz = 1.0/(2*dz);

    FLOAT_P *rhs;
    rhs = malloc(nz*sizeof(FLOAT_P));

    for (int i = 0; i < nz; i++)
    {
        #if ADVECTION_ON == 1 && THERMAL_DIFFUSIVITY_ON == 0
            //analytical = - 3.0*vz[i];
        #endif // ADVECTION_ON
        #if ADVECTION_ON == 1 && THERMAL_DIFFUSIVITY_ON == 1
            //analytical = - 3.0*vz[i];// + 
        #endif // ADVECTION_ON
        //numerical = rhs_ds1_dt_1D(bg, fg, grid_info, precalc, i);
        


        //rhs[i]
    }

    // Deallocate memory
    free(vz);
    free(s1);
    free(grad_s0);
    free(grad_T0);
    free(rhs);
    free(bg);
    free(fg);
    free(grid_info);
    free(precalc);

    return 1;
}




void test_rhs_functions_1D()
{
    // Test the rhs functions in 1D

    test_rhs_ds1_dt_1D();
}
