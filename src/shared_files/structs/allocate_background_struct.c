#include "structs.h"
#include <stdlib.h>
#include <stdio.h>

void allocate_background_struct(struct BackgroundVariables **bg, int nz_full)
{
    *bg = (struct BackgroundVariables *)malloc(sizeof(struct BackgroundVariables));

    allocate_1D_array(&(*bg)->r, nz_full);
    allocate_1D_array(&(*bg)->T0, nz_full);
    allocate_1D_array(&(*bg)->rho0, nz_full);
    allocate_1D_array(&(*bg)->one_over_rho0, nz_full);
    allocate_1D_array(&(*bg)->p0, nz_full);
    allocate_1D_array(&(*bg)->grad_s0, nz_full);
    allocate_1D_array(&(*bg)->g, nz_full);
    allocate_1D_array(&(*bg)->grad_g, nz_full);
    allocate_1D_array(&(*bg)->grad_rho0, nz_full);
    allocate_1D_array(&(*bg)->eta_over_four_pi_rho0_T0, nz_full);
}