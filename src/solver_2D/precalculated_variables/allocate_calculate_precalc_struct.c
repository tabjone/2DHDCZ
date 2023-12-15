#include "precalculated_variables.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void allocate_calculate_precalc_struct(struct PrecalculatedVariables **pv, struct BackgroundVariables *bg, struct GridInfo2D *grid_info)
{
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P dy = grid_info->dy;

    // Allocate memory for the struct
    *pv = (struct PrecalculatedVariables *)malloc(sizeof(struct PrecalculatedVariables));
    // HANDLE MALLOC FAILURE

    // Allocate memory for the arrays
    allocate_1D_array(&(*pv)->one_over_rho0, nz_full);
    allocate_1D_array(&(*pv)->grad_g, nz_full);
    allocate_1D_array(&(*pv)->grad_rho0, nz_full);
    allocate_1D_array(&(*pv)->eta_over_four_pi_rho0_T0, nz_full);
    allocate_1D_array(&(*pv)->grad_T0, nz_full);
    allocate_1D_array(&(*pv)->VIS_COEFF_over_rho0, nz_full);
    allocate_1D_array(&(*pv)->VIS_COEFF_over_T0_rho0, nz_full);
    allocate_1D_array(&(*pv)->THERM_COEFF_over_T0_rho0, nz_full);

    (*pv)->j_plus = (int *)malloc(ny * sizeof(int));
    (*pv)->j_minus = (int *)malloc(ny * sizeof(int));
    (*pv)->j_plus2 = (int *)malloc(ny * sizeof(int));
    (*pv)->j_minus2 = (int *)malloc(ny * sizeof(int));
    // HANDLE MALLOC FAILURE

    // Pre-calculate 1/rho0 and eta/(4*pi*rho0*T0), and one_over_rho0_T0_const
    for (int i = 0; i < nz_full; i++)
    {
        (*pv)->one_over_rho0[i] = 1.0/bg->rho0[i];
        (*pv)->eta_over_four_pi_rho0_T0[i] = ETA/(4*M_PI*bg->rho0[i]*bg->T0[i]);
        (*pv)->THERM_COEFF_over_T0_rho0[i] = THERMAL_DIFFUSIVITY_COEFF/(bg->T0[i]*bg->rho0[i]);

        (*pv)->VIS_COEFF_over_T0_rho0[i] = VISCOSITY_COEFF/(bg->T0[i]*bg->rho0[i]);
        (*pv)->VIS_COEFF_over_rho0[i] = VISCOSITY_COEFF/bg->rho0[i];

    }

    // Pre-calculate gradients
    // First for endpoints
    (*pv)->grad_g[0] = (bg->g[1]-bg->g[0])/dz;
    (*pv)->grad_rho0[0] = (bg->rho0[1]-bg->rho0[0])/dz;
    (*pv)->grad_T0[0] = (bg->T0[1]-bg->T0[0])/dz;

    (*pv)->grad_g[nz_full-1] = (bg->g[nz_full-1]-bg->g[nz_full-2])/dz;
    (*pv)->grad_rho0[nz_full-1] = (bg->rho0[nz_full-1]-bg->rho0[nz_full-2])/dz;
    (*pv)->grad_T0[nz_full-1] = (bg->T0[nz_full-1]-bg->T0[nz_full-2])/dz;

    // Then for the rest of the grid using central derivative second order
    for (int i = 1; i < nz_full-1; i++)
    {
        (*pv)->grad_g[i] = (bg->g[i+1]-bg->g[i-1])/(2.0*dz);
        (*pv)->grad_rho0[i] = (bg->rho0[i+1]-bg->rho0[i-1])/(2.0*dz);
        (*pv)->grad_T0[i] = (bg->T0[i+1]-bg->T0[i-1])/(2.0*dz);
    }

    (*pv)->two_VIS_COEFF = 2.0*VISCOSITY_COEFF;
    
    (*pv)->one_over_8dydydz = 1.0/(8.0*dy*dy*dz);
    (*pv)->one_over_8dydzdz = 1.0/(8.0*dy*dz*dz);
    (*pv)->one_over_2dy = 1.0/(2.0*dy);
    (*pv)->one_over_2dz = 1.0/(2.0*dz);
    (*pv)->one_over_dydy = 1.0/(dy*dy);
    (*pv)->one_over_dzdz = 1.0/(dz*dz);
    (*pv)->one_over_dz = 1.0/dz;
    (*pv)->one_over_dy = 1.0/dy;
    (*pv)->one_over_4dydz = 1.0/(4.0*dy*dz);

    // Precalculate periodic boundary indices
    for (int j = 0; j < ny; j++)
    {
        (*pv)->j_plus[j] = periodic_boundary(j+1, ny);
        (*pv)->j_minus[j] = periodic_boundary(j-1, ny);
        (*pv)->j_plus2[j] = periodic_boundary(j+2, ny);
        (*pv)->j_minus2[j] = periodic_boundary(j-2, ny);
    }
}