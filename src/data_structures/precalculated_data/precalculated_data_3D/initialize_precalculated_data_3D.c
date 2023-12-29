#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "precalculated_data_struct_3D.h"
#include "global_float_precision.h"
#include "global_constants.h"

void initialize_precalculated_data_3D(struct PrecalculatedVariables3D *pv, struct BackgroundVariables *bg, struct GridInfo3D *grid_info)
{
    int nz_full = grid_info->nz_full;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dx = grid_info->dx;

    // Pre-calculate 1/rho0 and eta/(4*pi*rho0*T0), and one_over_rho0_T0_const
    for (int i = 0; i < nz_full; i++)
    {
        pv->one_over_rho0[i] = 1.0/bg->rho0[i];
        pv->eta_over_four_pi_rho0_T0[i] = ETA/(4*M_PI*bg->rho0[i]*bg->T0[i]);
        pv->THERM_COEFF_over_T0_rho0[i] = THERMAL_DIFFUSIVITY_COEFF/(bg->T0[i]*bg->rho0[i]);

        pv->VIS_COEFF_over_T0_rho0[i] = VISCOSITY_COEFF/(bg->T0[i]*bg->rho0[i]);
        pv->VIS_COEFF_over_rho0[i] = VISCOSITY_COEFF/bg->rho0[i];
    }

    // Pre-calculate gradients
    // First for endpoints
    pv->grad_g[0] = (bg->g[1]-bg->g[0])/dz;
    pv->grad_rho0[0] = (bg->rho0[1]-bg->rho0[0])/dz;
    pv->grad_T0[0] = (bg->T0[1]-bg->T0[0])/dz;

    pv->grad_g[nz_full-1] = (bg->g[nz_full-1]-bg->g[nz_full-2])/dz;
    pv->grad_rho0[nz_full-1] = (bg->rho0[nz_full-1]-bg->rho0[nz_full-2])/dz;
    pv->grad_T0[nz_full-1] = (bg->T0[nz_full-1]-bg->T0[nz_full-2])/dz;

    // Then for the rest of the grid using central derivative second order
    for (int i = 1; i < nz_full-1; i++)
    {
        pv->grad_g[i] = (bg->g[i+1]-bg->g[i-1])/(2.0*dz);
        pv->grad_rho0[i] = (bg->rho0[i+1]-bg->rho0[i-1])/(2.0*dz);
        pv->grad_T0[i] = (bg->T0[i+1]-bg->T0[i-1])/(2.0*dz);
    }

    pv->two_VIS_COEFF = 2.0*VISCOSITY_COEFF;
    
    pv->one_over_8dxdxdy = 1.0/(8.0*dx*dx*dy);
    pv->one_over_8dxdydy = 1.0/(8.0*dx*dy*dy);
    pv->one_over_8dxdydz = 1.0/(8.0*dx*dy*dz);
    pv->one_over_8dxdzdz = 1.0/(8.0*dx*dz*dz);
    pv->one_over_8dydydz = 1.0/(8.0*dy*dy*dz);
    pv->one_over_8dydzdz = 1.0/(8.0*dy*dz*dz);
    pv->one_over_2dx = 1.0/(2.0*dx);
    pv->one_over_2dy = 1.0/(2.0*dy);
    pv->one_over_2dz = 1.0/(2.0*dz);
    pv->one_over_dxdx = 1.0/(dx*dx);
    pv->one_over_dydy = 1.0/(dy*dy);
    pv->one_over_dzdz = 1.0/(dz*dz);
    pv->one_over_dx = 1.0/dx;
    pv->one_over_dy = 1.0/dy;
    pv->one_over_dz = 1.0/dz;
    pv->one_over_4dxdy = 1.0/(4.0*dx*dy);
    pv->one_over_4dydz = 1.0/(4.0*dy*dz);
    pv->one_over_4dxdz = 1.0/(4.0*dx*dz);
}