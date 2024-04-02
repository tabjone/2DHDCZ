#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "precalculated_data_struct_2D.h"
#include "global_float_precision.h"
#include "global_constants.h"
#include "solver/boundary/boundary_2D/boundary_2D.h"
#include "MPI_module/mpi_info_struct.h"

void initialize_precalculated_data_2D(struct PrecalculatedVariables2D *pv, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    int nz_full = grid_info->nz_full;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P dy = grid_info->dy;

    calculate_damping_2D(pv->damping_factor, bg, grid_info, mpi_info);

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
    
    pv->one_over_8dydydz = 1.0/(8.0*dy*dy*dz);
    pv->one_over_8dydzdz = 1.0/(8.0*dy*dz*dz);
    pv->one_over_2dydydy = 1.0/(2.0*dy*dy*dy);
    pv->one_over_2dzdzdz = 1.0/(2.0*dz*dz*dz);
    pv->one_over_2dy = 1.0/(2.0*dy);
    pv->one_over_2dz = 1.0/(2.0*dz);
    pv->one_over_dydy = 1.0/(dy*dy);
    pv->one_over_dzdz = 1.0/(dz*dz);
    pv->one_over_dz = 1.0/dz;
    pv->one_over_dy = 1.0/dy;
    pv->one_over_4dydz = 1.0/(4.0*dy*dz);
}