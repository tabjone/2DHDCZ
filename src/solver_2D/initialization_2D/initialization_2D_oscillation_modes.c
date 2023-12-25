#include "initialization_2D.h"

void initialization_2D_oscillation_modes(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    With T=0
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P z0 = grid_info->z0;

    FLOAT_P Lz = R_SUN * (R_END - R_START);
    FLOAT_P Ly = R_SUN * Y_SIZE;

    initialize_foreground_zeros(fg, grid_info);

    FLOAT_P z, y;
    FLOAT_P sum_modes;

    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    int oscillation_modes_n[IC_OSCILLATION_MODES_N_NUM] = IC_OSCILLATION_MODES_N;
    int oscillation_modes_m[IC_OSCILLATION_MODES_M_NUM] = IC_OSCILLATION_MODES_M;
    FLOAT_P oscillation_modes_phi[IC_OSCILLATION_MODES_M_NUM] = IC_OSCILLATION_MODES_PHI;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        z = (i-nz_ghost) * dz + z0;
        for (int j = 0; j < ny; j++)
        {
            y = j * dy;
            sum_modes = 0.0;
            for (int k = 0; k < IC_OSCILLATION_MODES_N_NUM; k++)
            {
                for (int l = 0; l < IC_OSCILLATION_MODES_M_NUM; l++)
                {
                    sum_modes += sin(2.0 * M_PI * oscillation_modes_n[k] * z/Lz) * sin(2.0 * M_PI * oscillation_modes_m[l] * y/Ly + oscillation_modes_phi[l]);
                }
            }
            fg->s1[i][j] = 1.0e1 * sum_modes;
            // Calculating p1 from first law of thermodynamics
            fg->p1[i][j] = -bg->p0[i] * fg->s1[i][j]/c_p;
        }
    }
    update_vertical_boundary_ghostcells_2D(fg->p1, grid_info, mpi_info);
    update_vertical_boundary_ghostcells_2D(fg->s1, grid_info, mpi_info);

    equation_of_state(fg, bg, grid_info);  
}