#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/MPI_module.h"
#include "solver/equation_of_state/equation_of_state_3D/equation_of_state_3D.h"
#include "global_initialization.h"
#include "initialize_foreground_3D.h"
#include "global_constants.h"
#include "global_parameters.h"
#include <math.h>

void initialize_foreground_oscillation_modes_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    With T=0
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    int ny = grid_info->ny;
    int nx = grid_info->nx;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P z0 = grid_info->z0;

    FLOAT_P Lz = R_SUN * (R_END - R_START);
    FLOAT_P Ly = R_SUN * Y_SIZE;
    FLOAT_P Lx = R_SUN * X_SIZE;

    initialize_foreground_zeros_3D(fg, grid_info);

    FLOAT_P z, y, x;
    FLOAT_P sum_modes;

    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    int oscillation_modes_n[IC_OSCILLATION_MODES_N_NUM] = IC_OSCILLATION_MODES_N;
    int oscillation_modes_m[IC_OSCILLATION_MODES_M_NUM] = IC_OSCILLATION_MODES_M;
    int oscillation_modes_l[IC_OSCILLATION_MODES_L_NUM] = IC_OSCILLATION_MODES_L;
    FLOAT_P oscillation_modes_phi_m[IC_OSCILLATION_MODES_M_NUM] = IC_OSCILLATION_MODES_PHI_M;
    FLOAT_P oscillation_modes_phi_l[IC_OSCILLATION_MODES_L_NUM] = IC_OSCILLATION_MODES_PHI_L;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        z = (i-nz_ghost) * dz + z0;
        for (int j = 0; j < ny; j++)
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                sum_modes = 0.0;
                for (int n = 0; n < IC_OSCILLATION_MODES_N_NUM; n++)
                {
                    for (int m = 0; m < IC_OSCILLATION_MODES_M_NUM; m++)
                    {
                        for (int l = 0; l < IC_OSCILLATION_MODES_L_NUM; l++)
                        {
                            sum_modes += sin(2.0 * M_PI * oscillation_modes_n[n] * z/Lz) * sin(2.0 * M_PI * oscillation_modes_m[m] * y/Ly + oscillation_modes_phi_m[m]) * sin(2.0 * M_PI * oscillation_modes_l[l] * x/Lx + oscillation_modes_phi_l[l]);
                        }
                    }
                }
                fg->s1[i][j][k] = 1.0e1 * sum_modes;
                // Calculating p1 from first law of thermodynamics
                fg->p1[i][j][k] = -bg->p0[i] * fg->s1[i][j][k]/c_p;
            }
        }   
    communicate_3D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost, ny, nx);
    communicate_3D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost, ny, nx);

    equation_of_state_3D(fg, bg, grid_info);  
}