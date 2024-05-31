#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/MPI_module.h"
#include "solver/equation_of_state/equation_of_state_2D/equation_of_state_2D.h"
#include "global_initialization.h"
#include "initialize_foreground_2D.h"
#include "global_constants.h"
#include "global_parameters.h"
#include <math.h>

void initialize_foreground_oscillation_modes_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    With T=0
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;

    FLOAT_P Lz = R_SUN * (R_END - R_START);
    #if COORDINATES == 0
        FLOAT_P Ly = R_SUN * Y_SIZE;
    #elif COORDINATES == 1
        FLOAT_P Ly = Y_SIZE;
    #endif

    initialize_foreground_zeros_2D(fg, grid_info);

    FLOAT_P z, y;
    FLOAT_P sum_modes;

    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    int oscillation_modes_n[IC_OSCILLATION_MODES_N_NUM] = IC_OSCILLATION_MODES_N;
    int oscillation_modes_m[IC_OSCILLATION_MODES_M_NUM] = IC_OSCILLATION_MODES_M;
    FLOAT_P oscillation_modes_phi_m[IC_OSCILLATION_MODES_M_NUM] = IC_OSCILLATION_MODES_PHI_M;


    FLOAT_P max_s1 = 0.0;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        z = bg->r[i] - R_START*R_SUN;
        for (int j = 0; j < ny; j++)
        {
            y = j * dy;
            sum_modes = 0.0;
            for (int k = 0; k < IC_OSCILLATION_MODES_N_NUM; k++)
            {
                for (int l = 0; l < IC_OSCILLATION_MODES_M_NUM; l++)
                {
                    sum_modes += sin(2.0 * M_PI * oscillation_modes_n[k] * z/Lz) * sin(2.0 * M_PI * oscillation_modes_m[l] * y/Ly + oscillation_modes_phi_m[l]);
                }
            }
            fg->s1[i][j] = sum_modes;
            if (fg->s1[i][j] > max_s1)
            {
                max_s1 = fabs(fg->s1[i][j]);
            }
        }
    }

    FLOAT_P max_s1_global;
    MPI_Allreduce(&max_s1, &max_s1_global, 1, MPI_FLOAT_P, MPI_MAX, MPI_COMM_WORLD);
    max_s1 = max_s1_global;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            // Normalizing s1
            fg->s1[i][j] = fg->s1[i][j] / max_s1 * 0.001 * c_p;
            // Calculating p1 from first law of thermodynamics
            fg->p1[i][j] = -bg->p0[i] * fg->s1[i][j]/c_p;
        }
    }

    communicate_2D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost, ny);

    equation_of_state_2D(fg, bg, grid_info, mpi_info); 
}