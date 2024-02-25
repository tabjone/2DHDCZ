#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/MPI_module.h"
#include "solver/equation_of_state/equation_of_state_1D/equation_of_state_1D.h"
#include "global_initialization.h"
#include "initialize_foreground_1D.h"
#include "global_constants.h"
#include "global_parameters.h"
#include <math.h>

void initialize_foreground_oscillation_modes_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    With T=0
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P z0 = grid_info->z0;

    FLOAT_P Lz = R_SUN * (R_END - R_START);

    initialize_foreground_zeros_1D(fg, grid_info);

    FLOAT_P z;
    FLOAT_P sum_modes;

    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    int oscillation_modes_n[IC_OSCILLATION_MODES_N_NUM] = IC_OSCILLATION_MODES_N;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        z = (i-nz_ghost) * dz + z0;
        sum_modes = 0.0;
        for (int k = 0; k < IC_OSCILLATION_MODES_N_NUM; k++)
        { 
            sum_modes += sin(2.0 * M_PI * oscillation_modes_n[k] * (z-z0)/Lz);
        }
        fg->s1[i] = 1.0e1 * sum_modes;
        // Calculating p1 from first law of thermodynamics
        fg->p1[i] = -bg->p0[i] * fg->s1[i]/c_p;
    }

    communicate_1D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost);
    communicate_1D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost);

    equation_of_state_1D(fg, bg, grid_info);  
}