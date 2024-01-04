#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/MPI_module.h"
#include "solver/equation_of_state/equation_of_state_1D/equation_of_state_1D.h"
#include "initialize_foreground_1D.h"
#include "global_initialization.h"
#include "global_parameters.h"
#include "global_constants.h"
#include <math.h>

FLOAT_P gaussian_1D(FLOAT_P x, FLOAT_P x0, FLOAT_P sigma_x, FLOAT_P A)
{
    return A * exp(-(x - x0) * (x - x0) / (2 * sigma_x * sigma_x));
}

void initialize_foreground_entropy_perturbations_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the grid with a small entropy perturbation and setting T1=0.

    Parameters
    ----------
    fg : ForegroundVariables1D
        A pointer to the ForegroundVariables1D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo1D
        A pointer to the GridInfo1D struct.
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    FLOAT_P z_offset = grid_info->z_offset;
    FLOAT_P dz = grid_info->dz;

    FLOAT_P amplitude = 1.0e1;
    FLOAT_P sigma_z = 0.1*dz*NZ;

    FLOAT_P centre_z[IC_N_ENTROPY_PERTURBATION] = IC_ENTROPY_CENTRE_Z;

    for (int i = 0; i < IC_N_ENTROPY_PERTURBATION; i++)
    {
        centre_z[i] *= dz * NZ;
    }

    initialize_foreground_zeros_1D(fg, grid_info); // Sets everything to zero so boundary and ghost cells are automatically zero

    // Spesific heat at constant pressure
    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    // Inside the grid
    for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
    {
        for (int n = 0; n < IC_N_ENTROPY_PERTURBATION; n++)
        {
            // Entropy perturbation
            fg->s1[i] += gaussian_1D((i-nz_ghost)*dz+z_offset, centre_z[n], sigma_z, amplitude);
        }

        // Calculating p1 from first law of thermodynamics
        fg->p1[i] = -bg->p0[i] * fg->s1[i]/c_p;
    }

    communicate_1D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost);
    communicate_1D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost);

    equation_of_state_1D(fg, bg, grid_info); // Getting rho1 from equation of state
}