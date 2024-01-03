#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/MPI_module.h"
#include "solver/equation_of_state/equation_of_state_2D/equation_of_state_2D.h"
#include "initialize_foreground_2D.h"
#include "global_initialization.h"
#include "global_parameters.h"
#include "global_constants.h"
#include <math.h>

FLOAT_P gaussian_2D(FLOAT_P x, FLOAT_P y, FLOAT_P x0, FLOAT_P y0, FLOAT_P sigma_x, FLOAT_P sigma_y, FLOAT_P A) 
{
    return A * exp(-(x - x0) * (x - x0) / (2 * sigma_x * sigma_x) 
                   -(y - y0) * (y - y0) / (2 * sigma_y * sigma_y));
}

void initialize_foreground_entropy_perturbations_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the grid with a small entropy perturbation and setting T1=0.

    Parameters
    ----------
    fg : ForegroundVariables2D
        A pointer to the ForegroundVariables2D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo2D
        A pointer to the GridInfo2D struct.
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    int ny = grid_info->ny;
    FLOAT_P z_offset = grid_info->z_offset;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    FLOAT_P amplitude = 1.0e1;
    FLOAT_P sigma_z = 0.1*dz*NZ;
    FLOAT_P sigma_y = 0.1*dy*NY;

    FLOAT_P centre_z[IC_N_ENTROPY_PERTUBATION] = IC_ENTROPY_CENTRE_Z;
    FLOAT_P centre_y[IC_N_ENTROPY_PERTUBATION] = IC_ENTROPY_CENTRE_Y;

    for (int i = 0; i < IC_N_ENTROPY_PERTUBATION; i++)
    {
        centre_z[i] *= dz * NZ;
        centre_y[i] *= dy * NY;
    }

    initialize_foreground_zeros_2D(fg, grid_info); // Sets everything to zero so boundary and ghost cells are automatically zero

    // Spesific heat at constant pressure
    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    // Inside the grid
    for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int n = 0; n < IC_N_ENTROPY_PERTUBATION; n++)
            {
                // Entropy perturbation
                fg->s1[i][j] += gaussian_2D((i-nz_ghost)*dz+z_offset, j*dy, centre_z[n], centre_y[n], sigma_z, sigma_y, amplitude);
            }

            // Calculating p1 from first law of thermodynamics
            fg->p1[i][j] = -bg->p0[i] * fg->s1[i][j]/c_p;
        }
    }

    communicate_2D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost, ny);

    equation_of_state_2D(fg, bg, grid_info); // Getting rho1 from equation of state
}