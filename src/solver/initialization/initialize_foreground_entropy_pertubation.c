#include "initialization.h"

void initialize_foreground_entropy_pertubation(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info)
{
    /*
    Initializes the grid with a small entropy pertubation and setting T1=0.

    Parameters
    ----------
    fg : ForegroundVariables
        A pointer to the ForegroundVariables struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    // Getting grid info
    int nz = grid_info->nz;
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    FLOAT_P amplitude = 1.0e2;
    FLOAT_P centre_z = 0.5*dz*nz;
    FLOAT_P sigma_z = 0.1*dz*nz;
    FLOAT_P centre_y = 0.5*dy*ny;
    FLOAT_P sigma_y = 0.1*dy*ny;

    initialize_foreground_zeros(fg, grid_info); // Sets everything to zero so boundary and ghost cells are automatically zero

    // Spesific heat at constant pressure
    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    // Inside the grid
    for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            // Entropy pertubation
            fg->s1[i][j] = gaussian_2D((i-nz_ghost)*dz, j*dy, centre_z, centre_y, sigma_z, sigma_y, amplitude);;
            // Calculating p1 from first law of thermodynamics
            fg->p1[i][j] = -bg->p0[i] * fg->s1[i][j]/c_p;
        }
    }
    equation_of_state(fg, bg, grid_info); // Getting rho1 from equation of state
}