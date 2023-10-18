#include "initialization.h"

void initialize_foreground_density_pertubation(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info)
{
    /*
    Initializes the foreground struct with a density pertubation with entropy set to zero.

    Parameters
    ----------
    fg : ForegroundVariables
        A pointer to the ForegroundVariables struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    mpi_info : MpiInfo
        A pointer to the MpiInfo struct.
    */

    // Setting properties of gaussian and getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    FLOAT_P dz = grid_info->dz;
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    
    // Gaussian parameters
    FLOAT_P amplitude = 1.0e-5;
    FLOAT_P centre_z = 0.5*dz*nz;
    FLOAT_P sigma_z = 0.1*dz*nz;
    FLOAT_P centre_y = 0.5*dy*ny;
    FLOAT_P sigma_y = 0.1*dy*ny;

    initialize_foreground_zeros(fg, grid_info); // Sets everything to zero so boundary and ghost cells are automatically zero
    
    // Inside the grid
    for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->rho1[i][j] = gaussian_2D((i-nz_ghost)*dz, j*dy, centre_z, centre_y, sigma_z, sigma_y, amplitude);
            fg->p1[i][j] = GAMMA*bg->p0[i]*fg->rho1[i][j]/bg->rho0[i];
            fg->T1[i][j] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i][j]/bg->p0[i]);
        }
    }

    #if VERTICAL_BOUNDARY_TYPE == 2
        periodic_boundary_2D(fg->rho1, grid_info);
        periodic_boundary_2D(fg->p1, grid_info);
        periodic_boundary_2D(fg->T1, grid_info);
    #endif // VERTICAL_BOUNDARY_TYPE == 2
}