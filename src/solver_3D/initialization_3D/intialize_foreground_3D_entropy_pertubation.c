#include "initialization_3D.h"

void initialize_foreground_3D_entropy_pertubation(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info)
{
    /*
    Initializes the grid with a small entropy pertubation and setting T1=0.

    Parameters
    ----------
    fg : ForegroundVariables3D
        A pointer to the ForegroundVariables3D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo2D
        A pointer to the GridInfo2D struct.
    */

    // Getting grid info
    int nz = grid_info->nz;
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    int nx = grid_info->nx;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    FLOAT_P amplitude = 1.0e2;
    FLOAT_P centre_z = 0.5*dz*nz;
    FLOAT_P sigma_z = 0.1*dz*nz;
    FLOAT_P centre_y = 0.5*dy*ny;
    FLOAT_P sigma_y = 0.1*dy*ny;
    FLOAT_P centre_x = 0.5*dx*nx;
    FLOAT_P sigma_x = 0.1*dx*nx;

    // First setting all values to zero
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                fg->p1[i][j][k] = 0.0;
                fg->s1[i][j][k] = 0.0;
                fg->T1[i][j][k] = 0.0;
                fg->vx[i][j][k] = 0.0;
                fg->vy[i][j][k] = 0.0;
                fg->vz[i][j][k] = 0.0;
            }
        }
    }

    // Spesific heat at constant pressure
    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    // Inside the grid
    for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                // Entropy pertubation
                fg->s1[i][j][k] = gaussian_3D((i-nz_ghost)*dz, j*dy, k*dz, centre_z, centre_y, centre_x, sigma_z, sigma_y, sigma_x, amplitude);
                // Calculating p1 from first law of thermodynamics
                fg->p1[i][j][k] = -bg->p0[i] * fg->s1[i][j][k]/c_p;
            }
            
        }
    }

    extrapolate_3D_array_down(fg->p1, nz_ghost, ny, nx);
    extrapolate_3D_array_up(fg->p1, nz_full, nz_ghost, ny, nx);
    extrapolate_3D_array_down(fg->s1, nz_ghost, ny, nx);
    extrapolate_3D_array_up(fg->s1, nz_full, nz_ghost, ny, nx);

    equation_of_state_3D(fg, bg, grid_info); // Getting rho1 from equation of state
}