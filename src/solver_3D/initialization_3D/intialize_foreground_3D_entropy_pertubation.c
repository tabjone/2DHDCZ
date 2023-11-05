#include "initialization_3D.h"

void initialize_foreground_3D_entropy_pertubation(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info)
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

    // Getting mpi info
    FLOAT_P z_offset = grid_info->z_offset;

    // Setting up gaussians
    FLOAT_P amplitude = 1.0e2;
    FLOAT_P sigma_z = 0.1*dz*nz;
    FLOAT_P sigma_y = 0.1*dy*ny;
    FLOAT_P sigma_x = 0.1*dx*nx;

    // Getting centre of gaussians from initial conditions file
    FLOAT_P centre_z[IC_N_ENTROPY_PERTUBATION] = IC_ENTROPY_CENTRE_Z;
    FLOAT_P centre_y[IC_N_ENTROPY_PERTUBATION] = IC_ENTROPY_CENTRE_Y;
    FLOAT_P centre_x[IC_N_ENTROPY_PERTUBATION] = IC_ENTROPY_CENTRE_X;

    for (int i = 0; i < IC_N_ENTROPY_PERTUBATION; i++)
    {
        centre_z[i] *= dz * NZ;
        centre_y[i] *= dy * NY;
        centre_x[i] *= dx * NX;
    }

    initialize_foreground_zeros_3D(fg, grid_info); // Sets everything to zero so boundary and ghost cells are automatically zero

    // Spesific heat at constant pressure
    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    // Inside the grid
    for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                for (int n = 0; n < IC_N_ENTROPY_PERTUBATION; n++)
                {
                    // Entropy pertubation
                    fg->s1[i][j][k] += gaussian_3D((i-nz_ghost)*dz+z_offset, j*dy, k*dz, centre_z[n], centre_y[n], centre_x[n], sigma_z, sigma_y, sigma_x, amplitude);
                }
                // Calculating p1 from first law of thermodynamics
                fg->p1[i][j][k] = -bg->p0[i] * fg->s1[i][j][k]/c_p;
            }
            
        }
    }

    update_vertical_boundary_ghostcells_3D(fg->s1, grid_info, mpi_info);
    update_vertical_boundary_ghostcells_3D(fg->p1, grid_info, mpi_info);

    equation_of_state_3D(fg, bg, grid_info); // Getting rho1 from equation of state
}