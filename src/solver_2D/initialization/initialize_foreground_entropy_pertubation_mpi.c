#include "initialization.h"

void initialize_foreground_entropy_pertubation_mpi(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the grid with a small entropy pertubation and setting T1=0.

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
    int ny = grid_info->ny;
    FLOAT_P z_offset = grid_info->z_offset;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    FLOAT_P amplitude = 1.0e2;
    FLOAT_P centre_z = 0.5*dz*NZ;
    FLOAT_P sigma_z = 0.1*dz*NZ;
    FLOAT_P centre_y = 0.5*dy*ny;
    FLOAT_P sigma_y = 0.1*dy*ny;

    initialize_foreground_zeros(fg, grid_info); // Sets everything to zero so boundary and ghost cells are automatically zero

    FLOAT_P d[nz_full];
    calculate_damping(d, bg, grid_info);

    // Spesific heat at constant pressure
    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    // Inside the grid
    for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            // Entropy pertubation
            fg->s1[i][j] = d[i]*gaussian_2D((i-nz_ghost)*dz+z_offset, j*dy, centre_z, centre_y, sigma_z, sigma_y, amplitude);;
            // Calculating p1 from first law of thermodynamics
            fg->p1[i][j] = -bg->p0[i] * fg->s1[i][j]/c_p;
        }
    }

    #if VERTICAL_BOUNDARY_TYPE == 2
        periodic_boundary_2D(fg->p1, grid_info);
        periodic_boundary_2D(fg->s1, grid_info);
    #else
        extrapolate_2D_array_down(fg->p1, nz_ghost, ny);
        extrapolate_2D_array_up(fg->p1, nz_full, nz_ghost, ny);
        extrapolate_2D_array_down(fg->s1, nz_ghost, ny);
        extrapolate_2D_array_up(fg->s1, nz_full, nz_ghost, ny);
    #endif // VERTICAL_BOUNDARY_TYPE == 2

    equation_of_state(fg, bg, grid_info); // Getting rho1 from equation of state

    // Mpi
    communicate_2D_ghost_above_below(fg->p1, grid_info, mpi_info);
    communicate_2D_ghost_above_below(fg->s1, grid_info, mpi_info);
    communicate_2D_ghost_above_below(fg->rho1, grid_info, mpi_info);
}