#include "initialization_2D.h"

void sod_shock_horizontal_mpi(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the foreground struct with a Sod Shock Tube test.

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

    FLOAT_P left_rho, left_p;
    if (mpi_info->rank == 0)
    {
        printf("Initializing Sod Shock Tube\n");
        // For now just letting the left side be one of the top values from process 1
        left_rho = bg->rho0[nz] * 1.0e-8;
        left_p = bg->p0[nz] * 1.0e-8;
    }
    MPI_Bcast(&left_rho, 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);
    MPI_Bcast(&left_p, 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);

    FLOAT_P right_rho = left_rho * 0.5;
    FLOAT_P right_p = left_p * 0.5;

    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    initialize_foreground_zeros(fg, grid_info);

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (j < ny/2)
            {
                fg->p1[i][j] = left_p;
                fg->rho1[i][j] = left_rho;
            }
            else
            {
                fg->p1[i][j] = right_p;
                fg->rho1[i][j] = right_rho;
            }
            // Getting T1 from equation of state
            fg->T1[i][j] = bg->T0[i] * (fg->p1[i][j]/bg->p0[i] - fg->rho1[i][j]/bg->rho0[i]);
            // Getting entropy from first law of thermodynamics
            fg->s1[i][j] = c_p * (fg->T1[i][j]/bg->T0[i] - fg->p1[i][j]/bg->p0[i]);
        }
    }

    #if VERTICAL_BOUNDARY_TYPE == 2
        update_vertical_boundary_ghostcells_2D(fg->p1, grid_info, mpi_info);
        update_vertical_boundary_ghostcells_2D(fg->rho1, grid_info, mpi_info);
        update_vertical_boundary_ghostcells_2D(fg->T1, grid_info, mpi_info);
        update_vertical_boundary_ghostcells_2D(fg->s1, grid_info, mpi_info);
        update_vertical_boundary_ghostcells_2D(fg->vy, grid_info, mpi_info);
        update_vertical_boundary_ghostcells_2D(fg->vz, grid_info, mpi_info);
    #else
        // Setting top and bottom boundaries to zero
        // Not yet implemented
    #endif // VERTICAL_BOUNDARY_TYPE
}