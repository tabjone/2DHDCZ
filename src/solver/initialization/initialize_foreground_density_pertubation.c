#include "initialization.h"

void initialize_foreground_density_pertubation(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the foreground struct with a density pertubation.

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

    bool has_neighbor_above = mpi_info->has_neighbor_above;
    bool has_neighbor_below = mpi_info->has_neighbor_below;

    initialize_foreground_zeros(fg, grid_info); // Sets everything to zero so boundary and ghost cells are automatically zero

    #if DIMENSIONS == 1
        // Getting grid info
        int nz = grid_info->nz;
        int nz_full = grid_info->nz_full;
        int nz_ghost = grid_info->nz_ghost;
        double dz = grid_info->dz;

        // Inside the grid
        for (int i = nz_ghost; i < nz_full-nz_ghost)
        {
            fg->rho1[i] = gaussian_1D((i-nz_ghost)*dz, 0.5*dz*nz, 0.1*dz*nz, 1.0e-5);
            fg->p1[i] = GAMMA*bg->p0[i]*fg->rho1[i]/bg->rho0[i];
            fg->T1[i] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i]/bg->p0[i]);
        }

        // Top boundary
        if (has_neighbor_above) // This is false for MPI_ON = 0
        {
            // Boundary
            fg->rho1[nz_full-nz_ghost-1] = gaussian_1D((nz_full-nz_ghost-1-nz_ghost)*dz, 0.5*dz*nz, 0.1*dz*nz, 1.0e-5);
            fg->p1[nz_full-nz_ghost-1] = GAMMA*bg->p0[nz_full-nz_ghost-1]*fg->rho1[nz_full-nz_ghost-1]/bg->rho0[nz_full-nz_ghost-1];
            fg->T1[nz_full-nz_ghost-1] = bg->T0[nz_full-nz_ghost-1]*((GAMMA-1)/GAMMA * fg->p1[nz_full-nz_ghost-1]/bg->p0[nz_full-nz_ghost-1]);

            // Calculating ghost cells from neighbor above
            for (int i = nz_full-nz_ghost; i < nz_full; i++)
            {
                fg->rho1[i] = gaussian_1D((i-nz_ghost)*dz, 0.5*dz*nz, 0.1*dz*nz, 1.0e-5);
                fg->p1[i] = GAMMA*bg->p0[i]*fg->rho1[i]/bg->rho0[i];
                fg->T1[i] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i]/bg->p0[i]);
            }
        }

        // Bottom boundary
        if (has_neighbor_below) // This is false for MPI_ON = 0
        {
            // Boundary
            fg->rho1[nz_ghost] = gaussian_1D(0.0, 0.5*dz*nz, 0.1*dz*nz, 1.0e-5);
            fg->p1[nz_ghost] = GAMMA*bg->p0[nz_ghost]*fg->rho1[nz_ghost]/bg->rho0[nz_ghost];
            fg->T1[nz_ghost] = bg->T0[nz_ghost]*((GAMMA-1)/GAMMA * fg->p1[nz_ghost]/bg->p0[nz_ghost]);

            // Calculating ghost cells from neighbor below
            for (int i = 0; i < nz_ghost; i++)
            {
                fg->rho1[i] = gaussian_1D((i-nz_ghost)*dz, 0.5*dz*nz, 0.1*dz*nz, 1.0e-5);
                fg->p1[i] = GAMMA*bg->p0[i]*fg->rho1[i]/bg->rho0[i];
                fg->T1[i] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i]/bg->p0[i]);
            }
        }

    #elif DIMENSIONS == 2
        // Getting grid info
        int ny = grid_info->ny;
        int nz = grid_info->nz;
        int nz_ghost = grid_info->nz_ghost;
        int nz_full = grid_info->nz_full;
        double dz = grid_info->dz;
        double dy = grid_info->dy;

        // Inside the grid
        for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                fg->rho1[i][j] = gaussian_2D((i-nz_ghost)*dz, j*dy, 0.5*dz*nz, 0.5*dy*ny, 0.1*dz*nz, 0.1*dy*ny, 1.0e-5);
                fg->p1[i][j] = GAMMA*bg->p0[i]*fg->rho1[i][j]/bg->rho0[i];
                fg->T1[i][j] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i][j]/bg->p0[i]);
            }
        }

        // Top boundary
        if (has_neighbor_above) // This is false for MPI_ON = 0
        {
            // Boundary
            for (int j = 0; j < ny; j++)
            {
                fg->rho1[nz_full-nz_ghost-1][j] = gaussian_2D((nz_full-nz_ghost-1-nz_ghost)*dz, j*dy, 0.5*dz*nz, 0.5*dy*ny, 0.1*dz*nz, 0.1*dy*ny, 1.0e-5);
                fg->p1[nz_full-nz_ghost-1][j] = GAMMA*bg->p0[nz_full-nz_ghost-1]*fg->rho1[nz_full-nz_ghost-1][j]/bg->rho0[nz_full-nz_ghost-1];
                fg->T1[nz_full-nz_ghost-1][j] = bg->T0[nz_full-nz_ghost-1]*((GAMMA-1)/GAMMA * fg->p1[nz_full-nz_ghost-1][j]/bg->p0[nz_full-nz_ghost-1]);
            }

            // Calculating ghost cells from neighbor above
            for (int i = nz_full-nz_ghost; i < nz_full; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    fg->rho1[i][j] = gaussian_2D((i-nz_ghost)*dz, j*dy, 0.5*dz*nz, 0.5*dy*ny, 0.1*dz*nz, 0.1*dy*ny, 1.0e-5);
                    fg->p1[i][j] = GAMMA*bg->p0[i]*fg->rho1[i][j]/bg->rho0[i];
                    fg->T1[i][j] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i][j]/bg->p0[i]);
                }
            }
        }

        // Bottom boundary
        if (has_neighbor_below) // This is false for MPI_ON = 0
        {
            // Boundary
            for (int j = 0; j < ny; j++)
            {
                fg->rho1[nz_ghost][j] = gaussian_2D(0.0, j*dy, 0.5*dz*nz, 0.5*dy*ny, 0.1*dz*nz, 0.1*dy*ny, 1.0e-5);
                fg->p1[nz_ghost][j] = GAMMA*bg->p0[nz_ghost]*fg->rho1[nz_ghost][j]/bg->rho0[nz_ghost];
                fg->T1[nz_ghost][j] = bg->T0[nz_ghost]*((GAMMA-1)/GAMMA * fg->p1[nz_ghost][j]/bg->p0[nz_ghost]);
            }

            // Calculating ghost cells from neighbor below
            for (int i = 0; i < nz_ghost; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    fg->rho1[i][j] = gaussian_2D((i-nz_ghost)*dz, j*dy, 0.5*dz*nz, 0.5*dy*ny, 0.1*dz*nz, 0.1*dy*ny, 1.0e-5);
                    fg->p1[i][j] = GAMMA*bg->p0[i]*fg->rho1[i][j]/bg->rho0[i];
                    fg->T1[i][j] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i][j]/bg->p0[i]);
                }
            }
        }

    #elif DIMENSIONS == 3
        // Getting grid info
        int ny = grid_info->ny;
        int nx = grid_info->nx;
        int nz = grid_info->nz;
        int nz_ghost = grid_info->nz_ghost;
        int nz_full = grid_info->nz_full;
        double dz = grid_info->dz;
        double dy = grid_info->dy;
        double dx = grid_info->dx;

        // Inside the grid
        for (int i = nz_ghost; i < nz_full-nz_ghost; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    fg->rho1[i][j][k] = gaussian_3D((i-nz_ghost)*dz, j*dy, k*dx, 0.5*dz*nz, 0.5*dy*ny, 0.5*dx*nx, 0.1*dz*nz, 0.1*dy*ny, 0.1*dx*nx, 1.0e-5);
                    fg->p1[i][j][k] = GAMMA*bg->p0[i]*fg->rho1[i][j][k]/bg->rho0[i];
                    fg->T1[i][j][k] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i][j][k]/bg->p0[i]);
                }
            }
        }

        // Top boundary
        if (has_neighbor_above) // This is false for MPI_ON = 0
        {
            // Boundary
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    fg->rho1[nz_full-nz_ghost-1][j][k] = gaussian_3D((nz_full-nz_ghost-1-nz_ghost)*dz, j*dy, k*dx, 0.5*dz*nz, 0.5*dy*ny, 0.5*dx*nx, 0.1*dz*nz, 0.1*dy*ny, 0.1*dx*nx, 1.0e-5);
                    fg->p1[nz_full-nz_ghost-1][j][k] = GAMMA*bg->p0[nz_full-nz_ghost-1]*fg->rho1[nz_full-nz_ghost-1][j][k]/bg->rho0[nz_full-nz_ghost-1];
                    fg->T1[nz_full-nz_ghost-1][j][k] = bg->T0[nz_full-nz_ghost-1]*((GAMMA-1)/GAMMA * fg->p1[nz_full-nz_ghost-1][j][k]/bg->p0[nz_full-nz_ghost-1]);
                }
            }

            // Calculating ghost cells from neighbor above
            for (int i = nz_full-nz_ghost; i < nz_full; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int k = 0; k < nx; k++)
                    {
                        fg->rho1[i][j][k] = gaussian_3D((i-nz_ghost)*dz, j*dy, k*dx, 0.5*dz*nz, 0.5*dy*ny, 0.5*dx*nx, 0.1*dz*nz, 0.1*dy*ny, 0.1*dx*nx, 1.0e-5);
                        fg->p1[i][j][k] = GAMMA*bg->p0[i]*fg->rho1[i][j][k]/bg->rho0[i];
                        fg->T1[i][j][k] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i][j][k]/bg->p0[i]);
                    }
                }
            }
        }

        // Bottom boundary
        if (has_neighbor_below) // This is false for MPI_ON = 0
        {
            // Boundary
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    fg->rho1[nz_ghost][j][k] = gaussian_3D(0.0, j*dy, k*dx, 0.5*dz*nz, 0.5*dy*ny, 0.5*dx*nx, 0.1*dz*nz, 0.1*dy*ny, 0.1*dx*nx, 1.0e-5);
                    fg->p1[nz_ghost][j][k] = GAMMA*bg->p0[nz_ghost]*fg->rho1[nz_ghost][j][k]/bg->rho0[nz_ghost];
                    fg->T1[nz_ghost][j][k] = bg->T0[nz_ghost]*((GAMMA-1)/GAMMA * fg->p1[nz_ghost][j][k]/bg->p0[nz_ghost]);
                }
            }

            // Calculating ghost cells from neighbor below
            for (int i = 0; i < nz_ghost; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int k = 0; k < nx; k++)
                    {
                        fg->rho1[i][j][k] = gaussian_3D((i-nz_ghost)*dz, j*dy, k*dx, 0.5*dz*nz, 0.5*dy*ny, 0.5*dx*nx, 0.1*dz*nz, 0.1*dy*ny, 0.1*dx*nx, 1.0e-5);
                        fg->p1[i][j][k] = GAMMA*bg->p0[i]*fg->rho1[i][j][k]/bg->rho0[i];
                        fg->T1[i][j][k] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i][j][k]/bg->p0[i]);
                    }
                }
            }
        }
    #endif // DIMENSIONS
}