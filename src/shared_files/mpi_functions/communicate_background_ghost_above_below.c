#include "mpi_functions.h"

void communicate_background_ghost_above_below(struct BackgroundVariables *bg, struct GridInfo *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Communicates the background variables above and below the current process.

    Parameters
    ----------
    bg : struct BackgroundVariables
        Background variables
    grid_info : struct GridInfo
        Grid parameters
    mpi_info : struct MPIVariables
        MPI variables
    */

    // Getting grid info
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;

    // Getting MPI info
    int rank = mpi_info->rank;
    int size = mpi_info->size;

    communicate_1D_ghost_above_below(bg->rho0, grid_info, mpi_info);
    communicate_1D_ghost_above_below(bg->p0, grid_info, mpi_info);
    communicate_1D_ghost_above_below(bg->g, grid_info, mpi_info);
    communicate_1D_ghost_above_below(bg->T0, grid_info, mpi_info);
    communicate_1D_ghost_above_below(bg->grad_s0, grid_info, mpi_info);
    communicate_1D_ghost_above_below(bg->one_over_rho0, grid_info, mpi_info);
    communicate_1D_ghost_above_below(bg->eta_over_four_pi_rho0_T0, grid_info, mpi_info);
    communicate_1D_ghost_above_below(bg->grad_g, grid_info, mpi_info);
    communicate_1D_ghost_above_below(bg->grad_rho0, grid_info, mpi_info);
}