#include "mpi_functions.h"

void communicate_background_ghost_above_below(struct BackgroundVariables *bg, struct MpiInfo *mpi_info, int nz_full, int nz, int nz_ghost)
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

    communicate_1D_ghost_above_below(bg->rho0, mpi_info, nz_full, nz, nz_ghost);
    communicate_1D_ghost_above_below(bg->p0, mpi_info, nz_full, nz, nz_ghost);
    communicate_1D_ghost_above_below(bg->g, mpi_info, nz_full, nz, nz_ghost);
    communicate_1D_ghost_above_below(bg->T0, mpi_info, nz_full, nz, nz_ghost);
    communicate_1D_ghost_above_below(bg->grad_s0, mpi_info, nz_full, nz, nz_ghost);
    communicate_1D_ghost_above_below(bg->one_over_rho0, mpi_info, nz_full, nz, nz_ghost);
    communicate_1D_ghost_above_below(bg->eta_over_four_pi_rho0_T0, mpi_info, nz_full, nz, nz_ghost);
    communicate_1D_ghost_above_below(bg->grad_g, mpi_info, nz_full, nz, nz_ghost);
    communicate_1D_ghost_above_below(bg->grad_rho0, mpi_info, nz_full, nz, nz_ghost);
    communicate_1D_ghost_above_below(bg->r, mpi_info, nz_full, nz, nz_ghost);

}