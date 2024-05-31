#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/MPI_module.h"
#include "solver/equation_of_state/equation_of_state_2D/equation_of_state_2D.h"
#include "global_initialization.h"
#include "initialize_foreground_2D.h"
#include "global_constants.h"
#include "global_parameters.h"

void initialize_foreground_sod_shock_horizontal_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
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
    FLOAT_P dy = grid_info->dy;

    // We let 1/6 of the grid from the middle to the left be in the "up" state, 1/6 from the middle to the right be in the "up" state, and the rest in the "down" state. Meaning 2/6 = 1/3 of the domain is up and 2/3 is down

    FLOAT_P rho_up, p_up, rho_down, p_down;

    if (mpi_info->rank == 0)
    {
        printf("Initializing Sod Shock Tube\n");
        // For now just letting the left side be one of the top values from process 1
        rho_up = bg->rho0[nz];
        p_up = bg->p0[nz];
        rho_down = -rho_up;
        p_down = -p_up;
    }
    MPI_Bcast(&rho_up, 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p_up, 1, MPI_FLOAT_P, 0, MPI_COMM_WORLD);

    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);

    initialize_foreground_zeros_2D(fg, grid_info);

    FLOAT_P y;
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            y = j * dy;
            if (y > 0.33 * ny * dy && y < 0.66 * ny * dy)
            {
                fg->p1[i][j] = p_up;
                fg->rho1[i][j] = rho_up;
            }
            else
            {
                fg->p1[i][j] = p_down;
                fg->rho1[i][j] = rho_down;
            }
            // Getting T1 from equation of state
            fg->T1[i][j] = bg->T0[i] * (fg->p1[i][j]/bg->p0[i] - fg->rho1[i][j]/bg->rho0[i]);
            // Getting entropy from first law of thermodynamics
            fg->s1[i][j] = c_p * (fg->T1[i][j]/bg->T0[i] - fg->p1[i][j]/bg->p0[i]);
        }
    }

    communicate_2D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->rho1, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->T1, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->vy, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->vz, mpi_info, nz, nz_ghost, ny);
}