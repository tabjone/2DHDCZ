#ifndef ITERATIVE_SOLVER_2D_H_
#define ITERATIVE_SOLVER_2D_H_

#include "global_float_precision.h"
#include "../../MPI_module/mpi_info_struct.h"

void jacobi_2D(FLOAT_P **rhs, FLOAT_P **current_solution, FLOAT_P **previous_solution, int nz, int nz_ghost, int ny, FLOAT_P dz, FLOAT_P dy, struct MpiInfo *mpi_info);
void gauss_seidel_2D_(FLOAT_P **rhs, FLOAT_P **current_solution, FLOAT_P **previous_solution, int nz, int nz_ghost, int ny, FLOAT_P dz, FLOAT_P dy, struct MpiIinfo *mpi_info);
void iterative_solver_2D(FLOAT_P **rhs, FLOAT_P **final_solution, FLOAT_P **initial_guess, FLOAT_P, int nz, int nz_ghost, int ny, FLOAT_P dz, FLOAT_P dy, struct MpiInfo *mpi_info);

#endif // ITERATIVE_SOLVER_2D_H_