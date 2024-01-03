#ifndef ITERATIVE_SOLVER_3D_H_
#define ITERATIVE_SOLVER_3D_H_

#include "global_float_precision.h"
#include "MPI_module/mpi_info_struct.h"

void iterative_solver_3D(FLOAT_P ***rhs, FLOAT_P ***final_solution, FLOAT_P ***initial_guess, int nz, int nz_ghost, int ny, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P dx, struct MpiInfo *mpi_info);
void jacobi_3D(FLOAT_P ***rhs, FLOAT_P ***current_solution, FLOAT_P ***previous_solution, int nz, int nz_ghost, int ny, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P dx, struct MpiInfo *mpi_info);
void gauss_seidel_3D(FLOAT_P ***rhs, FLOAT_P ***current_solution, FLOAT_P ***previous_solution, int nz, int nz_ghost, int ny, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P dx, struct MpiInfo *mpi_info);

#endif // ITERATIVE_SOLVER_3D_H_