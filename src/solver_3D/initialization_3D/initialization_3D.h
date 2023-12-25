#ifndef INITIALIZATION_3D_H__
#define INITIALIZATION_3D_H__

#include "global_parameters.h"
#include "shared_files.h"
#include "../equations_3D/equations_3D.h"
#include "global_initialization.h"
#include "../boundary_3D/boundary_3D.h"

FLOAT_P gaussian_3D(FLOAT_P x, FLOAT_P y, FLOAT_P z, FLOAT_P x0, FLOAT_P y0, FLOAT_P z0, FLOAT_P sigma_x, FLOAT_P sigma_y, FLOAT_P sigma_z, FLOAT_P A);

void initialize_foreground_zeros_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info);

void initialize_foreground_3D_entropy_pertubation(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

#endif // INITIALIZATION_3D_H__