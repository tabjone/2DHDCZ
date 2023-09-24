#ifndef SHARED_FILES_SHARED_FILES_H_
#define SHARED_FILES_SHARED_FILES_H_

int periodic_boundary(int i, int limit);

#include "hdf5.h"

#include "./structs/structs.h"
#include "./array_memory_management/array_memory_management.h"
#include "./spacial_derivatives/spacial_derivatives.h"
#include "./snapshot_io_operations/snapshot_io_operations.h"
#include "./solar_s_initialization/solar_s_initialization.h"
#include "./interpolation/interpolation.h"
#include "./extrapolation/extrapolation.h"
#include "./array_copy/array_copy.h"

#endif // SHARED_FILES_SHARED_FILES_H_