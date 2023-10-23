#ifndef SNAPSHOT_IO_OPERATIONS_H__
#define SNAPSHOT_IO_OPERATIONS_H__

#include <hdf5.h>
#include <string.h>

#include "shared_files.h"

hid_t create_write_dataset(hid_t group, const char* name, hid_t datatype, hid_t dataspace, void* data, const char* unit);
void add_string_attribute(hid_t location, const char* attr_name, const char* attr_value);

void save_background(struct BackgroundVariables *bg, struct MpiInfo *mpi_info, int nz_full, int nz, int nz_ghost, FLOAT_P dz, FLOAT_P z0, FLOAT_P z1);

void save_mpi_info(struct MpiInfo *mpi_info);

#endif // SNAPSHOT_IO_OPERATIONS_H__