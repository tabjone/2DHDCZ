#ifndef SNAPSHOT_LOADING_SNAPSHOT_LOADING_FUNCTIONS_H
#define SNAPSHOT_LOADING_SNAPSHOT_LOADING_FUNCTIONS_H

#include "global_parameters.h"

void read_dataset_3D(hid_t file_id, const char* datasetname, hsize_t* dims, FLOAT_P*** data);
void read_dataset_2D(hid_t file_id, const char* datasetname, hsize_t* dims, FLOAT_P** data);
void read_dataset_1D(hid_t file_id, const char* datasetname, hsize_t* dims, FLOAT_P* data);

#endif // SNAPSHOT_LOADING_SNAPSHOT_LOADING_FUNCTIONS_H