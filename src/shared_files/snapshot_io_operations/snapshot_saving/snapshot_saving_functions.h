#ifndef SNAPSHOT_SAVING_SNAPSHOT_SAVING_FUNCTIONS_H
#define SNAPSHOT_SAVING_SNAPSHOT_SAVING_FUNCTIONS_H

#include "global_parameters.h"


void create_and_write_dataset_3D(hid_t file_id, const char* datasetname, hsize_t* dims, FLOAT_P*** data);
void create_and_write_dataset_2D(hid_t file_id, const char* datasetname, hsize_t* dims, FLOAT_P** data);
void create_and_write_dataset_1D(hid_t file_id, const char* datasetname, hsize_t* dims, FLOAT_P* data);

#endif // SNAPSHOT_SAVING_SNAPSHOT_SAVING_FUNCTIONS_H