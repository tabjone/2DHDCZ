#ifndef SNAPSHOT_LOADING_SNAPSHOT_LOADING_FUNCTIONS_H
#define SNAPSHOT_LOADING_SNAPSHOT_LOADING_FUNCTIONS_H

void read_dataset_3D(hid_t file_id, const char* datasetname, hsize_t* dims, double*** data);
void read_dataset_2D(hid_t file_id, const char* datasetname, hsize_t* dims, double** data);
void read_dataset_1D(hid_t file_id, const char* datasetname, hsize_t* dims, double* data);

#endif // SNAPSHOT_LOADING_SNAPSHOT_LOADING_FUNCTIONS_H