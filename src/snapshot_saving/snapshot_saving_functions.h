#ifndef SNAPSHOT_SAVING_SNAPSHOT_SAVING_FUNCTIONS_H
#define SNAPSHOT_SAVING_SNAPSHOT_SAVING_FUNCTIONS_H

void create_and_write_dataset_3D(hid_t file_id, const char* datasetname, hsize_t* dims, double*** data);
void create_and_write_dataset_2D(hid_t file_id, const char* datasetname, hsize_t* dims, double** data);
void create_and_write_dataset_1D(hid_t file_id, const char* datasetname, hsize_t* dims, double* data);

#endif // SNAPSHOT_SAVING_SNAPSHOT_SAVING_FUNCTIONS_H