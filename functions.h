void allocate_2D_array(double ***array_ptr, int ny, int nx);
void deallocate_2D_array(double **array_ptr);

void create_and_write_dataset_3D(hid_t file_id, const char* datasetname, hsize_t* dims, double*** data);
void create_and_write_dataset_2D(hid_t file_id, const char* datasetname, hsize_t* dims, double** data);
void create_and_write_dataset_1D(hid_t file_id, const char* datasetname, hsize_t* dims, double* data);