#include "global_parameters.h"
#include "io_functions.h"

void read_dataset_line(hid_t dataset, const char* name, hid_t datatype, FLOAT_P *data)
{
    {
        if (H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
            fprintf(stderr, "Failed to read dataset\n");
    }
}