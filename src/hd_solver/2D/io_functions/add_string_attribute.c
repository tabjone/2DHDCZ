#include "io_functions.h"

void add_string_attribute(hid_t location, const char* attr_name, const char* attr_value) 
{
    hid_t attr_space, attr_type, attr;
    hsize_t dims[1] = {1};

    // Create a datatype for the attribute (string)
    attr_type = H5Tcopy(H5T_C_S1);

    H5Tset_size(attr_type, strlen(attr_value) + 1);
    H5Tset_strpad(attr_type,H5T_STR_NULLTERM);

    // Create a dataspace for the attribute (scalar)
    attr_space = H5Screate_simple(1, dims, NULL);

    // Create the attribute on the given location (could be dataset or group)
    attr = H5Acreate2(location, attr_name, attr_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);

    // Write the string value to the attribute
    H5Awrite(attr, attr_type, attr_value);

    // Cleanup
    H5Aclose(attr);
    H5Sclose(attr_space);
    H5Tclose(attr_type);
}