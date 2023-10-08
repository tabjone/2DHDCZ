#ifndef SNAPSHOT_IO_OPERATIONS_H__
#define SNAPSHOT_IO_OPERATIONS_H__

#include <hdf5.h>
#include <string.h>

hid_t create_write_dataset(hid_t group, const char* name, hid_t datatype, hid_t dataspace, void* data, const char* unit);
void add_string_attribute(hid_t location, const char* attr_name, const char* attr_value);

#endif // SNAPSHOT_IO_OPERATIONS_H__