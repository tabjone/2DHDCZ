#ifndef IO_FUNCTIONS_H__
#define IO_FUNCTIONS_H__

#include "hdf5.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "shared_files.h"
#include "global_parameters.h"

void save_foreground(struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int snap_number, FLOAT_P time);
void save_background(struct BackgroundVariables *bg, struct GridInfo *grid_info);

hid_t create_write_dataset(hid_t group, const char* name, hid_t datatype, hid_t dataspace, void* data, const char* unit);

void load_grid_info(struct GridInfo **grid_info, const char *file_path);
void load_background(struct BackgroundVariables *bg, struct GridInfo *grid_info, const char *file_path);
FLOAT_P load_foreground(struct ForegroundVariables2D *fg, struct GridInfo *grid_info, const char *file_path);

void read_dataset_line(hid_t dataset, const char* name, hid_t datatype, FLOAT_P *data);

void add_string_attribute(hid_t location, const char* attr_name, const char* attr_value);

#endif // IO_FUNCTIONS_H__