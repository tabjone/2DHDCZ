#ifndef IO_FUNCTIONS_H__
#define IO_FUNCTIONS_H__

#include "hdf5.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "shared_files.h"

void save_foreground(struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int snap_number, double time);
void save_background(struct BackgroundVariables *bg, struct GridInfo *grid_info);

hid_t create_write_dataset(hid_t group, const char* name, hid_t datatype, hid_t dataspace, void* data);


#endif // IO_FUNCTIONS_H__