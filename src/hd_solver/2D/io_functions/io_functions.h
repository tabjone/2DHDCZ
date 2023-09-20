#ifndef IO_FUNCTIONS_H__
#define IO_FUNCTIONS_H__

#include "hdf5.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "shared_files.h"

void save_foreground(struct ForegroundVariables2D *fg, int snap_number, double time);
void save_background(struct BackgroundVariables *bg);

#endif // IO_FUNCTIONS_H__