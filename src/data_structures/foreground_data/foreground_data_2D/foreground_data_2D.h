#ifndef FOREGROUND_DATA_2D_H_
#define FOREGROUND_DATA_2D_H_

#include "foreground_variables_struct_2D.h"
#include "../grid_info/grid_info_2D/grid_info_2D.h"

void allocate_foreground_struct_2D(struct ForegroundVariables2D **fg, struct GridInfo2D *grid_info);
void deallocate_foreground_struct_2D(struct ForegroundVariables2D *fg);

#endif // FOREGROUND_DATA_2D_H_