#ifndef FOREGROUND_DATA_1D_H_
#define FOREGROUND_DATA_1D_H_

#include "foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_1D.h"

void allocate_foreground_struct_1D(struct ForegroundVariables1D **fg, struct GridInfo1D *grid_info);
void deallocate_foreground_struct_1D(struct ForegroundVariables1D *fg);

#endif // FOREGROUND_DATA_1D_H_