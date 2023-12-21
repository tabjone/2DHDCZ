#ifndef FOREGROUND_DATA_3D_H_
#define FOREGROUND_DATA_3D_H_

#include "foreground_variables_struct_3D.h"
#include "../grid_info/grid_info_3D/grid_info_3D.h"

void allocate_foreground_struct_3D(struct ForegroundVariables3D **fg, struct GridInfo3D *grid_info);
void deallocate_foreground_struct_3D(struct ForegroundVariables3D *fg);

#endif // FOREGROUND_DATA_3D_H_