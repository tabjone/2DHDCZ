#ifndef BACKGROUND_DATA_H_
#define BACKGROUND_DATA_H_

#include "background_variables_struct.h"

void allocate_background_struct(struct BackgroundVariables **bg, int nz_full);
void deallocate_background_struct(struct BackgroundVariables *bg);

#endif // BACKGROUND_DATA_H_