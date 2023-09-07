#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include "shared_files.h"

void extrapolate_2D(struct ForegroundVariables2D *fg);
void extrapolate_background(struct BackgroundVariables *bg);

#endif // EXTRAPOLATION_H