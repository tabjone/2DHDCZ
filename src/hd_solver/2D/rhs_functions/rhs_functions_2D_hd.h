#ifndef RHS_FUNCTIONS_2D_HD_H
#define RHS_FUNCTIONS_2D_HD_H

#include "../../../shared_files/structs/structs.h"

double rhs_dvx_dt_2D_hd(struct BackgroundVariables *background_variables, struct ForegroundVariables2D *foreground_variables, int i, int j, double dx, double dz, int nx);

double rhs_dvz_dt_2D_hd(struct BackgroundVariables *background_variables, struct ForegroundVariables2D *foreground_variables, int i, int j, double dx, double dz, int nx);

double rhs_ds1_dt_2D_hd(struct BackgroundVariables *background_variables, struct ForegroundVariables2D *foreground_variables, int i, int j, double dx, double dz, int nx);

#endif // RHS_FUNCTIONS_2D_HD_H