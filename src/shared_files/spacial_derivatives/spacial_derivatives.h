#ifndef SPACIAL_DERIVATIVES_H__
#define SPACIAL_DERIVATIVES_H__

#include "global_parameters.h"

FLOAT_P backward_first_derivative_first_order(FLOAT_P centre, FLOAT_P left, FLOAT_P dx_);
FLOAT_P backward_first_derivative_second_order(FLOAT_P centre, FLOAT_P left, FLOAT_P left2, FLOAT_P dx_);

FLOAT_P central_first_derivative_second_order(FLOAT_P left, FLOAT_P right, FLOAT_P dx_);
FLOAT_P central_second_derivative_second_order(FLOAT_P centre, FLOAT_P left, FLOAT_P right, FLOAT_P dx_);

FLOAT_P forward_first_derivative_first_order(FLOAT_P centre, FLOAT_P right, FLOAT_P dx_);
FLOAT_P forward_first_derivative_second_order(FLOAT_P centre, FLOAT_P right, FLOAT_P right2, FLOAT_P dx_);

FLOAT_P central_mixed_derivative_second_order
(FLOAT_P up_right, FLOAT_P down_right, FLOAT_P up_left, FLOAT_P down_left, FLOAT_P dx_, FLOAT_P dy_);
FLOAT_P central_mixed_derivative_fourth_order
(FLOAT_P up_right, FLOAT_P down_right, FLOAT_P up_left, FLOAT_P down_left, 
FLOAT_P up2_left, FLOAT_P up2_right, FLOAT_P down2_left, FLOAT_P down2_right, FLOAT_P dx_, FLOAT_P dy_);

#endif // SPACIAL_DERIVATIVES_H__