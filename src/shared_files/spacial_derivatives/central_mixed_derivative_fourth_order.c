#include "global_parameters.h"


FLOAT_P central_mixed_derivative_fourth_order
(FLOAT_P up_right, FLOAT_P down_right, FLOAT_P up_left, FLOAT_P down_left, 
FLOAT_P up2_left, FLOAT_P up2_right, FLOAT_P down2_left, FLOAT_P down2_right, FLOAT_P dx_, FLOAT_P dy_)
{
    return 1.0/(12.0*dx_*dy_)
           * ((-up2_right+8.0*up_right-8.0*down_right+down2_right)-(-up2_left+8.0*up_left-8.0*down_left+down2_left));
}