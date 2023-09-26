#include "global_parameters.h"


FLOAT_P central_mixed_derivative_second_order
(FLOAT_P up_right, FLOAT_P down_right, FLOAT_P up_left, FLOAT_P down_left, FLOAT_P dx_, FLOAT_P dy_)
{
    return 1.0/(4.0*dx_*dy_) * (up_right-up_left-down_right+down_left);
}