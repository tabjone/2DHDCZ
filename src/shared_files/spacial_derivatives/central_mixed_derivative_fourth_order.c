double central_mixed_derivative_fourth_order
(double up_right, double down_right, double up_left, double down_left, 
double up2_left, double up2_right, double down2_left, double down2_right, double dx_, double dy_)
{
    return 1.0/(12.0*dx_*dy_)
           * ((-up2_right+8.0*up_right-8.0*down_right+down2_right)-(-up2_left+8.0*up_left-8.0*down_left+down2_left))
}