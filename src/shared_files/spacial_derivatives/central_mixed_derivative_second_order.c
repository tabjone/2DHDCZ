double central_mixed_derivative_second_order
(double up_right, double down_right, double up_left, double down_left, double dx_, double dy_)
{
    return 1.0/(4.0*dx_*dy_) * (up_right-up_left-down_right+down_left);
}