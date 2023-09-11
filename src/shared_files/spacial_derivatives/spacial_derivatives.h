#ifndef SPACIAL_DERIVATIVES_H__
#define SPACIAL_DERIVATIVES_H__

double backward_first_derivative_first_order(double centre, double left, double dx_);
double backward_first_derivative_second_order(double centre, double left, double left2, double dx_);

double central_first_derivative_second_order(double left, double right, double dx_);
double central_second_derivative_second_order(double centre, double left, double right, double dx_);

double forward_first_derivative_first_order(double centre, double right, double dx_);
double forward_first_derivative_second_order(double centre, double right, double right2, double dx_);

double central_mixed_derivative_second_order
(double up_right, double down_right, double up_left, double down_left, double dx_, double dy_);
double central_mixed_derivative_fourth_order
(double up_right, double down_right, double up_left, double down_left, 
double up2_left, double up2_right, double down2_left, double down2_right, double dx_, double dy_);

#endif // SPACIAL_DERIVATIVES_H__