#ifndef SPACIAL_DERIVATIVES_SPACIAL_DERIVATIVE_FUNCTIONS_H_
#define SPACIAL_DERIVATIVES_SPACIAL_DERIVATIVE_FUNCTIONS_H_

double backward_first_derivative_first_order(double centre, double left, double dx_);
double backward_first_derivative_second_order(double centre, double left, double left2, double dx_);

double central_first_derivative_second_order(double left, double right, double dx_);
double central_second_derivative_second_order(double centre, double left, double right, double dx_);

double forward_first_derivative_first_order(double centre, double right, double dx_);
double forward_first_derivative_second_order(double centre, double right, double right2, double dx_);

#endif // SPACIAL_DERIVATIVES_SPACIAL_DERIVATIVE_FUNCTIONS_H_