void calculate_superadiabatic_parameter(double *p0, double *T0, double *superadiabacicity_parameter, double del_ad, int nz)
{
    /*
    Calculates the superadiabatic parameter for the solar model. The superadiabatic parameter is defined as:
    del - del_ad = dln(T)/dln(p)

    Parameters
    ----------
    p0 : double array
        Pointer to pressure array of background state
    T0 : double array
        Pointer to temperature array of background state
    superadiabacicity_parameter : double array
        Pointer to superadiabatic parameter array
    del_ad : double
        Adiabatic gradient
    nz : int
        Number of grid points
    */

    // Handle the end points using forward and backward difference
    superadiabacicity_parameter[0] = - del_ad + p0[0]/T0[0] * (T0[1] - T0[0]) / (p0[1] - p0[0]);
    superadiabacicity_parameter[nz-1] = - del_ad + p0[nz-1]/T0[nz-1] * (T0[nz-1] - T0[nz-2]) / (p0[nz-1] - p0[nz-2]);
    
    // Loop over the rest of the points using central difference
    for (int i = 1; i < nz-1; ++i) {
        superadiabacicity_parameter[i] = - del_ad + p0[i]/T0[i] * (T0[i+1] - T0[i-1]) / (p0[i+1] - p0[i-1]);
    }
}