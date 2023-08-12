void calculate_entropy_gradient(double *p0, double *rho0, double *T0, double *Gamma_1, double *H, double *superadiabacicity_parameter, double *entropy_gradient, int nz)
{
    /*
    Calculates the entropy gradient using eq. (72) in Lantz & Fan (1999)

    Parameters
    ----------
    p0 : double array
        Pointer to pressure array of background state
    rho0 : double array
        Pointer to density array of background state
    T0 : double array
        Pointer to temperature array of background state
    Gamma_1 : double array
        Pointer to adiabatic index array of background state
    H : double array
        Pointer to pressure scale height array of background state
    superadiabacicity_parameter : double array
        Pointer to superadiabacicity parameter array of background state
    entropy_gradient : double array
        Pointer to entropy gradient array of background state
    nz : int
        Number of grid points in z-direction
    */

    // Spesific heat with constant pressure
    double c_p;
    // Ideal gas constant normalized by the average mass per particle
    double r_star;

    for (int i = 0; i < nz; i++)
    {
        r_star = p0[i] / rho0[i] / T0[i];
        c_p = r_star / (1 - 1/Gamma_1[i]);

        entropy_gradient[i] = - (c_p / H[i]) * superadiabacicity_parameter[i];
    }
}