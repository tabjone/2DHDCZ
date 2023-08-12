void calculate_pressure_scale_height(double *r_over_R, double *p0, double *H, int nz)
{
    /*
    Calculates the pressure scale height using definition right after eq. (7) in Lantz & Fan (1999)

    Parameters
    ----------
    r_over_R : double array
        Pointer to radius array normalized by the solar radius
    p0 : double array
        Pointer to pressure array of background state
    H : double array
        Pointer to pressure scale height array of background state
    nz : int
        Number of grid points in z-direction
    */

    // Radius of the sun in cgs units
    double R_sun = 6.957e10;

    // Handle the end points using forward and backward difference
    H[0] = - R_sun * (r_over_R[1] - r_over_R[0]) * p0[0]/ (p0[1] - p0[0]);
    H[nz-1] = - R_sun * (r_over_R[nz-1] - r_over_R[nz-2]) * p0[nz-1]/ (p0[nz-1] - p0[nz-2]);

    // Loop over the rest of the points using central difference
    for (int i = 1; i < nz-1; ++i) {
        H[i] = - R_sun * (r_over_R[i+1] - r_over_R[i-1]) * p0[i]/ (p0[i+1] - p0[i-1]);
    }
}