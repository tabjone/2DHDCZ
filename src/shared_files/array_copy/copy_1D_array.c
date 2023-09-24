void copy_1D_array(double *src, double *dest, int nz)
{
    for (int i = 0; i < nz; i++)
    {
        dest[i] = src[i];
    }
}