void copy_3D_array(double ***src, double ***dest, int nz, int ny, int nx)
{
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                dest[i][j][k] = src[i][j][k];
            }
        }
    }
}