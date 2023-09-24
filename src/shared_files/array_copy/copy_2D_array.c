void copy_2D_array(double **src, double **dest, int i_start, int i_end, int j_start, int j_end)
{
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            dest[i][j] = src[i][j];
        }
    }
}