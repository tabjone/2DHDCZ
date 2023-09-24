#ifndef ARRAY_COPY_H__
#define ARRAY_COPY_H__

void copy_1D_array(double *src, double *dest, int nz);
void copy_2D_array(double **src, double **dest, int i_start, int i_end, int j_start, int j_end);
void copy_3D_array(double ***src, double ***dest, int nz, int ny, int nx);

#endif // ARRAY_COPY_H__