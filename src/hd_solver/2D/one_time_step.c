#include "hdf5.h"
#include "../../shared_files/shared_files.h"
#include "./hd_2D_functions.h"

void one_time_step_2D_hd(double **s1, double **p1, double **rho1, double **T1, double **vx, double **vz, double *grad_s0, double *p0, double *rho0, double *T0, double *g, int nz, int nx, double dx, double dz)
{

    // Should not have ghost cells but one-sided derivatives. Maybe I can write it so that if I put ghost cells it will work as one-sided derivatives.

    double rhs_dvx_dt, rhs_dvz_dt, rhs_ds1_dt;

    // Solve inside the grid
    for (int i=2; i < nz-2; i++)
    {
        for (int j=0; j < nx; j++)
        {
            rhs_dvx_dt = rhs_dvx_dt_2D_hd(p1, vx, vz, rho0, i, j, dx, dz, nx);
            rhs_dvz_dt = rhs_dvz_dt_2D_hd(rho1, p1, vx, vz, rho0, g, i, j, dx, dz, nx);
            rhs_ds1_dt = rhs_ds1_dt_2D_hd(s1, grad_s0, vx, vz, T0, rho0, i, j, dx, dz, nx);
        }
    }

    // Hard code periodic boundary conditions?

    for (int j = 1; j < nx-1; j++)
    {
        // Bottom boundary
        int i = 0;

        // Top boundary
        i = nz-1;
    }
}