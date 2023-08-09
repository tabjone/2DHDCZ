#include "hdf5.h"
#include "../../shared_files/shared_files.h"
#include "./hd_2D_functions.h"

void one_time_step_2D_hd(int nz, int nx)
{
    //z-direction on first axis. x-direction on second axis
    double **s1, **vx, **vz, **p1, **rho1;
    double *s0, *T0, *rho0, *g;
    double dx, dz;

    allocate_1D_array(&s0, nz);
    allocate_1D_array(&T0, nz);
    allocate_1D_array(&rho0, nz);
    allocate_1D_array(&g, nz);

    allocate_2D_array(&s1, nz, nx);
    allocate_2D_array(&vx, nz, nx);
    allocate_2D_array(&vz, nz, nx);
    allocate_2D_array(&p1, nz, nx);
    allocate_2D_array(&rho1, nz, nx);

    for (int i = 0; i < nz; i++)
    {
        s0[i] = 1.0;
        T0[i] = 1.0;
        rho0[i] = 1.0;
        g[i] = 1.0;
        for (int j = 0; j < nx; j++)
        {
            s1[i][j] = 1.0;
            vx[i][j] = 1.0;
            vz[i][j] = 1.0;
            p1[i][j] = 1.0;
            rho1[i][j] = 1.0;
        }
    }

    double rhs_dvx_dt, rhs_dvz_dt, rhs_ds1_dt;

    // Solve inside the grid
    for (int i=1; i < nz-1; i++)
    {
        for (int j=1; j < nx-1; j++)
        {
            rhs_dvx_dt = rhs_dvx_dt_2D_hd(p1, vx, vz, rho0, i, j, dx, dz);
            rhs_dvz_dt = rhs_dvz_dt_2D_hd(rho1, p1, vx, vz, rho0, g, i, j, dx, dz);
            rhs_ds1_dt = rhs_ds1_dt_2D_hd(s1, s0, vx, vz, T0, rho0, i, j, dx, dz);
        }
    }

    // Solve on the boundaries

    deallocate_1D_array(s0);
    deallocate_1D_array(T0);
    deallocate_1D_array(rho0);
    deallocate_1D_array(g);

    deallocate_2D_array(s1);
    deallocate_2D_array(vx);
    deallocate_2D_array(vz);
    deallocate_2D_array(p1);
    deallocate_2D_array(rho1);
}