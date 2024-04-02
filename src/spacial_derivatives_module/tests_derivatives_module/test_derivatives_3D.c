#include <stdio.h>
#include "math.h"
#include "global_float_precision.h"

#include "../derivatives_3D/derivatives_3D.h"
#include "../../array_utilities/array_memory_management/array_memory_management.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#define Lx (FLOAT_P)1.0e1
#define nx 400
#define Ly (FLOAT_P)1.0e1
#define ny 400
#define Lz (FLOAT_P)1.0e1
#define nz 400

#define epsilon (FLOAT_P)1.0e-2


int test_central_first_derivative_x_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_2dx = 1.0 / (2.0 * dx);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = 3.0*M_PI/Lx*cos(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value = central_first_derivative_x_3D(f, i, j, k, nx, one_over_2dx);
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_first_derivative_y_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_2dy = 1.0 / (2.0 * dy);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = 4.0*M_PI/Ly*sin(3.0*M_PI*x/Lx) * cos(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value = central_first_derivative_y_3D(f, i, j, k, ny, one_over_2dy);
                
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_first_derivative_z_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_2dz = 1.0 / (2.0 * dz);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = 5.0*M_PI/Lz*sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * cos(5.0*M_PI*z/Lz);
                computed_value = central_first_derivative_z_3D(f, i, j, k, one_over_2dz);
                
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_second_derivative_x_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_dxdx = 1.0 / (dx * dx);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = -9.0*M_PI*M_PI/(Lx*Lx)*sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value = central_second_derivative_x_3D(f, i, j, k, nx, one_over_dxdx);
                
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_second_derivative_y_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_dydy = 1.0 / (dy * dy);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = -16.0*M_PI*M_PI/(Ly*Ly)*sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value = central_second_derivative_y_3D(f, i, j, k, ny, one_over_dydy);
                
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_second_derivative_z_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_dzdz = 1.0 / (dz * dz);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = -25.0*M_PI*M_PI/(Lz*Lz)*sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value = central_second_derivative_z_3D(f, i, j, k, one_over_dzdz);
                
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_second_derivative_xy_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_4dxdy = 1.0 / (4.0 * dx * dy);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = 12.0*M_PI*M_PI/(Lx*Ly)*cos(3.0*M_PI*x/Lx) * cos(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value = central_second_derivative_xy_3D(f, i, j, k, nx, ny, one_over_4dxdy);
                
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_second_derivative_xz_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_4dxdz = 1.0 / (4.0 * dx * dz);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = 15.0*M_PI*M_PI/(Lx*Lz)*cos(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * cos(5.0*M_PI*z/Lz);
                computed_value = central_second_derivative_xz_3D(f, i, j, k, nx, one_over_4dxdz);
                
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_second_derivative_yz_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_4dydz = 1.0 / (4.0 * dy * dz);

    FLOAT_P x, y, z;
    for (int i = 1; i < nz-1; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 1; i < nz-1; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = 20.0*M_PI*M_PI/(Ly*Lz)*sin(3.0*M_PI*x/Lx) * cos(4.0*M_PI*y/Ly) * cos(5.0*M_PI*z/Lz);
                computed_value = central_second_derivative_yz_3D(f, i, j, k, ny, one_over_4dydz);
                
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_third_derivative_xxy_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;

    allocate_3D_array(&f, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_8dxdxdy = 1.0 / (8.0 * dx * dx * dy);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }

    FLOAT_P expected_value, computed_value;

    // Calculating derivative
    for (int i = 2; i < nz-2; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                expected_value = - 36.0*M_PI*M_PI*M_PI/(Lx*Lx*Ly) * sin(3.0*M_PI*x/Lx) * cos(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value = central_third_derivative_xxy_3D(f, i, j, k, nx, ny, one_over_8dxdxdy);
                
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }            
            }
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_third_derivative_xxz_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;
    allocate_3D_array(&f, nz, ny, nx);
    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_8dxdxdz = 1.0 / (8.0 * dx * dx * dz);
    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
            }
        }
    }
    FLOAT_P expected_value, computed_value;
    // Calculating derivative
    for (int i = 2; i < nz-2; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                expected_value = - 45.0*M_PI*M_PI*M_PI/(Lx*Lx*Lz) * sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * cos(5.0*M_PI*z/Lz);
                computed_value = central_third_derivative_xxz_3D(f, i, j, k, nx, one_over_8dxdxdz);
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }
            }       
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_third_derivative_xyy_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;
    allocate_3D_array(&f, nz, ny, nx);
    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_8dydydx = 1.0 / (8.0 * dy * dy * dx);
    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx)*sin(4.0*M_PI*y/Ly)*sin(5.0*M_PI*z/Lz);
            }
        }
    }
    FLOAT_P expected_value, computed_value;
    // Calculating derivative
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                expected_value = - 48.0*M_PI*M_PI*M_PI/(Ly*Ly*Lx) * cos(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value = central_third_derivative_xyy_3D(f, i, j, k, nx, ny, one_over_8dydydx);
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }
            }       
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_third_derivative_xyz_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;
    allocate_3D_array(&f, nz, ny, nx);
    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_8dxdydz = 1.0 / (8.0 * dx * dy * dz);
    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx)*sin(4.0*M_PI*y/Ly)*sin(5.0*M_PI*z/Lz);
            }
        }
    }
    FLOAT_P expected_value, computed_value;
    // Calculating derivative
    for (int i = 2; i < nz-2; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                expected_value = 60.0*M_PI*M_PI*M_PI/(Lx*Ly*Lz) * cos(3.0*M_PI*x/Lx) * cos(4.0*M_PI*y/Ly) * cos(5.0*M_PI*z/Lz);
                computed_value = central_third_derivative_xyz_3D(f, i, j, k, nx, ny, one_over_8dxdydz);
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }
            }       
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_third_derivative_xzz_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;
    allocate_3D_array(&f, nz, ny, nx);
    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_8dxdzdz = 1.0 / (8.0 * dx * dz * dz);
    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx)*sin(4.0*M_PI*y/Ly)*sin(5.0*M_PI*z/Lz);
            }
        }
    }
    FLOAT_P expected_value, computed_value;
    // Calculating derivative
    for (int i = 2; i < nz-2; i++) 
    {
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                expected_value = - 75.0*M_PI*M_PI*M_PI/(Lx*Lz*Lz) * sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * cos(5.0*M_PI*z/Lz);
                computed_value = central_third_derivative_xzz_3D(f, i, j, k, nx, one_over_8dxdzdz);
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }
            }       
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_third_derivative_yyz_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;
    allocate_3D_array(&f, nz, ny, nx);
    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_8dydydz = 1.0 / (8.0 * dy * dy * dz);
    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx)*sin(4.0*M_PI*y/Ly)*sin(5.0*M_PI*z/Lz);
            }
        }
    }
    FLOAT_P expected_value, computed_value;
    // Calculating derivative
    for (int i = 2; i < nz-2; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                expected_value = - 80.0*M_PI*M_PI*M_PI/(Ly*Ly*Lz) * sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * cos(5.0*M_PI*z/Lz);
                computed_value = central_third_derivative_yyz_3D(f, i, j, k, ny, one_over_8dydydz);
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }
            }       
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_central_third_derivative_yzz_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;
    allocate_3D_array(&f, nz, ny, nx);
    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_8dydzdz = 1.0 / (8.0 * dy * dz * dz);
    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx)*sin(4.0*M_PI*y/Ly)*sin(5.0*M_PI*z/Lz);
            }
        }
    }
    FLOAT_P expected_value, computed_value;
    // Calculating derivative
    for (int i = 2; i < nz-2; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                expected_value = - 100.0*M_PI*M_PI*M_PI/(Ly*Lz*Lz) * sin(3.0*M_PI*x/Lx) * cos(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value = central_third_derivative_yzz_3D(f, i, j, k, ny, one_over_8dydzdz);
                if (fabs(expected_value - computed_value)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    return 1;
                }
            }       
        }
    }
    deallocate_3D_array(f);
    return 0;
}

int test_upwind_first_derivative_x_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;
    FLOAT_P ***velocity_positive, ***velocity_negative;
    allocate_3D_array(&f, nz, ny, nx);
    allocate_3D_array(&velocity_positive, nz, ny, nx);
    allocate_3D_array(&velocity_negative, nz, ny, nx);
    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);
    FLOAT_P one_over_dx = 1.0 / dx;
    FLOAT_P one_over_2dx = 1.0 / (2.0 * dx);
    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx)*sin(4.0*M_PI*y/Ly)*sin(5.0*M_PI*z/Lz);
                velocity_positive[i][j][k] = 1.0;
                velocity_negative[i][j][k] = -1.0;
            }
        }
    }
    FLOAT_P expected_value;
    FLOAT_P computed_value_negative, computed_value_positive;
    // Calculating derivative
    for (int i = 0; i < nz; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                expected_value = 3.0*M_PI/Lx*cos(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value_negative = upwind_first_derivative_x_3D(f, velocity_negative, i, j, k, nx, one_over_dx, one_over_2dx);
                computed_value_positive = upwind_first_derivative_x_3D(f, velocity_positive, i, j, k, nx, one_over_dx, one_over_2dx);

                if (fabs(expected_value - computed_value_negative)/fabs(expected_value) > epsilon || fabs(expected_value - computed_value_positive)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    deallocate_3D_array(velocity_positive);
                    deallocate_3D_array(velocity_negative);
                    return 1;
                }
            }       
        }
    }
    deallocate_3D_array(f);
    deallocate_3D_array(velocity_positive);
    deallocate_3D_array(velocity_negative);
    return 0;
}

int test_upwind_first_derivative_y_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;
    FLOAT_P ***velocity_positive, ***velocity_negative;

    allocate_3D_array(&f, nz, ny, nx);
    allocate_3D_array(&velocity_positive, nz, ny, nx);
    allocate_3D_array(&velocity_negative, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_dy = 1.0 / dy;
    FLOAT_P one_over_2dy = 1.0 / (2.0 * dy);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx)*sin(4.0*M_PI*y/Ly)*sin(5.0*M_PI*z/Lz);
                velocity_positive[i][j][k] = 1.0;
                velocity_negative[i][j][k] = -1.0;
            }
        }
    }
    FLOAT_P expected_value;
    FLOAT_P computed_value_negative, computed_value_positive;
    // Calculating derivative
    for (int i = 0; i < nz; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                expected_value = 4.0*M_PI/Ly*sin(3.0*M_PI*x/Lx) * cos(4.0*M_PI*y/Ly) * sin(5.0*M_PI*z/Lz);
                computed_value_negative = upwind_first_derivative_y_3D(f, velocity_negative, i, j, k, ny, one_over_dy, one_over_2dy);
                computed_value_positive = upwind_first_derivative_y_3D(f, velocity_positive, i, j, k, ny, one_over_dy, one_over_2dy);

                if (fabs(expected_value - computed_value_negative)/fabs(expected_value) > epsilon || fabs(expected_value - computed_value_positive)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    deallocate_3D_array(velocity_positive);
                    deallocate_3D_array(velocity_negative);
                    return 1;
                }
            }
        }
    }
    deallocate_3D_array(f);
    deallocate_3D_array(velocity_positive);
    deallocate_3D_array(velocity_negative);
    return 0;
}

int test_upwind_first_derivative_z_3D()
{
    // Returns 0 if test passes, 1 if test fails
    FLOAT_P ***f;
    FLOAT_P ***velocity_positive, ***velocity_negative;

    allocate_3D_array(&f, nz, ny, nx);
    allocate_3D_array(&velocity_positive, nz, ny, nx);
    allocate_3D_array(&velocity_negative, nz, ny, nx);

    FLOAT_P dx = Lx / (nx-1);
    FLOAT_P dy = Ly / (ny-1);
    FLOAT_P dz = Lz / (nz-1);

    FLOAT_P one_over_dz = 1.0 / dz;
    FLOAT_P one_over_2dz = 1.0 / (2.0 * dz);

    FLOAT_P x, y, z;
    for (int i = 0; i < nz; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 0; k < nx; k++)
            {   
                x = k * dx;
                f[i][j][k] = sin(3.0*M_PI*x/Lx)*sin(4.0*M_PI*y/Ly)*sin(5.0*M_PI*z/Lz);
                velocity_positive[i][j][k] = 1.0;
                velocity_negative[i][j][k] = -1.0;
            }
        }
    }
    FLOAT_P expected_value;
    FLOAT_P computed_value_negative, computed_value_positive;
    // Calculating derivative
    for (int i = 2; i < nz-2; i++) 
    {   
        z = i * dz;
        for (int j = 0; j < ny; j++) 
        {   
            y = j * dy;
            for (int k = 1; k < nx-1; k++)
            {   
                x = k * dx;
                expected_value = 5.0*M_PI/Lz*sin(3.0*M_PI*x/Lx) * sin(4.0*M_PI*y/Ly) * cos(5.0*M_PI*z/Lz);
                computed_value_negative = upwind_first_derivative_z_3D(f, velocity_negative, i, j, k, one_over_dz, one_over_2dz);
                computed_value_positive = upwind_first_derivative_z_3D(f, velocity_positive, i, j, k, one_over_dz, one_over_2dz);

                if (fabs(expected_value - computed_value_negative)/fabs(expected_value) > epsilon || fabs(expected_value - computed_value_positive)/fabs(expected_value) > epsilon)
                {
                    deallocate_3D_array(f);
                    deallocate_3D_array(velocity_positive);
                    deallocate_3D_array(velocity_negative);
                    return 1;
                }
            }
        }
    }
    deallocate_3D_array(f);
    deallocate_3D_array(velocity_positive);
    deallocate_3D_array(velocity_negative);
    return 0;
}



int test_derivatives_3D()
{
    if(test_central_first_derivative_x_3D()){
        printf("FAIL: test_central_first_derivative_x_3D.\n");}
    else{
        printf("PASS: test_central_first_derivative_x_3D.\n");}
    if(test_central_first_derivative_y_3D()){
        printf("FAIL: test_central_first_derivative_y_3D.\n");}
    else{
        printf("PASS: test_central_first_derivative_y_3D.\n");}
    if(test_central_first_derivative_z_3D()){
        printf("FAIL: test_central_first_derivative_z_3D.\n");}
    else{
        printf("PASS: test_central_first_derivative_z_3D.\n");}
    if(test_central_second_derivative_x_3D()){
        printf("FAIL: test_central_second_derivative_x_3D.\n");}
    else{
        printf("PASS: test_central_second_derivative_x_3D.\n");}
    if(test_central_second_derivative_y_3D()){
        printf("FAIL: test_central_second_derivative_y_3D.\n");}
    else{
        printf("PASS: test_central_second_derivative_y_3D.\n");}
    if(test_central_second_derivative_z_3D()){
        printf("FAIL: test_central_second_derivative_z_3D.\n");}
    else{
        printf("PASS: test_central_second_derivative_z_3D.\n");}
    if(test_central_second_derivative_xy_3D()){
        printf("FAIL: test central_second_derivative_xy_3D.\n");}
    else{
        printf("PASS: test central_second_derivative_xy_3D.\n");}
    if(test_central_second_derivative_xz_3D()){
        printf("FAIL: test central_second_derivative_xz_3D.\n");}
    else{
        printf("PASS: test central_second_derivative_xz_3D.\n");}
    if(test_central_second_derivative_yz_3D()){
        printf("FAIL: test central_second_derivative_yz_3D.\n");}
    else{
        printf("PASS: test central_second_derivative_yz_3D.\n");}
    if(test_central_third_derivative_xxy_3D()){
        printf("FAIL: test central_third_derivative_xxy_3D.\n");}
    else{
        printf("PASS: test central_third_derivative_xxy_3D.\n");}
    if(test_central_third_derivative_xxz_3D()){
        printf("FAIL: test central_third_derivative_xxz_3D.\n");}
    else{
        printf("PASS: test central_third_derivative_xxz_3D.\n");}
    if(test_central_third_derivative_xyy_3D()){
        printf("FAIL: test central_third_derivative_xyy_3D.\n");}
    else{
        printf("PASS: test central_third_derivative_xyy_3D.\n");}
    if(test_central_third_derivative_xyz_3D()){
        printf("FAIL: test central_third_derivative_xyz_3D.\n");}
    else{
        printf("PASS: test central_third_derivative_xyz_3D.\n");}
    if(test_central_third_derivative_xzz_3D()){
        printf("FAIL: test central_third_derivative_xzz_3D.\n");}
    else{
        printf("PASS: test central_third_derivative_xzz_3D.\n");}
    if(test_central_third_derivative_yyz_3D()){
        printf("FAIL: test central_third_derivative_yyz_3D.\n");}
    else{
        printf("PASS: test central_third_derivative_yyz_3D.\n");}
    if(test_central_third_derivative_yzz_3D()){
        printf("FAIL: test central_third_derivative_yzz_3D.\n");}
    else{
        printf("PASS: test central_third_derivative_yzz_3D.\n");}
    if(test_upwind_first_derivative_x_3D()){
        printf("FAIL: test upwind_first_derivative_x_3D.\n");}
    else{
        printf("PASS: test upwind_first_derivative_x_3D.\n");}
    if(test_upwind_first_derivative_y_3D()){
        printf("FAIL: test upwind_first_derivative_y_3D.\n");}
    else{
        printf("PASS: test upwind_first_derivative_y_3D.\n");}
    if(test_upwind_first_derivative_z_3D()){
        printf("FAIL: test upwind_first_derivative_z_3D.\n");}
    else{
        printf("PASS: test upwind_first_derivative_z_3D.\n");}

    return 0;
}