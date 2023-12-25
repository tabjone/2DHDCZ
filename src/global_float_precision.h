#ifndef GLOBAL_FLOAT_PRECISION_H_
#define GLOBAL_FLOAT_PRECISION_H_

#include <mpi.h>

#define FLOAT_PRECISION 1 // 0 for float, 1 for double, 2 for long double

#if FLOAT_PRECISION == 0
    #define FLOAT_P float
    #define H5_FLOAT_P H5T_NATIVE_FLOAT
    #define MPI_FLOAT_P MPI_FLOAT
#elif FLOAT_PRECISION == 1
    #define FLOAT_P double
    #define H5_FLOAT_P H5T_NATIVE_DOUBLE
    #define MPI_FLOAT_P MPI_DOUBLE
#elif FLOAT_PRECISION == 2
    #define FLOAT_P long double
    #define H5_FLOAT_P H5T_NATIVE_LDOUBLE
    #define MPI_FLOAT_P MPI_LONG_DOUBLE
#else
    #error "Invalid FLOAT_PRECISION type specified."
#endif // FLOAT_PRECISION
#endif // GLOBAL_FLOAT_PRECISION_H_