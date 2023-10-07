#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include <stdlib.h>

void deallocate_foreground_struct(struct ForegroundVariables *fg)
{
    #if DIMENSIONS == 1
        deallocate_1D_array(fg->p1);
        deallocate_1D_array(fg->rho1);
        deallocate_1D_array(fg->T1);
        deallocate_1D_array(fg->s1);
        deallocate_1D_array(fg->vz);
        #if BFIELD_ON == 1
            deallocate_1D_array(fg->Bz);
        #endif // BFIELD_ON
    #elif DIMENSIONS == 2
        deallocate_2D_array(fg->p1);
        deallocate_2D_array(fg->rho1);
        deallocate_2D_array(fg->T1);
        deallocate_2D_array(fg->s1);
        deallocate_2D_array(fg->vz);
        deallocate_2D_array(fg->vx);
        #if BFIELD_ON == 1
            deallocate_2D_array(fg->Bz);
            deallocate_2D_array(fg->By);
        #endif // BFIELD_ON
    #elif DIMENSIONS == 3
        deallocate_3D_array(fg->p1);
        deallocate_3D_array(fg->rho1);
        deallocate_3D_array(fg->T1);
        deallocate_3D_array(fg->s1);
        deallocate_3D_array(fg->vz);
        deallocate_3D_array(fg->vx);
        deallocate_3D_array(fg->vy);
        #if BFIELD_ON == 1
            deallocate_3D_array(fg->Bz);
            deallocate_3D_array(fg->By);
            deallocate_3D_array(fg->Bx);
        #endif // BFIELD_ON
    #endif // DIMENSIONS

    free(fg);
}