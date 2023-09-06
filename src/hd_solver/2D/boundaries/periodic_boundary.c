#include "boundaries.h"

int periodic_boundary(int i, int limit) 
{
    /*
    Calculates the new index for a given index under periodic boundary conditions.

    Parameters
    ----------
    i : int
        Index
    limit : int
        Size of array

    Returns
    -------
    int
        New index i given periodic boundary conditions

    Ex.
    ----
    i = -1, limit = 5 -> 4 (or 5?)
    i = 5, limit = 5 -> 0
    */
    return (i + limit) % (limit);
}