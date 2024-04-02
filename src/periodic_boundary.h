#ifndef PERIODIC_BOUNDARY_H__
#define PERIODIC_BOUNDARY_H__

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

#endif // PERIODIC_BOUNDARY_H__