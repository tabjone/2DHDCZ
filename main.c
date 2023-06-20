#include "functions.h"

int main(int argc, char *argv[])
{
    double **radius;
    allocate_2D_array(&radius, 10, 10);
    deallocate_2D_array(radius);
}
