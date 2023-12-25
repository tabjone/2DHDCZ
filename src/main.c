#include "functions.h"
#include "global_testing.h"

int main(int argc, char *argv[])
{
    #if TESTS == 1
        main_tests(argc, argv);
    #else
        main_run(argc, argv);
    #endif // TESTS
}