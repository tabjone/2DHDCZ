Move derivatives into 2D part

Move structs into 2D part

Remember to change viscosity to some real value

Med B-felt kan det hende vi må ha hastighet i 3 retninger for både x og og z.

Se over RHS funksjoner etter fortegnsfeil og småbugs

Har antatt at rhs += -a - b gir samme resultat som rhs -= a + b

Husk hakket i rho0

Optimize Gauss-Seidel by pre-calculations

Det hadde vært nyttig å ha et dataset t som har [time, snapshot nr] slik at man kan se på snaps fra en viss tid uten å loade alt. Denne kan hete info.h5 og kanskje vi burde ha nx,nz osv her.

<b> TODO over time </b>

Test signbit trick in separate file. 1-(signbit(vx))*upwind + downwind*signbit(vx) eller liknende med funksjoner som gjør  6 flops (ca. second order). Dvs. double upwind(dx) og double downwind(dx) skal gjøre 6 FLOPS der 2 er divisions. Bare gjør en second order scheme likesågodt.

Plot radiative pressure / gas pressure for entire solar S to see if we need to add radiative pressure for the radiative zone if we want to go deeper.

<b>IMPORTANT</b>:

If code not running: Unit test derivative functions, test elliptic solver vs. analytical solution from here: https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html#gauss-seidel-method