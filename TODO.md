io_functions load_foreground all





Function that handles boundary: It takes in an array and extrapolates either extrapolates it or sends it to ghost cells. If periodic boundary we send it it to the ghost cells at the top/bottom ect. If not periodic we just extrapolate or send to ghost cells of mpi processes.

Function that applies damping ect: Takes the foreground struct and does whatever it should do on s1, vz.

Må skrive om periodic boundary greia. Den er dårlig nå.

Det fungerte bedre når jeg ikke regnet ut pressure, density, temperature imellom k1,k2,k3.

Add ghosts cells to all boundaries so we can have other than periodic. Then we don't need this j_minus, j_plus variable also.

Remember to change viscosity to some real value

Med B-felt kan det hende vi må ha hastighet i 3 retninger for både x og og z.

Husk hakket i rho0

Optimize Gauss-Seidel by pre-calculations

Det hadde vært nyttig å ha et dataset t som har [time, snapshot nr] slik at man kan se på snaps fra en viss tid uten å loade alt. Denne kan hete info.h5 og kanskje vi burde ha nx,nz osv her.

<b> TODO over time </b>

Test signbit trick in separate file. 1-(signbit(vx))*upwind + downwind*signbit(vx) eller liknende med funksjoner som gjør  6 flops (ca. second order). Dvs. double upwind(dx) og double downwind(dx) skal gjøre 6 FLOPS der 2 er divisions. Bare gjør en second order scheme likesågodt.

<b>IMPORTANT</b>:

If code not running: Unit test derivative functions, test elliptic solver vs. analytical solution from here: https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html#gauss-seidel-method