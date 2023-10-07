Hovedting nå: 
Se over RHS funksjoner etter fortegnsfeil og småbugs

Har antatt at rhs += -a - b gir samme resultat som rhs -= a + b

Skriv om: funksjonene i one_time_step, io_functions, initialization, equations

Ekstrapolar vy og vx slik at dvx/dz = 0 and dvy/dz = 0 (anti-symmetric ghost cells)
For nå er top/bottom boundary av vx = 0

Bytt fra Jacobi til Gauss-Seidel med en gang koden fungerer

Gjør overgangen til CZ mer smooth, aka. k-tallet kan være en lineær funksjon el. og se om det blir mer stabilt da.

Husk hakket i rho0

Test endre rekkefølge til elliptic -> ds1dt osv. -> T1, rho1

Prøv med litt annen initial, for eks pert i s1

Optimize Gauss-Seidel by pre-calculations

<b> TODO MPI </b>

Make sure dvx/dz is forced by the extrapolation, don't use constant extrapolation. We can do this by making sure the ghost cells are the negative of the "real" cells. Create a struct called mpi_info that keeps the number of processes and what process you are. Then we can use this to extrapolate if MPI_ON, but don't extrapolate one way or the other if MPI_OFF

<b> Latex </b>

Fiks link ut av siden. Fiks citep. Fiks ????, skriv inn likningsnavn

Skriv general, ideal anelastic. Elliptic.

Når ferdig: gå over og fiks.

Lag lite Latex dokument med alle likninger for møter.

Plot from Solar S GONG adiabatic gradient, temp gradient.

<b> Visualisation </b>

Sett vmin=-vmax

<b> TODO over time </b>

Test signbit trick in separate file. 1-(signbit(vx))*upwind + downwind*signbit(vx) eller liknende med funksjoner som gjør  6 flops (ca. second order). Dvs. double upwind(dx) og double downwind(dx) skal gjøre 6 FLOPS der 2 er divisions. Bare gjør en second order scheme likesågodt.

Plot radiative pressure / gas pressure for entire solar S to see if we need to add radiative pressure for the radiative zone if we want to go deeper.

use array copy function in gauss-seidel

<b> Info </b>

Why the anelastic model filters out soundwaves:

use dvrho/dt = ...
and something else = ...
(conservation of mass, conservation of momentum (I think))

Take the divergence of one and the time derivative of the other so you can insert these equations into eachother. Use the p/rho^gamma=const relation and you will end up with the wave equation. With the pertubations this gives zero or something.


<b> LIFE </b>

Long term: Fiks pult. 170x60 https://www.byggmax.no/benkeplate-furu-p09310

<b>IMPORTANT</b>:

If code not running: Unit test derivative functions, test elliptic solver vs. analytical solution from here: https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html#gauss-seidel-method

Redo pertubation theory to end in continuity equation from anelastic approximation. Write part on equation of state ect. Discretize equations. Upwind/downwind ect.

Re-write extrapolation functions to extrapolate_up, extrapolate_down sånn at det er bedre lagt opp for MPI senere

Husk å skrive ned hvorfor vi må ha superad_param = k > 0. Det er fordi i solar S så klarer ikke et numerisk program å produsere konveksjon med de små >0 greiene i solar S, dette er noe sola klarer i virkeligheten. Få GONG data av sondre for å vise dette i thesis.

AST5770: Convective turnover time 111 days from bottom of CZ.

grad_s0 is discontinous. Maybe this can cause a problem.

Mach number !(essential for approximation to hold)! and superadiabaticity parameter $\Delta\nabla$ should be low.

Alfvén Mach number $M_A=v_{ar}/c_{sr}$ should be comparable to the Mach number for the approximation to be valid. I.e magnetic energy should not greatly exceed kinetic energy (Lantz)

<b> Sondre spørsmål </b>

Hvordan ser upwind 3 schemet ut? Hva med upwind 4?

Er dette solar S paper? Ja + nettside

<b> Future improvements, my version, in order: </b>

Implement MPI

Improved and more versatile boundary conditions. (driver med dette, spunge)

Add differential rotation, coriolis effect, maybe solid body rotation below tachocline. This would create shear force.

Add fully copressible solver for outer layer of CZ and a simple atmosphere.

Add entropy source at radiation layer.

Include magnetic field.

Extend from 2D to 3D.

More realistic equation of state, include ionization of helium.

<b> Possible improvement of the output </b>

Tracer particles (LATER)