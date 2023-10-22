DET MÅ VÆRE ET PROBLEM I DE DERIVERTE

Hovedting nå: 

Remember to change viscosity to some real value

Gauss seidel må ta hensyn til periodic boundary i hver iterasjon

Må huske at i Gauss-Seidel så må den sette boundary på prosessen over

BUG: KODEN FÅR DT=0.0 ANNEN HVER KJØRING. WTF?


Not parallelized yet: Gauss-Seidel, Main, Solar S initialization, Saving/Loading

Move main_2D thing to main.c

Lot's of junk to remove in extrapolation folder

OPS OPS OPS! Jeg må tenke på offset i MPI. Det må være offset i initialization med pertubation (done). Det gjelder og for Solar S initialization. Kan endre z0 og z1 osv i bg struct. Det må gjøres endringer i pertubasjonsinitialiseringen for å ta hensyn til offset og det må gjøres endringer i solar s initialiseringen. Resten burde gå automatisk?

Med B-felt kan det hende vi må ha hastighet i 3 retninger for både x og og z. Se over likningene senere med d/dy=d/dx=0 for 1D. d/dx=0 for 2D.

Må bruke extrpolate_1D_array som velger constant eller hva enn den vil

Legg opp for implementering av B-felt

Se over RHS funksjoner etter fortegnsfeil og småbugs

Har antatt at rhs += -a - b gir samme resultat som rhs -= a + b

Skriv om: funksjonene i one_time_step, io_functions, extrapolation (nesten ferdig, bare litt gammelt igjen)
Extrapolation DONE, just need to remove old functions and get them out of the rk1 ect.

Put B-field term into rhs of vx,vy,vz and ellitpic equation rhs

Gjør overgangen til CZ mer smooth, aka. k-tallet kan være en lineær funksjon el. og se om det blir mer stabilt da.

Husk hakket i rho0

Prøv med litt annen initial, for eks pert i s1

Optimize Gauss-Seidel by pre-calculations

Det hadde vært nyttig å ha et dataset t som har [time, snapshot nr] slik at man kan se på snaps fra en viss tid uten å loade alt. Denne kan hete info.h5 og kanskje vi burde ha nx,nz osv her. (Dette kan vurderes, men her kunne og mpi_info vært)

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

Magnetic diffusivity, hvordan finner jeg ut hva jeg skal sette den til?

Finn ut om likningene er for CGS units. For jeg har jo variabel UNITS



Hva skal magnetic diffusivity settes til?

Hvordan ser upwind 3 schemet ut? Hva med upwind 4?

Er dette solar S paper? Ja + nettside (https://www.science.org/doi/10.1126/science.272.5266.1286)



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