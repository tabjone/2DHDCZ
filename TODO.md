Burde endre nx til ny slik at vi får mer lett leselig indeksering. Dvs j er alltid ny sin index

Se over RHS funksjonene etter bugs. Fant en bare tilfeldig

FINISHED: rhs_dvz_dt_1D, 2D, 3D

Check if this rhs -= vy[i][j]*upwind_ds1_dy + vz[i][j]*upwind_ds1_dz; is the same as
rhs += -vy[i][j]*upwind_ds1_dy - vz[i][j]*upwind_ds1_dz;

Har ikke gått gjennom rk funksjoner eller elliptic solver

NEED TO EXTRAPOLATE VY AND VX SO THAT dvx/dz = 0 and dvy/dz = 0 (this can be done by making the ghost points anti-symmetrical)
NEED TO EXTRAPOLATE s1 SO THAT ds1/dt = 0 at boundary. maybe I should force this instead

Bytt fra Jacobi til Gauss-Seidel med en gang koden fungerer

RHS Elliptic equation

For now I have set the bottom and top boundary of vx = 0

<b> NESTE MØTE </b>

Sondre går over Viscosity

Sondre går over Thermal diffusivity ( Vi har numerical diffusivity allerede så kanskje vi heller kan implementere et ledd som gjør numerical diffusivity -> Thermal diffusivity (Dette vet Juan mye om))

Vi må starte å gå over koden som er implementer. Neste møte kan vi gå over rk1,rk2,rk3 i 2D




Koden ser ut til å kjøre bra men får problemer med top boundary hvis max_dt er for høy. Ser ut til at hvor høyt i CZ vi setter top boundary er relatert til hvor stor max_dt vi kan ha. Ser iværtfall at selv med en liten dt for å få koden igang så velger den omtrent alltid max_dt uansett hvor høy den er pga at kravet vx/dx+vz/dz har så enormt stor dx,dz. Prøv å kjøre uten g(z) leddet. Det må vel være relatert til det?

Gjør overgangen til CZ mer smooth, aka. k-tallet kan være en lineær funksjon el. og se om det blir mer stabilt da.

Husk hakket i rho0

I ytterste while loop: ha en variabel dt_max som sendes videre inn i funksjonene

Noe liknende dette:
elif t<10000 dt_max = 1000
elif t<1000 dt_max = 100
elif t<100 dt_max = 10

Test endre rekkefølge til elliptic -> ds1dt osv. -> T1, rho1

Prøv med litt annen initial, for eks pert i s1

#define GRAVITY_ON 0, 1

<b> TODO MPI </b>
Make sure dvx/dz is forced by the extrapolation, don't use constant extrapolation. We can do this by making sure the ghost cells are the negative of the "real" cells. Create a struct called mpi_info that keeps the number of processes and what process you are. Then we can use this to extrapolate if MPI_ON, but don't extrapolate one way or the other if MPI_OFF

<b> TODO ASAP </b>

Create first Latex draft.

Completed:
- Initialization (but what is the one equation with ????? and why is \citep not working)
- Elliptic: Not finished, marked by TODO

Fix visualisation script. Diverging colormap i.e RdBu. Sett vmin=-vmax. Fix quivers

Test pertubation in small area, not random initialization.

Rk3 test it

Clean up code.

Write all equations into Latex document to have ready at meetings

<b> SIMULATIONS </b>

Small pertubation for rk1,2,3 with upwind 1,2



<b> TODO over time </b>

Optimize Gauss-Seidel (Almost ASAP)

GO OVER ENTIRE CODEBASE AND CLEAN UP. WRITE COMMENTS. WRITE PARAMETERS, RETURNS ECT. WRITE FUNCTIONS.

Test signbit trick in separate file. 1-(signbit(vx))*upwind + downwind*signbit(vx) eller liknende med funksjoner som gjør  6 flops (ca. second order).
Dvs. double upwind(dx) og double downwind(dx) skal gjøre 6 FLOPS der 2 er divisions. Bare gjør en second order scheme likesågodt.

In solar S, load file in there. Remove old Io operations

Plot radiative pressure / gas pressure for entire solar S to see if we need to add radiative pressure for the radiative zone if we want to go deeper.

Check rk2 inplementation. Do I pass the correct arguments when calculating k2's, dvs. jeg skal ha de gamle verdiene for s1,vx om jeg renger ut vz, men den nye kun for vz i den utregningen. (Dette kan kanskje gjøres ved å lage en struct som bare peker til de gamle verdiene for de andre, men til den nye for seg selv)

use array copy function in gauss-seidel

<b> Info </b>

Why the anelastic model filters out soundwaves:

use dvrho/dt = ...
and something else = ...
(conservation of mass, conservation of momentum (I think))

Take the divergence of one and the time derivative of the other so you can insert these equations into eachother. Use the p/rho^gamma=const relation and you will end up with the wave equation. With the pertubations this gives zero or something.


<b> LIFE </b>

Kaste batterier, kaste alluminium. Kjøp tynne plastsøpleposer.

Long term: Fiks pult. 170x60 https://www.byggmax.no/benkeplate-furu-p09310

<b>IMPORTANT</b>:

Produce plots for initialization (solar S + superad_param from GONG) and von neumann

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

<b> Improvement of gauss-seidel, problem </b>

Problem med periodic boundary i p1:

Orginal Gauss-Seidel:
for j in range(1, ny-1):
        for i in range(1, nx-1):
            pnew[i, j] = (0.25 * (p[i-1, j]+p[i+1, j]+p[i, j-1]
                       + p[i, j+1]-b[i, j]*dx**2))

Gauss-Seidel med 2x speedup:
for j in range(1, ny-1):
        for i in range(1, nx-1):
            pnew[i, j] = (0.25 * (pnew[i-1, j]+p[i+1, j]+pnew[i, j-1]
                       + p[i, j+1]-b[i, j]*dx**2))

2x speedup går ikke med periodic boundary fordi pnew[i, j-1] er ikke regnet ut når vi er i punkt [0][0]
Possible to find trick for the 2x speedup even with periodic boundary?

<b> Future improvements, my version, in order: </b>

Implement MPI

Improved and more versatile boundary conditions.

Add differential rotation, coriolis effect, maybe solid body rotation below tachocline. This would create shear force.

Add fully copressible solver for outer layer of CZ and a simple atmosphere.

Add entropy source at radiation layer.

Include magnetic field.

Extend from 2D to 3D.

More realistic equation of state, include ionization of helium.

<b> Possible improvement of the output </b>

Tracer particles (LATER)