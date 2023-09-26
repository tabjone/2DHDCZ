<b> SONDRE </b>

Skriv likningene i et Latex dokument

Se over A4 ark.

Sett opp pertubasjon "For å debugge blir kanskje det enkle å sette en usymmetrisk perturbasjon (med sentrum rundt x=Lx/4 ) og null allerede ved grensen. F.eks. sett bredden på gaussen til w_x=Lx/10  eller mindre slik at den blir rundt 0 på grensen. evt. velg noe som blir eksakt null utenfor en grense."

<b> AFTER CODE WORKS </b>

GO OVER ENTIRE CODEBASE AND CLEAN UP. WRITE COMMENTS. WRITE PARAMETERS, RETURNS ECT. WRITE FUNCTIONS.

In solar S, load file in there. Remove old Io operations

Test signbit trick in separate file. 1-(signbit(vx))*upwind + downwind*signbit(vx) eller liknende med funksjoner som gjør  6 flops (ca. second order).
Dvs. double upwind(dx) og double downwind(dx) skal gjøre 6 FLOPS der 2 er divisions. Bare gjør en second order scheme likesågodt.

Optimization: #define PRECISION = 0 // 0 for float, 1 for double, 2 for long double. This needs to be implemented ASAP to hold codebase relevant for future.

<b> IMORGEN / Neste dagene </b>

do #define LOAD_SNAP = 0
#define LOAD_SNAP ? = snap_to_load 

funksjoner for eos og termo 1st law
random field initiation

<b>IMPORTANT</b>:

Produce plots for initialization (solar S + superad_param from GONG) and von neumann

If code not running: Unit test derivative functions, test elliptic solver vs. analytical solution from here: https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html#gauss-seidel-method

Redo pertubation theory to end in continuity equation from anelastic approximation. Write part on equation of state ect. Discretize equations. Upwind/downwind ect.

Re-write extrapolation functions to extrapolate_up, extrapolate_down sånn at det er bedre lagt opp for MPI senere

Husk å skrive ned hvorfor vi må ha superad_param = k > 0. Det er fordi i solar S så klarer ikke et numerisk program å produsere konveksjon med de små >0 greiene i solar S, dette er noe sola klarer i virkeligheten. Få GONG data av sondre for å vise dette i thesis.

<b> Sondre spørsmål </b>

Kan vi øke toleransen for gauss-seidel med det nye kriteriet?

<b>THINGS TO REMEMBER</b>:

AST5770: Convective turnover time 111 days from bottom of CZ.

grad_s0 is discontinous. Maybe this can cause a problem.

Mach number !(essential for approximation to hold)! and superadiabaticity parameter $\Delta\nabla$ should be low.

Alfvén Mach number $M_A=v_{ar}/c_{sr}$ should be comparable to the Mach number for the approximation to be valid. I.e magnetic energy should not greatly exceed kinetic energy (Lantz)

<b> From debugging </b>

Fra Periodic boundary funksjon:

Periodic boundary j=-4: Returning 6 and nx=10
Periodic boundary j=-3: Returning 7 and nx=10
Periodic boundary j=-2: Returning 8 and nx=10
Periodic boundary j=-1: Returning 9 and nx=10
Periodic boundary j=0: Returning 0 and nx=10
Periodic boundary j=1: Returning 1 and nx=10
Periodic boundary j=2: Returning 2 and nx=10
Periodic boundary j=3: Returning 3 and nx=10
Periodic boundary j=4: Returning 4 and nx=10
Periodic boundary j=5: Returning 5 and nx=10
Periodic boundary j=6: Returning 6 and nx=10
Periodic boundary j=7: Returning 7 and nx=10
Periodic boundary j=8: Returning 8 and nx=10
Periodic boundary j=9: Returning 9 and nx=10
Periodic boundary j=10: Returning 0 and nx=10
Periodic boundary j=11: Returning 1 and nx=10
Periodic boundary j=12: Returning 2 and nx=10
Periodic boundary j=13: Returning 3 and nx=10

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

<b> Future improvements, my version, in order: </b>

Implement MPI

Improved and more versatile boundary conditions.

Add differential rotation, coriolis effect, maybe solid body rotation below tachocline. This would create shear force.

Add fully copressible solver for outer layer of CZ and a simple atmosphere.

More realistic equation of state, include ionization of helium.

Add entropy source at radiation layer.

Include magnetic field.

Extend from 2D to 3D.

<b> Possible improvement of the output </b>
Tracer particles. (what does this mean?)