Use refeq instead of ref

Link sticking out of page. Move to margen på bunn.

Theory about the sun.

Derivation of why the sound speed is filtered out

??? in background initialization

Should write (z,x,t) as (t,x,z). Always keep t first and then (x,y,z)

Boundary conditions

Plot P_radiation/P_gas to show why we only need to include gas pressure.

Plot from Solar S GONG adiabatic gradient, temp gradient. To show why we need k>0 in a numerical model -> Numerical convective instability must be bigger than the real sun.

Husk å skrive ned hvorfor vi må ha superad_param = k > 0. Det er fordi i solar S så klarer ikke et numerisk program å produsere konveksjon med de små >0 greiene i solar S, dette er noe sola klarer i virkeligheten. Få GONG data av sondre for å vise dette i thesis.

AST5770: Convective turnover time 111 days from bottom of CZ.

Mach number !(essential for approximation to hold)! and superadiabaticity parameter $\Delta\nabla$ should be low.

Alfvén Mach number $M_A=v_{ar}/c_{sr}$ should be comparable to the Mach number for the approximation to be valid. I.e magnetic energy should not greatly exceed kinetic energy (Lantz)

<b> Sondre spørsmål </b>

Magnetic diffusivity, hvordan finner jeg ut hva jeg skal sette den til?

Finn ut om likningene er for CGS units. For jeg har jo variabel UNITS

Hvordan ser upwind 3 schemet ut? Hva med upwind 4?

Er dette solar S paper? Ja + nettside (https://www.science.org/doi/10.1126/science.272.5266.1286)

<b> Future improvements, my version, in order: </b>

Improved and more versatile boundary conditions. (Maybe also ghosts on horizontal border)

Add differential rotation, coriolis effect, maybe solid body rotation below tachocline. This would create shear force.

Add fully copressible solver for outer layer of CZ and a simple atmosphere. (maybe not)

Add entropy source at radiation layer.

Include magnetic field. (yes next)

Extend from 2D to 3D. (okay)

More realistic equation of state, include ionization of helium. (maybe not)

MPI (yes after B-field included, before 3D)

<b> Possible improvement of the output </b>

Tracer particles (LATER)