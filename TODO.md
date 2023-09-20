<b>IMPORTANT</b>:

Produce plots for initialization (solar S + superad_param from GONG) and von neumann

If code not running: Unit test derivative functions, test elliptic solver vs. analytical solution from here: https://aquaulb.github.io/book_solving_pde_mooc/solving_pde_mooc/notebooks/05_IterativeMethods/05_01_Iteration_and_2D.html#gauss-seidel-method

Redo pertubation theory to end in continuity equation from anelastic approximation. Write part on equation of state ect. Discretize equations. Upwind/downwind ect.

Re-write extrapolation functions to extrapolate_up, extrapolate_down sånn at det er bedre lagt opp for MPI senere

Husk å skrive ned hvorfor vi må ha superad_param = k > 0. Det er fordi i solar S så klarer ikke et numerisk program å produsere konveksjon med de små >0 greiene i solar S, dette er noe sola klarer i virkeligheten. Få GONG data av sondre for å vise dette i thesis.

<b> Sondre spørsmål </b>

Må vi sette boundary conditions for p1 solver i vertikal retning? Vi har jo periodic boundary

<b>THINGS TO REMEMBER</b>:

grad_s0 is discontinous. Maybe this can cause a problem.

Mach number !(essential for approximation to hold)! and superadiabaticity parameter $\Delta\nabla$ should be low.

Alfvén Mach number $M_A=v_{ar}/c_{sr}$ should be comparable to the Mach number for the approximation to be valid. I.e magnetic energy should not greatly exceed kinetic energy (Lantz)

Ideas for model expansion:

Regular MHD solver with simple physics and simple atmosphere for 0.97 solar radii and upwards for looking at convective overshooting.

Add solid rotating interiour under tachocline. No mass transfer ect, make some shear force at the bottom of CZ.