<b>IMPORTANT</b>:

Implement boundary conditions

Finish intitialization.tex

Produce plots for initialization.tex

Use function for extrapolation in solar_s_background_initialization.c

Unit test derivative functions if we have bugs

Redo pertubation theory to end in continuity equation from anelastic approximation. Write part on equation of state ect. Discretize equations. Upwind/downwind ect.

<b>THINGS TO REMEMBER</b>:

grad_s0 is discontinous. Maybe this can cause a problem.

Mach number !(essential for approximation to hold)! and superadiabaticity parameter $\Delta\nabla$ should be low.

Alfvén Mach number $M_A=v_{ar}/c_{sr}$ should be comparable to the Mach number for the approximation to be valid. I.e magnetic energy should not greatly exceed kinetic energy (Lantz)

Change "current" in time derivatives to something else. Also this probably has to be changed entirely to account for RK methods.