TODO CODE:

Check that g(z) is the same in C-program as in python. Also see if the way this is calculated is correct. Og regn heller ut M(r) og bruk denne til å regne ut g(r)

Unit test derivative functions.

Husk å endre delta z i rhs funksjoner. Den er jo ikke konstant.

TODO PDFS:

Redo pertubation theory to end in continuity equation from anelastic approximation. Write part on equation of state ect. Discretize equations. Upwind/downwind ect.



TODO notes in anelastic approximation

<b>Code stuff to remember</b>

Keep track of Mach number and make sure this is low. Also that the superadiabaticity parameter $\Delta\nabla$ is low (this is also related to the Mach number by eq (9) in Lantz).

Keep track of "Alfvén Mach number" $M_A=v_{ar}/c_{sr}$. This needs to be comparable to the Mach number for the approximation to be valid. I.e magnetic energy should not greatly exceed kinetic energy (Lantz)

Change "current" in time derivatives to something else. Also this probably has to be changed entirely to account for RK methods.

Fix potential out of bounds problem in interpolate_solar_s.c
