Gauss seidel og vanlige greier må bruke en funksjon "send_ghost_cells_2D(**array, nz_ghost, nz_full, ny)" der et av arguementene er nz_ghost slik at man kan sende 1,2,3 osv ghost cells istedetfor en spesiell for gauss-seidel. Det gjør det bare vanskeligere å teste og vanskeligere å finne bugs


OG HUSK: VED PERIODIC BOUNDARY ER F[0] ===== F[-1]. DE ER DE SAMME!!
