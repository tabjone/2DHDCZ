TODO notes in anelastic approximation
Create folder for interpolation
Create folder for exterpolation
Create calculate_background_entropy in shared solar_s folder
Remove entropy calculations from main_2D_hd


<b>Code stuff to remember</b>

Keep track of Mach number and make sure this is low. Also that the superadiabaticity parameter $\Delta\nabla$ is low (this is also related to the Mach number by eq (9) in Lantz).

Keep track of "Alfvén Mach number" $M_A=v_{ar}/c_{sr}$. This needs to be comparable to the Mach number for the approximation to be valid. I.e magnetic energy should not greatly exceed kinetic energy (Lantz)

Change "current" in time derivatives to something else. Also this probably has to be changed entirely to account for RK methods.

Fix potential out of bounds problem in interpolate_solar_s.c

Structure from main repo:

src

    main.c

    functions.h

    hd_solver

        ...

    mhd_solver

        ...

    shared_files

        array_memory_management

            array_allocation/...

            array_deallocation/...

        derivatives

            spacial_derivatives

            time_derivatives (TBD)

        snapshot_io_operations

            snapshot_saving/...

            snapshot_loading/...

        solar_s_initialization

            ... loading from solar_s, interpolating to spesific z-values

    public

        ... HTML,JS, Server (If time allows)

    data

        ... Folder with run_name/snap_names.h5

    additional_files

        ... Ex. .py File for creating .h5 from solar s downloaded data

    visualisation

        ... Python files for visualising/saving movies, (If time, link to visualise tab on webpage)

    thesis

        drafts
            ... Keep everything on spesific topic in same folder