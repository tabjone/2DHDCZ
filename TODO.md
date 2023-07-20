<b>Code stuff to remember</b>

Keep track of Mach number and make sure this is low. Also that the superadiabaticity parameter $\Delta\nabla$ is low (this is also related to the Mach number by eq (9) in Lantz).

Create folder for solar s handling. Need interpolation.

Change "current" in time derivatives to something else. Also this probably has to be changed entirely to account for RK methods.

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
            tex
            pdf