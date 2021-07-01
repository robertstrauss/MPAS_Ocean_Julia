# MPAS Ocean in Julia

This is a shallow water model of the ocean capable of running on irregular, non-rectilinear, meshes.


## TODO
 * analyze performance and work to make GPU kernels more performant
     * memory management / thread communication
 * do boundary conditions in normal velocity tendency calculation
     * currently it may error with non-periodic borders, cellsOnEdge is not valid for both sides of a boundary edge
 * add nonlinear terms to tendency methods
 * restructure modules
     * divide mode_forward into a file for tendency calculations, and for integrators
     * modify MPAS_Ocean struct to have more useful fields, put prognostic vars, mesh data, in sub structure to keep seperate from meta information about the simulation
 * add more integrators