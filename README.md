# Finite Element Method

* Finite element method for two-dimensional thermal problem for parabolic equation in Cartesian coordinate system.
* The default implicit four-layer scheme is used for time approximation, but a two- or three-layer scheme can also be used.
* Hierarchical cubic functions on triangles are used as basis functions.
* Boundary conditions of the first and second types.
* Equation coefficients are constant within a finite element.
* The coefficients are set according to the materials of the regions specified in the input file with the grid.
* Matrix of the system of linear equations is formed in sparse row-column format.
* The conjugate gradient method (CGM), local optimum scheme (LOS) or bis-conjugate gradient method (BSG) can be used to solve SLAEs.
