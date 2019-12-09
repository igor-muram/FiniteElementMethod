# FiniteElementMethod
The finite element method for a two-dimensional thermal problem for an elliptic equation in a Cartesian coordinate system.<br>
The basic functions are cubic Lagrangian type on triangles.<br>
Boundary conditions of the first and second types.<br>
The coefficients of the equation are constant within the finite element.<br>
The coefficients are set in accordance with the materials of the areas indicated in the input file with the grid.<br>
The matrix of a system of linear equations is generated in a sparse row-column format.<br>
Conjugate gradient method with LU-preconditioning is used to solve SLAE, because the SLAE matrix is symmetric.
