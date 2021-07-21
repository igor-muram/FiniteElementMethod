# Finite Element Method

<ul>
  <li>Finite element method for two-dimensional thermal problem for parabolic equation in Cartesian coordinate system.</li>
  <li>The default implicit four-layer scheme is used for time approximation, but a two- or three-layer scheme can also be used.</li>
  <li>Hierarchical cubic functions on triangles are used as basis functions.</li>
  <li>Boundary conditions of the first and second types.</li>
  <li>Equation coefficients are constant within a finite element.</li>
  <li>The coefficients are set according to the materials of the regions specified in the input file with the grid.</li>
  <li>Matrix of the system of linear equations is formed in sparse row-column format.</li>
  <li>The conjugate gradient method (CGM), local optimum scheme (LOS) or bis-conjugate gradient method (BSG) can be used to solve SLAEs.</li>
</ul>
