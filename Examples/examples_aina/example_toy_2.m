L = 0.01; % a 1 cm domain
Nx = size(C1, 1); % number of cells
m = createMesh2D(Nx,Nx,L,L); % create the mesh
BC = createBC(m); % construct the BC structure (Neumann by default)
%%
% Now as you may remeber, we have to switch from Neumann to Dirichlet
% boundary conditions
BC.left.a(:) = 0; BC.left.b(:) = 1; BC.left.c(:) = 1; % Left boundary to Dirichlet
BC.right.a(:) = 0; BC.right.b(:) = 1; BC.left.c(:) = 1; % right boundary to Dirichlet
%%
% The next sep is to define the diffusivity coefficient. In this FVTool,
% the physical properties of the domain are defined for each cell, with the
% function createCellVariable
D = createCellVariable(m, 1e-5); % assign a constant value of 1e-5 to diffusivity value on each cell
%%
% However, the transfer coefficients must be known on the face of each cell.
% For this reason, we have a few averaging schemes implemented in the
% Utilities folder. For a 1D domain, we can use a harmonic mean scheme:
D_face = harmonicMean(D); % average diffusivity value on the cell faces
%%
% where M is the matrix of coefficient that is going to be calculated using
% this toolbox. The matrix of coefficient, _M_ has two parts. The diffusion
% equation and the boundary conditions. They are calculated by:
M_diff = diffusionTerm(D_face); % matrix of coefficients for diffusion term
[M_bc, RHS_bc] = boundaryCondition(BC); % matrix of coefficient and RHS vector for the boundary condition
%%
% A vector of right hand side values are always obtained during the
% discretization of the boundary conditions.
% Now that the PDE is discretized, we can solve it by a Matlab linear solver.
c = solvePDE(m, M_diff+M_bc, RHS_bc);
%%
% finally, the resut can be visualized:
visualizeCells(c);
