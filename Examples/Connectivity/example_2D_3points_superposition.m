% Phantom 2D small example. Compare potentials in 3 points. 
% Validating superposition property
% Author: Aina Frau-Pascual


% Input image
im = phantom(10);
im = im(2:end-1,3:end-2);

% Construct mesh structure
Lx = size(im, 1); % domain length
Nx = size(im, 1); % number of cells
Ly = size(im, 2); % domain length
Ny = size(im, 2); % number of cells
meshstruct = createMesh2D(Nx, Ny, Lx, Ly); % construct mesh
x = meshstruct.cellcenters.x; % extract the cell center positions

% Define the boundary condition:
BC = createBC(meshstruct); % all Neumann boundary condition structure

% Define the transfer coefficients:
D = createCellVariable(meshstruct, im);
Dave = harmonicMean(D); % convert a cell variable to face variable
Mdiff = diffusionTerm(Dave); % diffusion term
[Mbc, RHSbc0] = boundaryCondition(BC); % boundary condition discretization
M = Mdiff + Mbc; % matrix of coefficient for central scheme

% Define points that I will look at
p1 = sub2ind([Ny + 2, Nx + 2], 5, 7);
p2 = sub2ind([Ny + 2, Nx + 2], 2, 6);
p3 = sub2ind([Ny + 2, Nx + 2], 3, 8);

% Define currents (right hand side equation)
RHSbc = RHSbc0;
RHSbc(p1) = 1;
RHSbc(p2) = -1;

% Solve PDE
c = solvePDE(meshstruct, M, RHSbc); % solve for the central scheme

% Compute and compare potentials
potential1 = (c.value(p1) - c.value(p2));
potential2 = (c.value(p2) - c.value(p3));
potential3 = (c.value(p1) - c.value(p3));

% Compare that they are equal
potential1 + potential2
potential3
