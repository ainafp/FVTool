% Diffusion equation  $\nabla \cdot (-D \nabla \phi) = \gamma$  
% for a mesh of ones. Compute conductance between 2 points p1 and p2

% Input image
im = ones(100);

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
[Mbc, RHSbc] = boundaryCondition(BC); % boundary condition discretization
M = Mdiff + Mbc; % matrix of coefficient for central scheme

% Define currents (right hand side equation)
p1 = sub2ind([Ny + 2, Nx + 2], 72, 9);
p2 = sub2ind([Ny + 2, Nx + 2], 30, 80);
RHSbc(p1) = 1;
RHSbc(p2) = -1;

% Solve PDE
c = solvePDE(meshstruct, M, RHSbc); % solve for the central scheme

% Compute conductance
conductance = 1 / (c.value(p1) - c.value(p2));

% Plot results
figure; image(c.value,'CDataMapping','scaled'); colorbar;

