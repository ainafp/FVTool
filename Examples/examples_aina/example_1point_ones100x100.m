% Diffusion equation $ \nabla. (-D \nabla \phi) = \gamma $  
% for D = ones(145)
% one only point

% Input image
im = phantom(2);
im = [ 1, 2, 1; 2, 3, 2; 3, 4, 3];
im = [ 1, 2; 2, 1];


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
p2 = sub2ind([Ny + 2, Nx + 2], 70, 138);
RHSbc(p1) = 1;
RHSbc(p2) = -1;

% Solve PDE
c = solvePDE(meshstruct, M, RHSbc); % solve for the central scheme

% Compute conductance
conductance = 1 / (c.value(p1) - c.value(p2));

% Plot results
figure; image(c.value,'CDataMapping','scaled'); colorbar;
savefig('example_0.fig')
saveas(gcf, 'example_0.eps', 'epsc')

