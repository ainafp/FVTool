% Diffusion equation $ \nabla. (-D \nabla \phi) = \gamma $  
% for Phantom 2D small example
% Checking conductance in 3 points
% Validating superposition property

% Input image
im = phantom(10);
im = im(2:end-1,3:end-2);

im = rand(3);
%im = im(2:end-1,3:end-2);

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

% Define currents (right hand side equation)
%p1 = sub2ind([Ny + 2, Nx + 2],  37, 74);
%p2 = sub2ind([Ny + 2, Nx + 2], 104, 83);
%p3 = sub2ind([Ny + 2, Nx + 2],  55, 65);

%p1 = sub2ind([Ny + 2, Nx + 2], 5, 6);
p1 = sub2ind([Ny + 2, Nx + 2], 5, 7);
p2 = sub2ind([Ny + 2, Nx + 2], 2, 6);
p3 = sub2ind([Ny + 2, Nx + 2], 3, 8);

RHSbc = RHSbc0;
RHSbc(p1) = 1;
RHSbc(p3) = 0;
RHSbc(p2) = -1;

% Solve PDE
c = solvePDE(meshstruct, M, RHSbc); % solve for the central scheme

% Compute conductance
potential1 = (c.value(p1) - c.value(p2))
potential2 = (c.value(p2) - c.value(p3))
potential3 = (c.value(p1) - c.value(p3))

potential1 + potential2

B0 = c.value
figure; image(B0)

%conductance1 = 1 / (c.value(p1) - c.value(p2))
%conductance = 1 / (c.value(p2) - c.value(p3))
%conductance = 1 / (c.value(p1) - c.value(p3))


RHSbc = RHSbc0;
RHSbc(p1) = 1;
RHSbc(p3) = -1;

% Solve PDE
c = solvePDE(meshstruct, M, RHSbc); % solve for the central scheme

% Compute conductance
conductance = (c.value(p1) - c.value(p2))
conductance = (c.value(p2) - c.value(p3))
conductance2 = (c.value(p1) - c.value(p3))

B1 = c.value
figure; image(B1)

RHSbc = RHSbc0;
RHSbc(p2) = 1;
RHSbc(p3) = -1;

% Solve PDE
c = solvePDE(meshstruct, M, RHSbc); % solve for the central scheme

% Compute conductance
conductance = (c.value(p1) - c.value(p2))
conductance3 = (c.value(p2) - c.value(p3))
conductance = (c.value(p1) - c.value(p3))

B2 = c.value
figure; image(B2)

(B0-B1+B2)

potential1 + conductance2
conductance3
%%
% Plot results
%figure; image(im,'CDataMapping','scaled'); colorbar;
%savefig('example_phantom.fig')
%saveas(gcf, 'example_phantom.eps', 'epsc')
figure; image(c.value, 'CDataMapping', 'scaled'); colorbar;
savefig('example_phantom_0.fig')
saveas(gcf, 'example_phantom_0.eps', 'epsc')


 1/(1/conductance1 + 1/conductance2)
 conductance3
 
