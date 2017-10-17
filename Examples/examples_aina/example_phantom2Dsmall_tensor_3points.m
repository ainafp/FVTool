% Diffusion equation $ \nabla. (-D \nabla \phi) = \gamma $ 
% for Phantom example

% Input image

im = zeros(5,5,2,2) + 0.5;
im(:,:,1,1) = 0.7;
im(:,:,2,2) = 0.3;
im(3:5,3:5,1,1) = 0.5;
im(3:5,3:5,2,2) = 0.5;
im(5,5,:,:) = 0.5;
%eps = 0.000001;
%im = im + eps;

%im = rand(3,3,2,2);

% Construct mesh structure
Lx = size(im, 1); % domain length
Nx = size(im, 1); % number of cells
Ly = size(im, 2); % domain length
Ny = size(im, 2); % number of cells
Tx = size(im, 3); % size tensor

% Define the boundary condition:
meshstruct = createMesh2D(Nx, Ny, Lx, Ly); % construct mesh
BC = createBC(meshstruct); % all Neumann boundary condition structure
[Mbc, RHSbc0] = boundaryCondition(BC); % boundary condition discretization

%Mdiff = diffusionTerm(D);
[xvalue, yvalue] = extendBoundaryTensor2D(im);
D = FaceVariable(meshstruct, xvalue, yvalue, []);
[Mdiff, Mx, My] = diffusionTermTensor2D(D);
M = Mdiff + Mbc; % matrix of coefficient for central scheme

% Define currents (right hand side equation)
p1 = sub2ind([Ny+2, Nx+2], 3, 3);
p2 = sub2ind([Ny+2, Nx+2], 5, 5);
p3 = sub2ind([Ny+2, Nx+2], 4, 4);

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

potential1 + potential2 - potential3

