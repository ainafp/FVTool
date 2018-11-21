% Diffusion equation $ \nabla. (-D \nabla \phi) = \gamma $  
% for Phantom example

eps = 0.000001;

% Input image
%im = phantom(10);
%im = im(2:end-1,3:end-2);
%im = zeros(3);
%im(1, :) = 1;
%im(:, 1) = 1;
im = zeros(3,3,2,2) + 0.5;

%im(1, :) = 0.8;
%im(:, 1) = 0.8;
%im(:,:,2) = 1 - im(:,:,1);
%im(im==1) = 0;

%im = im + eps;

% Construct mesh structure
Lx = size(im, 1); % domain length
Nx = size(im, 1); % number of cells
Ly = size(im, 2); % domain length
Ny = size(im, 2); % number of cells
meshstruct = createMesh2D(Nx, Ny, Lx, Ly); % construct mesh
x = meshstruct.cellcenters.x; % extract the cell center positions

% Define the boundary condition:
BC = createBC(meshstruct); % all Neumann boundary condition structure
[Mbc, RHSbc0] = boundaryCondition(BC); % boundary condition discretization

% Define the transfer coefficients:
D = createFaceVariable(meshstruct, 0.0);
D.xvalue = im(:,:,:,1);
D.xvalue(size(D.xvalue,1)+1,:,:) = D.xvalue(end,:,:);
D.yvalue = im(:,:,:,2);
D.yvalue(:,size(D.yvalue,1)+1,:) = D.yvalue(:,end,:);
%Mdiff = diffusionTerm(D);
[Mdiff, Mx, My] = diffusionTermTensor2D(D);
M = Mdiff + Mbc; % matrix of coefficient for central scheme

% Define currents (right hand side equation)
%p1 = sub2ind([Ny + 2, Nx + 2],  37, 74);
%p2 = sub2ind([Ny + 2, Nx + 2], 104, 83);
%p3 = sub2ind([Ny + 2, Nx + 2],  55, 65);

%p1 = sub2ind([Ny + 2, Nx + 2], 5, 6);
p1 = sub2ind([Ny + 2, Nx + 2], 3, 3);
p2 = sub2ind([Ny + 2, Nx + 2], 1, 2);
p3 = sub2ind([Ny + 2, Nx + 2], 3, 1);

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

