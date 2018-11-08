function [M, RHSbc] = computeDiffusionMatrix3D(im, voxel_size)

% Construct mesh structure
Nx = size(im, 1);        % number of cells
Lx = Nx * voxel_size(1); % domain length
Ny = size(im, 2);        % number of cells
Ly = Ny * voxel_size(2); % domain length
Nz = size(im, 3); % number of cells
Lz = Nz * voxel_size(3); % domain length
Tx = size(im, 4); % size tensor
meshstruct = createMesh3D(Nx, Ny, Nz, Lx, Ly, Lz); 

% Define the boundary condition:
BC = createBC(meshstruct); % all Neumann boundary condition structure
[Mbc, RHSbc] = boundaryCondition(BC); % boundary condition discretization

% Define the transfer coefficients:
[xvalue, yvalue, zvalue] = extendBoundaryTensor3D_2v(im);
D = FaceVariable(meshstruct, xvalue, yvalue, zvalue);
Dave = arithmeticMeanTensor(D); % convert a cell variable to face variable
Mdiff = diffusionTermTensor3D_d(Dave);
M = Mdiff + Mbc; % matrix of coefficient for central scheme

end
