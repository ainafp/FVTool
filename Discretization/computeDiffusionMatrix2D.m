function [meshstruct, D, M, RHSbc] = computeDiffusionMatrix2D(im, voxel_size)
% Compute Diffusion Matrix M from mesh 2D of tensors
%
% Author: Aina Frau-Pascual


% Construct mesh structure
Nx = size(im, 1);        % number of cells
Lx = Nx * voxel_size(1); % domain length
Ny = size(im, 2);        % number of cells
Ly = Ny * voxel_size(2); % domain length
Tx = size(im, 3);        % size tensor
meshstruct = createMesh2D(Nx, Ny, Lx, Ly); 

% Define the boundary condition:
BC = createBC(meshstruct); % all Neumann boundary condition structure
[Mbc, RHSbc] = boundaryCondition(BC); % boundary condition discretization

% Define the transfer coefficients:
[xvalue, yvalue] = extendBoundaryTensor2D(im);
D = FaceVariable(meshstruct, xvalue, yvalue, []);
Dave = arithmeticMeanTensor(D); % convert a cell variable to face variable
Mdiff = diffusionTermTensor2D(Dave);
M = Mdiff + Mbc; % matrix of coefficient for central scheme

%figure(104); imagesc(Mbc); colormap jet; colorbar;

end
