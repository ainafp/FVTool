% Diffusion equation $\nabla \cdot (-D \nabla \phi) = \gamma$  
% for Phantom 2D small example or ones 3x3, computing connectivity matrix
% no superposition used

plot_potentials = false;
use_phantom = true;

% Input image
if use_phantom
    % input phantom in very low resolution
    im = phantom(10);
    im = im(2:end-1,3:end-2);
else
    % or a 3x3 image with some zeros and ones
    im = zeros(3);
    im(1, :) = 1;
    im(:, 1) = 1;
    im(3, 3) = 1;
end

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

% Define mask
Nx = Dave.domain.dims(1);
Ny = Dave.domain.dims(2);
G = reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
mnx = Nx*Ny;	mny = Nx*Ny;
rowx_index = reshape(G(2:Nx+1,2:Ny+1),mnx,1); % main diagonal x
rowy_index = reshape(G(2:Nx+1,2:Ny+1),mny,1); % main diagonal y

% Compute conductance
conductance = zeros(size(RHSbc, 1), size(RHSbc, 1));
for p1=rowx_index'
    for p2=rowy_index'
        if p2>p1
            RHSbc0 = RHSbc;
            RHSbc0(p1) = 1; % define current i
            RHSbc0(p2) = -1; % define current j
            c = solvePDE(meshstruct, M, RHSbc0); % solve for the central scheme
            conductance(p1, p2) = abs(1 / (c.value(p1) - c.value(p2)));
            conductance(p2, p1) = conductance(p1, p2);
            if ~isnan(conductance(p1, p2)) && ~isinf(conductance(p1, p2)) && plot_potentials
                figure(11); image(c.value(2:end-1,2:end-1), 'CDataMapping', 'scaled'); 
                colorbar; title('Potentials')
                pause(.2)
            end
        end
    end
end

% Plot results
conductance(isnan(conductance)) = 0;
figure(10); image(im,'CDataMapping','scaled'); colorbar; title('Image')
figure(12); image(conductance(rowx_index,rowy_index), 'CDataMapping', 'scaled'); 
colorbar; title('Conductance')