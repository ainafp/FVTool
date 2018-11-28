% Diffusion equation $\nabla \cdot (-D \nabla \phi) = \gamma$  
% for Phantom 2D small example or ones small example
% use superposition property and save potentials to construct
% connectivity matrix

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
figure(101); imagesc(im, [0,1]); colorbar; title('Image')

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

p0 = rowx_index(round(mnx/2)); % choose a sink

potentials = zeros(Nx, Ny, mny);
i = 1;
for p2=rowy_index'
    if p0~=p2
        RHSbc0 = RHSbc;
        RHSbc0(p0) = 1;                       % define current i
        RHSbc0(p2) = -1;                      % define current j
        c = solvePDE(meshstruct, M, RHSbc0);  % solve for the central scheme
        potentials(:,:,i) = c.value(2:end-1, 2:end-1);
        if plot_potentials
            figure(102); image(c.value(2:end-1, 2:end-1),'CDataMapping','scaled'); 
            colorbar; title('Potential'); pause(0.2);
        end
    end
    i = i + 1;
end

if plot_potentials
    figure(103); title('All potentials'); 
    for i=1:mny
        subplot(Nx, Ny, i); imagesc(potentials(:,:,i)); colorbar; 
    end
end

conduct0 = zeros(mnx, mny);
for p1=1:mnx
    for p2=1:mny
        if p2>p1
            pot = (potentials(:,:,p2) - potentials(:,:,p1));
            [a1, b1] = ind2sub(size(pot), p1);
            [a2, b2] = ind2sub(size(pot), p2);
            conduct0(p1, p2) = 1 / (pot(a1, b1) - pot(a2, b2));
            conduct0(p2, p1) = 1 / (pot(a2, b2) - pot(a1, b1));
        end
    end
end

% Plot results
conduct0(isnan(conduct0)) = 0;
figure(104); imagesc(abs(conduct0)); colorbar; title('Conductance')


%%
% To compare results with pointbypoint case
comparison = abs(conduct0) - abs(conductance(rowx_index,rowy_index));
sum(comparison(:))
figure(105); imagesc(comparison); colorbar;

