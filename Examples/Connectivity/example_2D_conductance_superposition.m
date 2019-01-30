% Phantom 2D small example or ones 3x3, compute conductance matrix
% use superposition property 
% Compare to point-by-point case generated with:
% example_2D_conductance_pointbypoint.m
% Author: Aina Frau-Pascual


% Set plot_potentials to true if you want to plot all potentials
% WARNING: it is time consuming
plot_potentials = false;   

% Change use_phantom flag to use a different input image
use_phantom = true;        

% Input image: choose among the two (default: phantom)
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
mnx = Nx*Ny;	
mny = Nx*Ny;
rowx_index = reshape(G(2:Nx+1,2:Ny+1),mnx,1); % main diagonal x
rowy_index = reshape(G(2:Nx+1,2:Ny+1),mny,1); % main diagonal y

% Compute conductance for all pairs of voxels source-sink: O(N)
% First you choose your sink p0
p0 = rowx_index(round(mnx/2)); 
% You will have to save all the potentials to be able to combine them
% It is a compromise time / memory used
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
            figure(101); image(c.value(2:end-1, 2:end-1),'CDataMapping','scaled'); 
            colorbar; title('Potentials'); xlabel('x'); ylabel('y');
            pause(0.2);
        end
    end
    i = i + 1;
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
conduct0(isnan(conduct0)) = 0;
conduct0 = abs(conduct0);


% Plot results
figure; 
set(gcf,'position', [300, 300, 1100, 400])
subplot(1,2,1), image(im,'CDataMapping','scaled'); 
colorbar; title('Original image'); xlabel('x'); ylabel('y');
subplot(1,2,2), 
image(conduct0, 'CDataMapping', 'scaled'); 
colorbar; title('Conductance matrix derived from original image');
xlabel('voxel i'); ylabel('voxel s');


%%
% To compare results with pointbypoint case
% generated using example_2D_conductance_pointbypoint.m
comparison = conduct0 - conductance(rowx_index,rowy_index);
sum(comparison(:))
figure; imagesc(comparison); colorbar; xlabel('voxel i'); ylabel('voxel s');
title('Comparison of point-by-point and using superposition')

