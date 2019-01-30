% Test superposition property in conductance computation, 
% in a 2D vector example
% Author: Aina Frau-Pascual


% Input image
im = zeros(3) + 0.0;
im(1, :) = 0.8;
im(:, 1) = 0.8;
im(:,:,2) = 1 - im(:,:,1);
im(im==1) = 0;

% Construct mesh structure
Lx = size(im, 1); % domain length
Nx = size(im, 1); % number of cells
Ly = size(im, 2); % domain length
Ny = size(im, 2); % number of cells
meshstruct = createMesh2D(Nx, Ny, Lx, Ly); % construct mesh
x = meshstruct.cellcenters.x; % extract the cell center positions

% Define the boundary condition:
BC = createBC(meshstruct); % all Neumann boundary condition structure
[Mbc, RHSbc] = boundaryCondition(BC); % boundary condition discretization

% Define the transfer coefficients:
D = createFaceVariable(meshstruct, 0.0);
D.xvalue = im(:,:,1);
D.xvalue(size(D.xvalue,1)+1,:,:) = D.xvalue(end,:,:);
D.yvalue = im(:,:,2);
D.yvalue(:,size(D.yvalue,1)+1,:) = D.yvalue(:,end,:);
Mdiff = diffusionTerm(D);
M = Mdiff + Mbc; % matrix of coefficient for central scheme

% Define mask
G = reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
mnx = Nx*Ny;	
mny = Nx*Ny;
rowx_index = reshape(G(2:Nx+1,2:Ny+1),mnx,1); % main diagonal x
rowy_index = reshape(G(2:Nx+1,2:Ny+1),mny,1); % main diagonal y

% Compute potentials and conductance
p0 = rowx_index(round(mnx/3)); % choose a sink
potentials = zeros(Nx, Ny, mny);
i = 1;
for p2=rowy_index'
    if p0~=p2
        RHSbc0 = RHSbc;
        RHSbc0(p0) = 1;                       % define current i
        RHSbc0(p2) = -1;                      % define current j
        c = solvePDE(meshstruct, M, RHSbc0);  % solve for the central scheme
        potentials(:,:,i) = c.value(2:end-1, 2:end-1);
    end
    i = i + 1;
end

conduct0 = zeros(mnx, mny);
for p1=1:mnx
    for p2=1:mny
        if p2>p1
            pot = potentials(:,:,p2) - potentials(:,:,p1);
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
subplot(1,2,1), image(sum(im,3),'CDataMapping','scaled'); 
colorbar; title('Original image (sum tensor dim)'); xlabel('x'); ylabel('y');
subplot(1,2,2), 
image(conductance(rowx_index,rowy_index), 'CDataMapping', 'scaled'); 
colorbar; title('Conductance matrix derived from original image');
xlabel('voxel i'); ylabel('voxel s');


%%
% Compare to point-by-point matrix
comparison = conduct0 - conductance(rowx_index,rowy_index);
sum(sum(comparison))
figure; image(comparison, 'CDataMapping', 'scaled'); 
colorbar; xlabel('voxel i'); ylabel('voxel s');
title('Comparison of point-by-point and using superposition')
