% Test point by point conductance computation, in a 2D vector example
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
D.xvalue(size(D.xvalue,1)+1,:) = D.xvalue(end,:);
D.yvalue = im(:,:,2);
D.yvalue(:,size(D.yvalue,1)+1) = D.yvalue(:,end);
Mdiff = diffusionTerm(D);
M = Mdiff + Mbc; % matrix of coefficient for central scheme

% Define mask
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
        end
    end
end

% Plot results
conductance(isnan(conductance)) = 0;
figure; image(conductance(rowx_index,rowy_index), 'CDataMapping', 'scaled'); colorbar
