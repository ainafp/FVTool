% Diffusion equation $\nabla \cdot (-D \nabla \phi) = \gamma$  
% Compute conductance from vectors in a 3D example

plot_input = false;

% Input image
%im = zeros(2,2,2,3,3,3) + 0.5;
im = zeros(3,3,3);
im(:,1,1) = 0.8;
im(:,:,1) = 0.8;
im(:,1,:) = 0.8;
im(3,3,3) = 0.8;
im(:,:,:,2) = 1 - im - 0.1;
im(:,:,:,3) = 0.1;
im(im>0.8) = 0;
size(im)

if plot_input
    figure; imagesc(im(:,:,1,1), [0,1]); colorbar;
    figure; imagesc(im(:,:,2,1), [0,1]); colorbar;
    figure; imagesc(im(:,:,3,1), [0,1]); colorbar;
    figure; imagesc(im(:,:,1,2), [0,1]); colorbar;
    figure; imagesc(im(:,:,2,2), [0,1]); colorbar;
    figure; imagesc(im(:,:,3,2), [0,1]); colorbar;
end

%%
% Construct mesh structure
Lx = size(im, 1);                   % domain length
Nx = size(im, 1);                   % number of cells
Ly = size(im, 2);                   % domain length
Ny = size(im, 2);                   % number of cells
Lz = size(im, 3);                   % domain length
Nz = size(im, 3);                   % number of cells
meshstruct = createMesh3D(Nx, Ny, Nz, Lx, Ly, Lz);         % construct mesh
x = meshstruct.cellcenters.x;           % extract the cell center positions

% Define the boundary condition:
BC = createBC(meshstruct);       % all Neumann boundary condition structure
[Mbc, RHSbc] = boundaryCondition(BC);   % boundary condition discretization

%%
% Define the transfer coefficients:
D = createFaceVariable(meshstruct, 0.0);
D.xvalue = im(:,:,:,1);
D.xvalue(size(D.xvalue,1)+1,:,:) = D.xvalue(end,:,:);
D.yvalue = im(:,:,:,2);
D.yvalue(:,size(D.yvalue,1)+1,:) = D.yvalue(:,end,:);
D.zvalue = im(:,:,:,3);
D.zvalue(:,:,size(D.zvalue,1)+1) = D.zvalue(:,:,end);
Mdiff = diffusionTerm(D);
M = Mdiff + Mbc;                 % matrix of coefficient for central scheme

%%
% Define mask
G = reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);
mnx = Nx*Ny*Nz;	mny = Nx*Ny*Nz; mnz = Nx*Ny*Nz;
rowx_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnx,1); % main diagonal x
rowy_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mny,1); % main diagonal y
rowz_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnz,1); % main diagonal z

%%
% Compute conductance
conductance = zeros(size(RHSbc, 1), size(RHSbc, 1));
for p1=rowx_index'
    for p2=rowy_index'
        if p2>=p1
            RHSbc0 = RHSbc;
            RHSbc0(p1) = 1;                          % define current i
            RHSbc0(p2) = -1;                         % define current j
            c = solvePDE(meshstruct, M, RHSbc0);     % solve for the central scheme
            conductance(p1, p2) = abs(1 / (c.value(p1) - c.value(p2)));
            conductance(p2, p1) = conductance(p1, p2);
        end
    end
end

%%
% Plot results
figure; image(conductance(rowx_index,rowy_index), 'CDataMapping', 'scaled'); colorbar
