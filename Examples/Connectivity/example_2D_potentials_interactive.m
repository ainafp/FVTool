% Diffusion equation $\nabla \cdot (-D \nabla \phi) = \gamma$  
% for Phantom example

% Input image
im = phantom(100);

compute_connmatrix = false;

% Construct mesh structure
Lx = size(im, 1);                % domain length
Nx = size(im, 1);                % number of cells
Ly = size(im, 2);                % domain length
Ny = size(im, 2);                % number of cells
meshstruct = createMesh2D(Nx, Ny, Lx, Ly); % construct mesh
x = meshstruct.cellcenters.x;    % extract the cell center positions

% Define the boundary condition:
BC = createBC(meshstruct);       % all Neumann boundary condition structure

% Define the transfer coefficients:
D = createCellVariable(meshstruct, im);
Dave = harmonicMean(D);          % convert a cell variable to face variable
Mdiff = diffusionTerm(Dave);     % diffusion term
[Mbc, RHSbc] = boundaryCondition(BC);   % boundary condition discretization
M = Mdiff + Mbc;                 % matrix of coefficient for central scheme

% Define mask
Nx = Dave.domain.dims(1);
Ny = Dave.domain.dims(2);
G = reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
mnx = Nx*Ny;	
mny = Nx*Ny;
rowx_index = reshape(G(2:Nx+1,2:Ny+1),mnx,1);      % main diagonal x
rowy_index = reshape(G(2:Nx+1,2:Ny+1),mny,1);      % main diagonal y

% Compute conductance
p0 = rowx_index(round(mnx/2));

potentials = zeros(Nx, Ny, mny);
i = 1;
for p2=rowy_index'
    if p0~=p2
        RHSbc0 = RHSbc;
        RHSbc0(p0) = 1;                      % define current i
        RHSbc0(p2) = -1;                     % define current j
        c = solvePDE(meshstruct, M, RHSbc0); % solve for the central scheme
        potentials(:,:,i) = c.value(2:end-1, 2:end-1);
    end
    i = i + 1;
end

% Connectivity matrix
if compute_connmatrix
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
    figure(100); imagesc(abs(conduct0), [-1,1]); colorbar;
end

%%

% Conductance map
figure(101); image(im,'CDataMapping','scaled'); colorbar;

% Interactive conductance map
for i=1:100
    figure(101);
    [x, y] = ginput(1);
    a1 = round(y);
    b1 = round(x);
    p1 = sub2ind([Nx, Ny], a1, b1);
    
    c = zeros(Nx, Ny);
    for p2=1:mny
        [a2, b2] = ind2sub([Nx, Ny], p2);
        c(a2, b2) = abs(1 / ((potentials(a1,b1,p2) - potentials(a1,b1,p1)) ... 
                         - (potentials(a2,b2,p2) - potentials(a2,b2,p1))));
    end
    figure(102); imagesc((c), [0,0.5]); colormap jet; colorbar;

end
