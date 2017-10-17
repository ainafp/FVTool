% Diffusion equation $ \nabla. (-D \nabla \phi) = \gamma $ 
% for Phantom example

eps = 0.000001;

% Input image
file_name = '/autofs/space/atlas_001/users/afrau/Data/Fibercup/acq-averaged_b-1500.nii.gz.src.gz.012fz.dti.fib.gz'; 
[fa tensor odf_vertices odf_faces] = read_fib(file_name);
im = tensor;

% Construct mesh structure
Lx = size(im, 1); % domain length
Nx = size(im, 1); % number of cells
Ly = size(im, 2); % domain length
Ny = size(im, 2); % number of cells
Lz = size(im, 3); % domain length
Nz = size(im, 3); % number of cells
Tx = size(im, 4); % size tensor
meshstruct = createMesh3D(Nx, Ny, Nz, Lx, Ly, Lz); 

% Define the boundary condition:
BC = createBC(meshstruct); % all Neumann boundary condition structure
[Mbc, RHSbc0] = boundaryCondition(BC); % boundary condition discretization

% Define the transfer coefficients:
[xvalue, yvalue, zvalue] = extendBoundaryTensor3D(im);
D = FaceVariable(meshstruct, xvalue, yvalue, zvalue);
[Mdiff, Mx, My, Mz] = diffusionTermTensor3D(D);
M = Mdiff + Mbc; % matrix of coefficient for central scheme

%% Conductance
% Define mask
Nx = D.domain.dims(1)+2;
Ny = D.domain.dims(2)+2;
Nz = D.domain.dims(3)+2;
Nt = (Nx+2)*(Ny+2)*(Nz+2);
G=reshape((1:Nt), Nx+2, Ny+2, Nz+2);
mnx = (Nx-2)*(Ny-2)*(Nz-2);
rowx_index = reshape(G(3:Nx,3:Ny,3:Nz),mnx,1); % main diagonal x

% Compute conductance
p0 = rowx_index(round(mnx/2));

potentials = zeros(Nx-2, Ny-2, Nz-2, mnx);
i = 1;
for p2=rowx_index'
    p2
    i
    if p0~=p2
        RHSbc0 = RHSbc;
        RHSbc0(p0) = 1;                      % define current i
        RHSbc0(p2) = -1;                     % define current j
        c = solvePDE(meshstruct, M, RHSbc0); % solve for the central scheme
        potentials(:,:,:,i) = c.value(2:end-1, 2:end-1, 2:end-1);
    end
    i = i + 1;
end

%% Connectivity matrix
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
    if saveims
        savefig(strcat(name, '_connectivity_matrix.fig'))
        saveas(gcf, strcat(name, '_connectivity_matrix.eps'), 'epsc')
    end
end

%% Interactive conductance map

% Conductance map
figure(101); image(im,'CDataMapping','scaled'); colorbar;
if saveims
    savefig(strcat(name, '.fig'));
    saveas(gcf, strcat(name, '.eps'), 'epsc');
end

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
    if saveims
        savefig(strcat(name, '_conductance_p', int2str(p1), '.fig'));
        saveas(gcf, strcat(name, '_conductance_p', int2str(p1), '.eps'), 'epsc');
    end
end

%% Check with 3 points

% Define currents (right hand side equation)
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

%dsistudio_trk = "/autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/DSI-Studio/0.9/bin";
%dsistudio = "/autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/DSI-Studio/20170531";
