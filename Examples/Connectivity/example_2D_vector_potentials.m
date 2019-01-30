% Phantom 2D small example using vectors instead of scalars
% computing potential and conductance between 2 points
% Author: Aina Frau-Pascual


% Construct mesh structure and boundary conditions
Lx = 10; % domain length
Nx = 10; % number of cells
Ly = 10; % domain length
Ny = 10; % number of cells
meshstruct = createMesh2D(Nx, Ny, Lx, Ly); % construct mesh
x = meshstruct.cellcenters.x; % extract the cell center positions
BC = createBC(meshstruct); % all Neumann boundary condition structure
[Mbc, RHSbc] = boundaryCondition(BC); % boundary condition discretization

% In this example, I input tensors here manually
D = createFaceVariable(meshstruct, 0.0);
D.yvalue = [1,2,3,4,5,6,7,8,9,8,7; 7,6,5,4,3,2,1,2,3,4,7;
            5,6,7,8,9,8,7,6,5,4,7; 3,2,1,2,3,4,5,6,7,8,7;
            9,8,7,6,5,4,3,2,1,2,7; 1,2,3,4,5,6,7,8,9,8,7;
            7,6,5,4,3,2,1,2,3,4,7; 5,6,7,8,9,8,7,6,5,4,7;
            3,2,1,2,3,4,5,6,7,8,7; 9,8,7,6,5,4,3,2,1,2,7];
D.xvalue = [7,6,5,4,3,2,1,2,3,4;   5,6,7,8,9,8,7,6,5,4;
            3,2,1,2,3,4,5,6,7,8;   9,8,7,6,5,4,3,2,1,2;
            1,2,3,4,5,6,7,8,9,8;   7,6,5,4,3,2,1,2,3,4;
            5,6,7,8,9,8,7,6,5,4;   3,2,1,2,3,4,5,6,7,8;
            9,8,7,6,5,4,3,2,1,2;   1,2,3,4,5,6,7,8,9,8;
            9,8,7,6,5,4,3,2,1,2];
Mdiff = diffusionTerm(D);
M = Mdiff + Mbc; % matrix of coefficient for central scheme

% Define currents (right hand side equation)
p1 = sub2ind([Ny + 2, Nx + 2], 7, 9);
p2 = sub2ind([Ny + 2, Nx + 2], 4, 3);
RHSbc(p1) = 1;
RHSbc(p2) = -1;

% Solve PDE
c = solvePDE(meshstruct, M, RHSbc); % solve for the central scheme

% Compute conductance
conductance = 1 / (c.value(p1) - c.value(p2));

% Plot results
figure; image(c.value, 'CDataMapping', 'scaled'); colorbar;

% Plot results
figure; 
set(gcf,'position', [300, 300, 1100, 400])
subplot(1,2,1), image(D.xvalue + D.yvalue,'CDataMapping','scaled'); 
colorbar; title('Original image (sum tensor dim)'); xlabel('x'); ylabel('y');
subplot(1,2,2), 
image(conductance(rowx_index,rowy_index), 'CDataMapping', 'scaled'); 
colorbar; title('Conductance matrix derived from original image');
xlabel('voxel i'); ylabel('voxel s');

