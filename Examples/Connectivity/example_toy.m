% You can solve a diffusion equation
% $ \nabla. (-D \nabla \phi) = \gamma $ 

L = 145;  % domain length
im = I;
im = ones(145);
Nx = size(im, 1); % number of cells
Ny = size(im, 2); % number of cells
meshstruct = createMesh2D(Nx,Ny, L,L);
%m = createMesh3D(Nx,Nx,Nx, L,L,L);
%m = createMesh1D(Nx, L);
x = meshstruct.cellcenters.x; % extract the cell center positions
% The next step is to define the boundary condition:
BC = createBC(meshstruct); % all Neumann boundary condition structure
%BC.left.a = 0; BC.left.b=1; % switch the left boundary to Dirichlet
%BC.left.c=0; % value = 0 at the left boundary
%BC.right.a = 0; BC.right.b=1; % switch the right boundary to Dirichlet
%BC.right.c=1; % value = 1 at the right boundary
%BC2 = createBC(meshstruct); % all Neumann boundary condition structure
%BC2.left.a(:) = 0; BC2.left.b(:)=1; BC2.left.c(:)=1; % Dirichlet for the left boundary
%BC2.right.a(:) = 0; BC2.right.b(:)=1; BC2.right.c(:)=0; % right boundary

% Now we define the transfer coefficients:
%D_val = 1.0; % diffusion coefficient value
%D = createCellVariable(meshstruct, D_val); % assign dif. coef. to all the cells    
D = createCellVariable(meshstruct, im); %, BC2);ex0

Dave = harmonicMean(D); % convert a cell variable to face variable
Mdiff = diffusionTerm(Dave); % diffusion term
[Mbc, RHSbc] = boundaryCondition(BC); % boundary condition discretization
M = Mdiff+Mbc; % matrix of coefficient for central scheme


RHSbc(sub2ind([Ny - 2,Nx-2], 72, 9)) = 1;
result = zeros(size(RHSbc));
for i=1:1 %size(RHSbc)
    RHSbc(sub2ind([Ny,Nx], 70, 138)) = -1;
    c = solvePDE(meshstruct, M, RHSbc); % solve for the central scheme
    result(i) = c.value(1);
    visualizeCells(c)
    aux = c.value(c.value~=0)';
    aux = aux(~isnan(aux))'
end
%result
%visualizeCells(c); % visualize the results
figure; image(reshape(log(result), [Nx + 2, Ny + 2]),'CDataMapping','scaled'); colorbar;

%%
D_val = C1; % diffusion coefficient value
D = createCellVariable(meshstruct, D_val); % assign dif. coef. to all the cells
Dave = harmonicMean(D); % convert a cell variable to face variable
Mdiff = diffusionTerm(Dave); % diffusion term
BC = createBC(meshstruct); % all Neumann boundary condition structure
BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % Dirichlet for the left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
[Mbc, RHSbc] = boundaryCondition(BC); % boundary condition discretization
M = Mdiff+Mbc; % matrix of coefficient for central scheme
c = solvePDE(meshstruct, M, RHSbc); % solve for the central scheme
visualizeCells(c); % visualize the results

figure; image(c.value,'CDataMapping','scaled'); colorbar;
