% Diffusion equation $ \nabla. (-D \nabla \phi) = \gamma $  
% for Phantom 2D small example or ones small example
% use superposition property and save potentials to construct
% connectivity matrix

eps = 0.000001;

% Input image
%im = phantom(10) + eps;
%im = im(2:end-1,3:end-2); % + eps;
im = zeros(3);
im(1, :) = 1;
im(:, 1) = 1;
im(3, 3) = 1;
figure(101); imagesc(im, [0,1]); colorbar
saveas(gcf, strcat('example_phantom_3x3.eps'), 'epsc')

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
%cond = zeros(size(RHSbc, 1), size(RHSbc, 1));

p0 = rowx_index(round(mnx/2))
%for p0=rowx_index' %(end)
    potentials = zeros(Nx, Ny, mny);
    i = 1;
    for p2=rowy_index'
        if p0~=p2
            RHSbc0 = RHSbc;
            RHSbc0(p0) = 1;                       % define current i
            RHSbc0(p2) = -1;                      % define current j
            c = solvePDE(meshstruct, M, RHSbc0);  % solve for the central scheme
            potentials(:,:,i) = c.value(2:end-1, 2:end-1);
            %figure(100); image(c.value(2:end-1, 2:end-1), 'CDataMapping', 'scaled'); colorbar;
            %pause(0.5);
        end
        i = i + 1;
    end

    %figure(103); 
    %for i=1:mny
    %    subplot(Nx, Ny, i); imagesc(potentials(:,:,i)); colorbar;
    %end
    
    conduct0 = zeros(mnx, mny);
    for p1=1:mnx
        for p2=1:mny
            if p2>p1
                pot = (potentials(:,:,p2) - potentials(:,:,p1));
                %figure(102); subplot(Nx, Ny, i); imagesc(pot);                
                [a1, b1] = ind2sub(size(pot), p1);
                [a2, b2] = ind2sub(size(pot), p2);
                
                %im = ones(Nx, Ny) ./ (pot - pot(a2, b2) * ones(Nx, Ny))
                %if sum(im>0)<(Nx*Ny/2)
                %    im = -im;
                %end
                %figure(102); imagesc(im, [0,1]); colorbar;
                %pause(0.5);
                
                conduct0(p1, p2) = 1 / (pot(a1, b1) - pot(a2, b2));
                conduct0(p2, p1) = 1 / (pot(a2, b2) - pot(a1, b1));
            end
        end
    end

    % Plot results
    figure(101); imagesc(abs(conduct0), [-1,1]); colorbar
    saveas(gcf, strcat('example_phantom_3x3_connmatrix.eps'), 'epsc')
    aa = abs(conduct0) - abs(conductance(rowx_index,rowy_index));
    sum(sum(aa))
    %figure(105); imagesc(aa); colorbar
    %pause(0.5);
%end

%savefig(strcat('example_phantom_8x9_sup_',int2str(p0),'.fig'))
%saveas(gcf, strcat('example_phantom_8x9_sup_',int2str(p0),'.eps'), 'epsc')

%figure; image(im,'CDataMapping','scaled'); colorbar;
%savefig('example_phantom_8x9.fig')
%saveas(gcf, 'example_phantom_8x9.eps', 'epsc')

