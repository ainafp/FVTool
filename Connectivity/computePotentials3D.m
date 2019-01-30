function [potentials,mask_index0] = computePotentials3D(D, meshstruct, M, RHSbc)
% Compute potentials 3D

% Define mask
Nx = D.domain.dims(1);
Ny = D.domain.dims(2);
Nz = D.domain.dims(3);
G=reshape((1:(Nx+2)*(Ny+2)*(Nz+2)), Nx+2, Ny+2, Nz+2);
mnx = Nx*Ny*Nz;
rowx_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnx,1); % main diagonal x

% Compute conductance
aux = sum(D.xvalue(2:end-1,2:end-1,2:end-1,:),4) + ...
      sum(D.yvalue(2:end-1,2:end-1,2:end-1,:),4) + ...
      sum(D.zvalue(2:end-1,2:end-1,2:end-1,:),4);
mask = reshape(aux,mnx,1);
mask_index = rowx_index(mask>0);
mask_index0 = find(mask>0);
p0 = mask_index(round(length(mask_index)/2));

potentials = zeros(Nx, Ny, Nz, size(mask_index,1));
for i=1:size(mask_index,1) % mnx
    %p2=rowx_index(i);
    p2=mask_index(i);
    if p0~=p2
        RHSbc0 = RHSbc;
        RHSbc0(p0) = 1;                          % define current i
        RHSbc0(p2) = -1;                         % define current j
        pot_0 = solvePDE(meshstruct, M, RHSbc0); % solve for the central scheme
        potentials(:,:,:,i) = pot_0.value(2:end-1, 2:end-1, 2:end-1);
    end
end