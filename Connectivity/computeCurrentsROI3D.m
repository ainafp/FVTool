function RHSbcM = computeCurrentsROI3D(im_sizes, RHSbc, atlas, mask)
% Compute potentials 3D
%
% Author: Aina Frau-Pascual


% Define mask
Nx = im_sizes(1);
Ny = im_sizes(2);
Nz = im_sizes(3);
G = reshape((1:(Nx+2)*(Ny+2)*(Nz+2)), Nx+2, Ny+2, Nz+2);
rowx_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1), Nx*Ny*Nz, 1);

% Compute ROI list
roi_list0 = unique(atlas);                    % consider 0 is background
roi_list = roi_list0(roi_list0>0);            

% Define sink in ROI=2 of WM
%p0b = rowx_index(atlas==2);                  % and 2 is WM
p0b = rowx_index(mask>0);
p0 = p0b(round(length(p0b)/10*4));
clear p0b roi_list0 G Nx Ny Nz im_sizes

% Create currents for all ROI
RHSbcM = sparse(zeros([size(RHSbc,1),size(roi_list,1)]));
for i=1:size(roi_list,1)
    p2 = rowx_index(atlas==roi_list(i));
    RHSbc0 = RHSbc;
    RHSbc0(p0) = 1 / length(p0);             % define current i
    RHSbc0(p2) = RHSbc0(p2) - 1/length(p2);  % define current j
    RHSbcM(:,i) = RHSbc0;
end

end

