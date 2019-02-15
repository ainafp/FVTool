% Example using the phantom Fibercup (original version)
% Data can be downloaded in http://www.tractometer.org/original_fibercup/data/
%
% Author: Aina Frau-Pascual


% Load Input image, generated with DSI Studio
%dsi_studio_run --action=src --source=data.nii.gz --bval=bvals --bvec=bvecs --output=data.src.gz'
%dsi_studio_run --action=rec --source=data.src.gz --mask=nodif_brain_mask.nii.gz --method=1 --output_tensor=1 --output_dif=1 --output=data.src.gz.012fy.dti.fib.gz';

% Alternatively, you could generate your own tensors and load them as
% im = (im_x, im_y, im_z, tensor_x, tensor_y)
folder = 'your_folder';
file_name = fullfile(folder, 'data.nii.gz.src.gz.012fz.dti.fib.gz'); 
[fa0, md, im0, voxel_size] = read_fib(file_name);

% Choose to use the whole image '0' (default), a medium size image '1' 
% or a small size image '2'
image_size = 0;

% Load the image you chose, with its corresponding FA map
% We are here removing one dimension a bit roughly, only for testing
if image_size==2
    ax = 32:36; ay = 22:26;
elseif image_size==1
    ax = 31:50; ay = 13:33;
else
    ax = 1:size(im0,1); ay = 1:size(im0,2);
end
im = squeeze(im0(ax, ay, :, :, :));
fa = fa0(ax, ay, :);

% Compute potentials
Nx = size(im, 1); 
Ny = size(im, 2); 
Nz = size(im, 3); 
[meshstruct, D, M, RHSbc] = computeDiffusionMatrix3D(im, size(im));
[potentials, mask_index] = computePotentials3D(D, meshstruct, M, RHSbc);

% Plot conductance interactively
dimension = 3;
slice = 2; 
interactiveConductance3Dslice(fa, [Nx,Ny,Nz], potentials, mask_index, dimension, slice);

