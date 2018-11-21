function computeConductanceFromTensors(data_folder)
% Compute global connectivity from tensors computed with DSIstudio


%%% Load folders and code

% Add code folders to path
folder_nii='NIFTI';
addpath(genpath(folder_nii));


%%% Prepare inputs

% Input tensors. In this example from DSI Studio
file_name = fullfile(data_folder, 'data.src.gz.012fy.dti.fib.gz');
[fa, md, im, voxel_size]= read_fib(file_name);

% Parcellation
atlas_name = fullfile(data_folder, 'your_atlas.nii.gz');
atlas_obj = load_nii(atlas_name);
atlas = double(atlas_obj.img);
atlas = rot90(atlas);
clear atlas_obj

% Masks
% mask of whole-brain from nodif_brain_mask
mask_name = fullfile(data_folder, 'nodif_brain_mask.nii.gz');
mask_obj = load_nii(mask_name);
mask0 = double(mask_obj.img);
mask0 = rot90(mask0);
% mask of WM from wm_mask
mask_name = fullfile(data_folder, 'wm_mask.nii.gz');
mask_obj = load_nii(mask_name);
mask1 = double(mask_obj.img);
mask1 = rot90(mask1);
% mask of GM from your parcellation
mask2 = atlas0;
mask2(mask2>0) = 1;
% mask of WM and GM together
mask = mask0 .* (mask1 + mask2);
mask(mask>0) = 1;
clear mask_obj mask0 mask1 mask2
% CAUTION: make sure tensors and parcellation correspond


%%% Compute potentials

% Compute diffusion matrix
[M, RHSbc] = computeDiffusionMatrix3D(im.*mask, voxel_size);
sprintf('diffusion matrix computed')
clear meshstruct im

% Compute currents per ROI
im_size = size(im);
RHSbcM = computeCurrentsROI3D(im_size(1:3), sparse(RHSbc), atlas.*mask, mask);

% Compute inversion
x = M\RHSbcM;
sprintf('inversion computed')
clear M RHSbc RHSbcM

% Reshape result
phival = reshape(full(x), [im_size(1:3)+2 size(x,2)]);
sprintf('reshape computed')
clear x

% Take potentials final result
potentials = phival(2:end-1, 2:end-1, 2:end-1, :);
clear phival
sprintf('potentials computed')

% Compute conductance
conn = computeConductanceMatrix(potentials, atlas.*mask);

% Save results
conn_fn = fullfile(data_folder, 'conductance_matrix.mat');
save(conn_fn, 'conn');
filename = fullfile(data_folder, 'conductance.mat');
save(filename, '-v7.3') % change according to Matlab version
sprintf('file saved')

end
