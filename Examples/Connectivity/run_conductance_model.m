function run_conductance_model(data_folder)
% Compute conductance maps from tensors computed with DSIstudio
% This is a generic example.
% Adapt paths and file names to your own!
%
% Example: run_conductance_model(your_path)
%
% Author: Aina Frau-Pascual


%%% Load folders and code

% Add code folders to path
folder_nii='NIFTI';
addpath(genpath(folder_nii));
folder_nii='FVT4DWI';
addpath(genpath(folder_nii));


%%% Prepare inputs

% Input tensors
file_name = fullfile(data_folder, 'tensors_dsistudio.fib.gz');
[fa, md, im, voxel_size]= read_fib(file_name);
im_size = size(im);

% Parcellation
atlas_name = fullfile(data_folder, 'your_atlas.nii.gz');
atlas_obj = load_nii(atlas_name);
atlas = double(atlas_obj.img);
atlas = rot90(2,atlas); % WARNING! check that orientation tensors/atlas match
% mask of GM from your parcellation
mask2 = atlas;
mask2(mask2>0) = 1;
clear atlas_obj

% Mask
% mask of whole-brain from nodif_brain_mask
mask_name = fullfile(data_folder, 'your_brain_mask.nii.gz');
mask_obj = load_nii(mask_name);
mask0 = double(mask_obj.img);
mask0 = rot90(mask0,2);
% mask of WM from wm_mask
mask_name = fullfile(data_folder, 'your_wm_mask.nii.gz');
mask_obj = load_nii(mask_name);
mask1 = double(mask_obj.img);
mask1 = rot90(mask1,2);
% mask of WM and GM together
mask = mask0 .* (mask1 + mask2);
mask(mask>0) = 1;


%%% Compute conductance

% Compute diffusion matrix
[meshstruct, D, M, RHSbc] = computeDiffusionMatrix3D(im.*mask, voxel_size);
sprintf('diffusion matrix computed')
clear meshstruct im D

tic 

% Compute currents per ROI
RHSbcM = computeCurrentsROI3D(im_size(1:3), M, sparse(RHSbc), atlas.*mask);

% Compute inversion
% WARNING! Make sure you have enough RAM to do the inversion
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

% Compute conductance matrix
conn = computeConductanceMatrix(potentials, atlas.*mask);
conn_fn = fullfile(data_folder, 'conductance_matrix.mat');
save(conn_fn, 'conn');

toc

% Save results if you want to look at the potentials
% WARNING! The file will be big
filename = fullfile(data_folder, 'conductance.mat');
save(filename, '-v7.3')
sprintf('file saved')

end
