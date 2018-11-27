function run_global_connectivity_hcp(varargin)
% Compute global connectivity from tensors computed with DSIstudio


%%% Load folders and code

% Add code folders to path
folder='Code/FVTool';
addpath(genpath(folder));
folder_nii='Code/NIFTI';
addpath(genpath(folder_nii));

% Subject name
args = varargin{:};
subject = num2str(args)

% Data folder
hcp_folder = 'Data/HCP_WashU-UMN/';
subject_folder = [subject '_3T_Diffusion_preproc/' subject '/T1w/'];
data_folder = [hcp_folder subject_folder '/Diffusion/'];


%%% Prepare inputs

% Input image
file_name = [data_folder 'data.src.gz.012fy.dti.fib.gz'];
[fa, md, im, voxel_size]= read_fib(file_name);
im_size = size(im);

% Parcellation
atlas_name = [data_folder './../wmparc_1.25.nii.gz']; % from FreeSurfer
atlas_obj = load_nii(atlas_name);
atlas = double(atlas_obj.img);
atlas = rot90(atlas,2);
clear atlas_obj
% Consider all WM in 2 regions lh and rh
atlas((atlas>2999) & (atlas<3036) | (atlas==5001)) = 2;
atlas((atlas>3999) & (atlas<4036) | (atlas==5002)) = 41;

% Mask
mask_name = [data_folder 'nodif_brain_mask.nii.gz'];
mask_obj = load_nii(mask_name);
mask0 = double(mask_obj.img);
mask0 = rot90(mask0,2);
% Remove ventricles from mask
mask2 = atlas;
mask2((atlas==0) | (atlas==4) | (atlas==5) | (atlas==14) | (atlas==15) ...
   | (atlas==24) | (atlas==43) | (atlas==44) | (atlas==72) | (atlas==80)) = 0;
mask2(mask2>0) = 1;
mask = mask0.*mask2;
clear mask_obj mask0 mask2


%%% Compute global connectivity

% Compute diffusion matrix
[M, RHSbc] = computeDiffusionMatrix3D(im.*mask, voxel_size);
sprintf('diffusion matrix computed')
clear meshstruct im

tic 

% Compute currents per ROI
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

% Take potentials final result
conn = computeConductanceMatrix(potentials, atlas.*mask);
conn_fn = sprintf([data_folder 'conductance_matrix.mat']);
save(conn_fn, 'conn');

toc

% Save results
filename = sprintf([data_folder 'conductance.mat']);
save(filename, '-v7.3')
sprintf('file saved')

end

