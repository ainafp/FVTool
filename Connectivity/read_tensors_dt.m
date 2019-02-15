function [fa, tensor, voxel_size]= read_tensors_dt(folder, basename)
% Read files generated with TrackVis Diffusion Toolkit and 
% return FA, tensors and voxel size
% MD could be computed as (e1+e2+e3)/3
%
% Author: Aina Frau-Pascual


%%% Prepare inputs

% Input tensors from FSL dtifit outputs
fa = load_nii(fullfile(folder, [basename '_fa.nii'])); 
fa = fa.img;
I = load_nii(fullfile(folder, [basename '_tensor.nii']));
tensor = I.img;
voxel_size = I.hdr.dime.pixdim(2:4);

end