function [fa, md, tensor, voxel_size]= read_tensors_fsl(folder, basename)
% Read files generated with FSL and return FA, MD, tensors and 
% voxel size
%
% Author: Aina Frau-Pascual


%%% Prepare inputs

% Input tensors from FSL dtifit outputs
fa = load_nii(fullfile(folder, [basename '_FA.nii.gz'])); 
fa = fa.img;
md = load_nii(fullfile(folder, [basename '_MD.nii.gz']));
md = md.img;
I = load_nii(fullfile(folder, [basename '_tensor.nii.gz']));
tensor = I.img;
voxel_size = I.hdr.dime.pixdim(2:4);

end