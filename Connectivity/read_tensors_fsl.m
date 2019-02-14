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
im = I.img;
voxel_size = I.hdr.dime.pixdim;

tensor(:,:,:,1,1) = im(:,:,:,1); % xx
tensor(:,:,:,1,2) = im(:,:,:,2); % xy
tensor(:,:,:,1,3) = im(:,:,:,3); % xz
tensor(:,:,:,2,1) = im(:,:,:,2); % xy
tensor(:,:,:,2,2) = im(:,:,:,4); % yy
tensor(:,:,:,2,3) = im(:,:,:,5); % yz
tensor(:,:,:,3,1) = im(:,:,:,3); % xz
tensor(:,:,:,3,2) = im(:,:,:,5); % yz
tensor(:,:,:,3,3) = im(:,:,:,6); % zz

end