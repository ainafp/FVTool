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
im = I.img;
voxel_size = I.hdr.dime.pixdim(2:4);

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