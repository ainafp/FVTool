function [fa, md, tensor, voxel_size]= read_fib(file_name)
% Read .fib generated with DSI Studio and return FA, MD, tensors and 
% voxel size
%
% Author: Aina Frau-Pascual, modified from public code in 
% http://dsi-studio.labsolver.org/Manual/export-data-to-matlab


if ~exist('file_name')
    file_name =  uigetfile('*.fib.gz');
end
if file_name == 0
    image = [];
    return
end
file_name_unzip = gunzip(file_name);
[pathstr, name0, ext] = fileparts(file_name);
[pathstr, name, ext] = fileparts(file_name_unzip{1});
file_mat = strcat(pathstr, '/', name, '.mat');
movefile(strcat(pathstr,'/',name0),file_mat);
fib = load(file_mat);

max_fib = 0;
for i = 1:10
    if isfield(fib,strcat('fa',int2str(i-1)))
        max_fib = i;
    else
        break;
    end
end

fa = zeros([fib.dimension max_fib]);
md = zeros([fib.dimension max_fib]);
index = zeros([fib.dimension max_fib]);

for i = 1:max_fib
    eval(strcat('fa(:,:,:,i) = reshape(fib.fa',int2str(i-1),',fib.dimension);'));
    eval(strcat('md(:,:,:,i) = reshape(fib.md,fib.dimension);'));
end

tensor(:,:,:,1,1) = reshape(fib.txx,fib.dimension);
tensor(:,:,:,1,2) = reshape(fib.txy,fib.dimension);
tensor(:,:,:,1,3) = reshape(fib.txz,fib.dimension);
tensor(:,:,:,2,1) = reshape(fib.txy,fib.dimension);
tensor(:,:,:,2,2) = reshape(fib.tyy,fib.dimension);
tensor(:,:,:,2,3) = reshape(fib.tyz,fib.dimension);
tensor(:,:,:,3,1) = reshape(fib.txz,fib.dimension);
tensor(:,:,:,3,2) = reshape(fib.tyz,fib.dimension);
tensor(:,:,:,3,3) = reshape(fib.tzz,fib.dimension);

odf_vertices = fib.odf_vertices;
odf_faces = fib.odf_faces;
voxel_size = fib.voxel_size;
%delete(strcat(name,'.mat'));
end