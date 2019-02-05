function [xvalue, yvalue, zvalue] = extendBoundaryTensor3D(im)
% This function extends the 3D tensors in the boundary. When using tensors,
% you need to use more neighbors to compute discretized matrix
%
% Author: Aina Frau-Pascual


Nx = size(im, 1); % number of cells
Ny = size(im, 2); % number of cells
Nz = size(im, 3); % number of cells
Tx = size(im, 4); % size tensor

% Extend image
xvalue = zeros(Nx + 2, Ny + 2, Nz + 2, Tx);
xvalue(2:end-1, 2:end-1, 2:end-1, :) = im(:,:,:,:,1);
xvalue(1,2:end-1,2:end-1,:) = im(1,:,:,:,1);
xvalue(end,2:end-1,2:end-1,:) = im(end,:,:,:,1);
xvalue(:,1,:,:) = xvalue(:,2,:,:);
xvalue(:,end,:,:) = xvalue(:,end-1,:,:);
xvalue(:,:,1,:) = xvalue(:,:,2,:);
xvalue(:,:,end,:) = xvalue(:,:,end-1,:);

yvalue = zeros(Nx + 2, Ny + 2, Nz + 2, Tx);
%yvalue(2:end-1, 2:end-1, 2:end-1, :) = im(:,:,:,:,2);
%yvalue(1,2:end-1,2:end-1,:) = im(1,:,:,:,2);
%yvalue(end,2:end-1,2:end-1,:) = im(end,:,:,:,2);
yvalue(:,1,:,:) = yvalue(:,2,:,:);
yvalue(:,end,:,:) = yvalue(:,end-1,:,:);
yvalue(:,:,1,:) = yvalue(:,:,2,:);
yvalue(:,:,end,:) = yvalue(:,:,end-1,:);

zvalue = zeros(Nx + 2, Ny + 2, Nz + 2, Tx);
%zvalue(2:end-1, 2:end-1, 2:end-1, :) = im(:,:,:,:,3);
%zvalue(1,2:end-1,2:end-1,:) = im(1,:,:,:,3);
%zvalue(end,2:end-1,2:end-1,:) = im(end,:,:,:,3);
zvalue(:,1,:,:) = zvalue(:,2,:,:);
zvalue(:,end,:,:) = zvalue(:,end-1,:,:);
zvalue(:,:,1,:) = zvalue(:,:,2,:);
zvalue(:,:,end,:) = zvalue(:,:,end-1,:);
