function [xvalue, yvalue, zvalue] = extendBoundaryTensor3D(im)

Nx = size(im, 1); % number of cells
Ny = size(im, 2); % number of cells
Nz = size(im, 3); % number of cells
Tx = size(im, 4); % size tensor

% Extend image
xvalue = zeros(Nx + 4, Ny + 4, Nz + 4, Tx);
xvalue(3:end-2, 3:end-2, 3:end-2, :) = im(:,:,:,:,1);
xvalue(1,3:end-2,3:end-2,:) = im(1,:,:,:,1);
xvalue(2,3:end-2,3:end-2,:) = im(1,:,:,:,1);
xvalue(end,3:end-2,3:end-2,:) = im(end,:,:,:,1);
xvalue(end-1,3:end-2,3:end-2,:) = im(end,:,:,:,1);
xvalue(:,1,:,:) = xvalue(:,3,:,:);
xvalue(:,2,:,:) = xvalue(:,3,:,:);
xvalue(:,end,:,:) = xvalue(:,end-2,:,:);
xvalue(:,end-1,:,:) = xvalue(:,end-2,:,:);
xvalue(:,:,1,:) = xvalue(:,:,3,:);
xvalue(:,:,2,:) = xvalue(:,:,3,:);
xvalue(:,:,end,:) = xvalue(:,:,end-2,:);
xvalue(:,:,end-1,:) = xvalue(:,:,end-2,:);

yvalue = zeros(Nx + 4, Ny + 4, Nz + 4, Tx);
yvalue(3:end-2, 3:end-2, 3:end-2, :) = im(:,:,:,:,2);
yvalue(1,3:end-2,3:end-2,:) = im(1,:,:,:,2);
yvalue(2,3:end-2,3:end-2,:) = im(1,:,:,:,2);
yvalue(end,3:end-2,3:end-2,:) = im(end,:,:,:,2);
yvalue(end-1,3:end-2,3:end-2,:) = im(end,:,:,:,2);
yvalue(:,1,:,:) = yvalue(:,3,:,:);
yvalue(:,2,:,:) = yvalue(:,3,:,:);
yvalue(:,end,:,:) = yvalue(:,end-2,:,:);
yvalue(:,end-1,:,:) = yvalue(:,end-2,:,:);
yvalue(:,:,1,:) = yvalue(:,:,3,:);
yvalue(:,:,2,:) = yvalue(:,:,3,:);
yvalue(:,:,end,:) = yvalue(:,:,end-2,:);
yvalue(:,:,end-1,:) = yvalue(:,:,end-2,:);

zvalue = zeros(Nx + 4, Ny + 4, Nz + 4, Tx);
zvalue(3:end-2, 3:end-2, 3:end-2, :) = im(:,:,:,:,3);
zvalue(1,3:end-2,3:end-2,:) = im(1,:,:,:,3);
zvalue(2,3:end-2,3:end-2,:) = im(1,:,:,:,3);
zvalue(end,3:end-2,3:end-2,:) = im(end,:,:,:,3);
zvalue(end-1,3:end-2,3:end-2,:) = im(end,:,:,:,3);
zvalue(:,1,:,:) = zvalue(:,3,:,:);
zvalue(:,2,:,:) = zvalue(:,3,:,:);
zvalue(:,end,:,:) = zvalue(:,end-2,:,:);
zvalue(:,end-1,:,:) = zvalue(:,end-2,:,:);
zvalue(:,:,1,:) = zvalue(:,:,3,:);
zvalue(:,:,2,:) = zvalue(:,:,3,:);
zvalue(:,:,end,:) = zvalue(:,:,end-2,:);
zvalue(:,:,end-1,:) = zvalue(:,:,end-2,:);
