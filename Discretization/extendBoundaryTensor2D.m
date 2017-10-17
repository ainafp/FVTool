function [xvalue, yvalue] = extendBoundaryTensor2D(im)

Nx = size(im, 1); % number of cells
Ny = size(im, 2); % number of cells
Tx = size(im, 3); % size tensor

% Extend image
xvalue = zeros(Nx + 4, Ny + 4, Tx);
xvalue(3:end-2, 3:end-2, :) = im(:,:,:,1);
xvalue(1,3:end-2,:) = im(1,:,:,1);
xvalue(2,3:end-2,:) = im(1,:,:,1);
xvalue(end,3:end-2,:) = im(end,:,:,1);
xvalue(end-1,3:end-2,:) = im(end,:,:,1);
xvalue(:,1,:) = xvalue(:,3,:);
xvalue(:,2,:) = xvalue(:,3,:);
xvalue(:,end,:) = xvalue(:,end-2,:);
xvalue(:,end-1,:) = xvalue(:,end-2,:);

yvalue = zeros(Nx + 4, Ny + 4, Tx);
yvalue(3:end-2, 3:end-2, :) = im(:,:,:,2);
yvalue(1,3:end-2,:) = im(1,:,:,2);
yvalue(2,3:end-2,:) = im(1,:,:,2);
yvalue(end,3:end-2,:) = im(end,:,:,2);
yvalue(end-1,3:end-2,:) = im(end,:,:,2);
yvalue(:,1,:) = yvalue(:,3,:);
yvalue(:,2,:) = yvalue(:,3,:);
yvalue(:,end,:) = yvalue(:,end-2,:);
yvalue(:,end-1,:) = yvalue(:,end-2,:);
