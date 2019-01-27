function phiFaceAverage = arithmeticMeanTensor(phi)
% This function gets the value of the field variable phi defined
% over the MeshStructure and calculates the arithmetic average on
% the cell faces, for a uniform mesh.
%
% Author: Aina Frau-Pascual, adapting code from arithmeticMean.m
% for its use with tensors


% extract data from the mesh structure

d = phi.domain.dimension;
if (d == 2)
    Nx = phi.domain.dims(1);
    Ny = phi.domain.dims(2);
    dx = repmat(phi.domain.cellsize.x, 1, Ny);
    dy = repmat(phi.domain.cellsize.y', Nx, 1);
	% if I dont average cross terms
    xvalue=phi.xvalue;
    yvalue=phi.yvalue;
    % diagonal terms
    xvalue(2:end,2:end-1,1)=(dx(1:end-1,:).*phi.xvalue(1:end-1,2:end-1,1)+...
        dx(2:end,:).*phi.xvalue(2:end,2:end-1,1))./(dx(2:end,:)+dx(1:end-1,:));
    yvalue(2:end-1,2:end,2)=(dy(:,1:end-1).*phi.yvalue(2:end-1,1:end-1,2)+...
        dy(:,2:end).*phi.yvalue(2:end-1,2:end,2))./(dy(:,2:end)+dy(:,1:end-1));
    xvalue(2:end,1,1)=xvalue(2:end,2,1);
    xvalue(2:end,end,1)=xvalue(2:end,end-1,1);
    yvalue(1,2:end,2)=yvalue(2,2:end,2);
    yvalue(end,2:end,2)=yvalue(end-1,2:end,2);
    zvalue=[];
elseif (d == 3)
    Nx = phi.domain.dims(1);
    Ny = phi.domain.dims(2);
    Nz = phi.domain.dims(3);
    dx = repmat(phi.domain.cellsize.x, 1, Ny, Nz);
    dy = repmat(phi.domain.cellsize.y', Nx, 1, Nz);
    DZ = zeros(1,1,Nz+2);
    DZ(1,1,:) = phi.domain.cellsize.z;
    dz=repmat(DZ, Nx, Ny, 1);
    % if I dont average cross terms
    xvalue=phi.xvalue;
    yvalue=phi.yvalue;
    zvalue=phi.zvalue;
    % diagonal terms
    xvalue(2:end,2:end-1,2:end-1,1)=(dx(1:end-1,:,:).*phi.xvalue(1:end-1,2:end-1,2:end-1,1)+...
        dx(2:end,:,:).*phi.xvalue(2:end,2:end-1,2:end-1,1))./(dx(2:end,:,:)+dx(1:end-1,:,:));
    yvalue(2:end-1,2:end,2:end-1,2)=(dy(:,1:end-1,:).*phi.yvalue(2:end-1,1:end-1,2:end-1,2)+...
        dy(:,2:end,:).*phi.yvalue(2:end-1,2:end,2:end-1,2))./(dy(:,1:end-1,:)+dy(:,2:end,:));
    zvalue(2:end-1,2:end-1,2:end,3)=(dz(:,:,1:end-1).*phi.zvalue(2:end-1,2:end-1,1:end-1,3)+...
        dz(:,:,2:end).*phi.zvalue(2:end-1,2:end-1,2:end,3))./(dz(:,:,1:end-1)+dz(:,:,2:end));
    xvalue(2:end,1,:,1)=xvalue(2:end,2,:,1);
    xvalue(2:end,end,:,1)=xvalue(2:end,end-1,:,1);
    xvalue(2:end,:,1,1)=xvalue(2:end,:,2,1);
    xvalue(2:end,:,end,1)=xvalue(2:end,:,end-1,1);
    yvalue(1,2:end,:,2)=yvalue(2,2:end,:,2);
    yvalue(end,2:end,:,2)=yvalue(end-1,2:end,:,2);
    yvalue(:,2:end,1,2)=yvalue(:,2:end,2,2);
    yvalue(:,2:end,end,2)=yvalue(:,2:end,end-1,2);
    zvalue(1,:,2:end,3)=zvalue(2,:,2:end,3);
    zvalue(end,:,2:end,3)=zvalue(end-1,:,2:end,3);
    zvalue(:,1,2:end,3)=zvalue(:,2,2:end,3);
    zvalue(:,end,2:end,3)=zvalue(:,end-1,2:end,3);
end
phiFaceAverage=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
