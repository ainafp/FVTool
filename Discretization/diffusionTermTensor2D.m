function M = diffusionTermTensor2D(D)
% This function uses the central difference scheme to discretize a 2D
% diffusion term in the form \grad . (D \grad \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%
%
% PARAMETERS:
%
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% extract data from the mesh structure
Nx = D.domain.dims(1);
Ny = D.domain.dims(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);

DX = repmat(D.domain.cellsize.x, 1, Ny);
DY = repmat(D.domain.cellsize.y', Nx, 1);
dx = 0.5*(DX(1:end-1,:)+DX(2:end,:));
dy = 0.5*(DY(:,1:end-1)+DY(:,2:end));

% define the vectors to store the sparse matrix data
iix = zeros(3*(Nx+2)*(Ny+2),1);	iiy = zeros(3*(Nx+2)*(Ny+2),1);
jjx = zeros(3*(Nx+2)*(Ny+2),1);	jjy = zeros(3*(Nx+2)*(Ny+2),1);
sx = zeros(3*(Nx+2)*(Ny+2),1);	sy = zeros(3*(Nx+2)*(Ny+2),1);
mnx = Nx*Ny;	mny = Nx*Ny;

% reassign the east, west, north, and south velocity vectors for the
% code readability
rx_1 = 1:Nx;     ry_1 = 1:Ny;
rx0 = 2:Nx+1;    ry0 = 2:Ny+1;
rx1 = 3:Nx+2;    ry1 = 3:Ny+2;
De = D.xvalue(rx1,ry0,1)./(dx(rx0,:).*DX(rx0,:));
Dw = D.xvalue(rx0,ry0,1)./(dx(rx0,:).*DX(rx0,:));
Dn = D.yvalue(rx0,ry1,2)./(dy(:,ry0).*DY(:,ry0));
Ds = D.yvalue(rx0,ry0,2)./(dy(:,ry0).*DY(:,ry0));
Dex2 = D.xvalue(rx1,ry0,2)./(dx(rx0,:).*DX(rx0,:));
Dwx2 = D.xvalue(rx_1,ry0,2)./(dx(rx0,:).*DX(rx0,:));
Dny2 = D.yvalue(rx0,ry1,1)./(dy(:,ry0).*DY(:,ry0));
Dsy2 = D.yvalue(rx0,ry_1,1)./(dy(:,ry0).*DY(:,ry0));

% calculate the coefficients for the internal cells
AE = reshape(De,mnx,1);
AS = reshape(Ds,mny,1);
AN = reshape(Dn,mny,1);
AW = reshape(Dw,mnx,1);
APx = reshape(-(De+Dw),mnx,1);
APy = reshape(-(Dn+Ds),mny,1);
ASE = reshape(-(Dex2+Dsy2)/4,mnx,1);
ANE = reshape(+(Dex2+Dny2)/4,mnx,1);
ANW = reshape(-(Dwx2+Dny2)/4,mnx,1);
ASW = reshape(+(Dwx2+Dsy2)/4,mnx,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1,2:Ny+1),mnx,1); % main diagonal x
iix(1:3*mnx) = repmat(rowx_index,3,1);
rowy_index = reshape(G(2:Nx+1,2:Ny+1),mny,1); % main diagonal y
iiy(1:3*mny) = repmat(rowy_index,3,1);
jjx(1:3*mnx) = [reshape(G(1:Nx,2:Ny+1),mnx,1); reshape(G(2:Nx+1,2:Ny+1),mnx,1); reshape(G(3:Nx+2,2:Ny+1),mnx,1)];
jjy(1:3*mny) = [reshape(G(2:Nx+1,1:Ny),mny,1); reshape(G(2:Nx+1,2:Ny+1),mny,1); reshape(G(2:Nx+1,3:Ny+2),mny,1)];
iixy(1:4*mny) = repmat(rowx_index,4,1);
jjxy(1:4*mnx) = [reshape(G(3:Nx+2,1:Ny),mnx,1); reshape(G(1:Nx,3:Ny+2),mnx,1);...
                 reshape(G(3:Nx+2,3:Ny+2),mnx,1); reshape(G(1:Nx,1:Ny),mnx,1);];
sx(1:3*mnx) = [AW; APx; AE];
sy(1:3*mny) = [AS; APy; AN];
sxy(1:4*mnx) = [ASE; ANW; ANE; ASW];

% build the sparse matrix
kx = 3*mnx;
ky = 3*mny;
kxy = 4*mnx;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));
Mxy = sparse(iixy(1:kxy), jjxy(1:kxy), sxy(1:kxy), (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));
M = Mx + My + Mxy;
