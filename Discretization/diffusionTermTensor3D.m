function M = diffusionTermTensor3D(D)
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
Nz = D.domain.dims(3);
Nt = (Nx+2)*(Ny+2)*(Nz+2);
G=reshape((1:Nt), Nx+2, Ny+2, Nz+2);
DX = repmat(D.domain.cellsize.x, 1, Ny, Nz);
DY = repmat(D.domain.cellsize.y', Nx, 1, Nz);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = D.domain.cellsize.z;
DZ=repmat(DZ, Nx, Ny, 1);
dx = 0.5*(DX(1:end-1,:,:)+DX(2:end,:,:));
dy = 0.5*(DY(:,1:end-1,:)+DY(:,2:end,:));
dz = 0.5*(DZ(:,:,1:end-1)+DZ(:,:,2:end));

% define the vectors to stores the sparse matrix data
iix = zeros(3*Nt,1);  jjx = zeros(3*Nt,1);  sx = zeros(3*Nt,1);
iiy = zeros(3*Nt,1);  jjy = zeros(3*Nt,1);  sy = zeros(3*Nt,1);
iiz = zeros(3*Nt,1);  jjz = zeros(3*Nt,1);  sz = zeros(3*Nt,1);
iixy = zeros(4*Nt,1); jjxy = zeros(4*Nt,1); sxy = zeros(4*Nt,1);
iixz = zeros(4*Nt,1); jjxz = zeros(4*Nt,1); sxz = zeros(4*Nt,1);
iiyz = zeros(4*Nt,1); jjyz = zeros(4*Nt,1); syz = zeros(4*Nt,1);
mnx = Nx*Ny*Nz;	      mny = Nx*Ny*Nz;       mnz = Nx*Ny*Nz;

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
Dx = D.xvalue;
Dy = D.yvalue;
Dz = D.zvalue;

% reassign the east, west, north, and south velocity vectors for the
% code readability
rx_1 = 1:Nx;     ry_1 = 1:Ny;     rz_1 = 1:Nz;
rx0 = 2:Nx+1;    ry0 = 2:Ny+1;    rz0 = 2:Nz+1;
rx1 = 3:Nx+2;    ry1 = 3:Ny+2;    rz1 = 3:Nz+2;
De = Dx(rx1,ry0,rz0,1)./(dx(rx0,:,:).*DX(rx0,:,:));
Dw = Dx(rx0,ry0,rz0,1)./(dx(rx0,:,:).*DX(rx0,:,:));
Dn = Dy(rx0,ry1,rz0,2)./(dy(:,ry0,:).*DY(:,ry0,:));
Ds = Dy(rx0,ry0,rz0,2)./(dy(:,ry0,:).*DY(:,ry0,:));
Df = Dz(rx0,ry0,rz1,3)./(dz(:,:,rz0).*DZ(:,:,rz0));
Db = Dz(rx0,ry0,rz0,3)./(dz(:,:,rz0).*DZ(:,:,rz0));
Dex2 = Dx(rx1,ry0,rz0,2)./(dx(rx0,:,:).*DX(rx0,:,:));
Dwx2 = Dx(rx_1,ry0,rz0,2)./(dx(rx0,:,:).*DX(rx0,:,:));
Dny2 = Dy(rx0,ry1,rz0,1)./(dy(:,ry0,:).*DY(:,ry0,:));
Dsy2 = Dy(rx0,ry_1,rz0,1)./(dy(:,ry0,:).*DY(:,ry0,:));
Dfz2 = Dz(rx0,ry0,rz1,1)./(dz(:,:,rz0).*DZ(:,:,rz0));
Dbz2 = Dz(rx0,ry0,rz_1,1)./(dz(:,:,rz0).*DZ(:,:,rz0));
Dfx2 = Dx(rx1,ry0,rz0,3)./(dx(rx0,:,:).*DX(rx0,:,:));
Dbx2 = Dx(rx_1,ry0,rz0,3)./(dx(rx0,:,:).*DX(rx0,:,:));
Dfny2 = Dy(rx0,ry1,rz0,3)./(dy(:,ry0,:).*DY(:,ry0,:));
Dbny2 = Dy(rx0,ry_1,rz0,3)./(dy(:,ry0,:).*DY(:,ry0,:));
Dfsz2 = Dz(rx0,ry0,rz1,2)./(dz(:,:,rz0).*DZ(:,:,rz0));
Dbsz2 = Dz(rx0,ry0,rz_1,2)./(dz(:,:,rz0).*DZ(:,:,rz0));

% calculate the coefficients for the internal cells
AE = reshape(De,mnx,1);
AW = reshape(Dw,mnx,1);
AN = reshape(Dn,mny,1);
AS = reshape(Ds,mny,1);
AF = reshape(Df,mnz,1);
AB = reshape(Db,mnz,1);
APx = reshape(-(De+Dw),mnx,1);
APy = reshape(-(Dn+Ds),mny,1);
APz = reshape(-(Df+Db),mnz,1);
ASE = reshape(-(Dex2+Dsy2)/4,mnx,1); % crossed xy
ANE = reshape(+(Dex2+Dny2)/4,mnx,1);
ANW = reshape(-(Dwx2+Dny2)/4,mnx,1);
ASW = reshape(+(Dwx2+Dsy2)/4,mnx,1);
ABE = reshape(-(Dfx2+Dbz2)/4,mnx,1); % crossed xz
AFE = reshape(+(Dfx2+Dfz2)/4,mnx,1);
AFW = reshape(-(Dbx2+Dfz2)/4,mnx,1);
ABW = reshape(+(Dbx2+Dbz2)/4,mnx,1);
ABN = reshape(-(Dfny2+Dbsz2)/4,mnx,1); % crossed yz
AFN = reshape(+(Dfny2+Dfsz2)/4,mnx,1);
AFS = reshape(-(Dbny2+Dfsz2)/4,mnx,1);
ABS = reshape(+(Dbny2+Dbsz2)/4,mnx,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnx,1); % main diagonal x
iix(1:3*mnx) = repmat(rowx_index,3,1);
rowy_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mny,1); % main diagonal y
iiy(1:3*mny) = repmat(rowy_index,3,1);
rowz_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnz,1); % main diagonal z
iiz(1:3*mnz) = repmat(rowz_index,3,1);
jjx(1:3*mnx) = [reshape(G(1:Nx,2:Ny+1,2:Nz+1),mnx,1); reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnx,1); reshape(G(3:Nx+2,2:Ny+1,2:Nz+1),mnx,1)];
jjy(1:3*mny) = [reshape(G(2:Nx+1,1:Ny,2:Nz+1),mny,1); reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mny,1); reshape(G(2:Nx+1,3:Ny+2,2:Nz+1),mny,1)];
jjz(1:3*mnz) = [reshape(G(2:Nx+1,2:Ny+1,1:Nz),mnz,1); reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnz,1); reshape(G(2:Nx+1,2:Ny+1,3:Nz+2),mnz,1)];
sx(1:3*mnx) = [AW; APx; AE];
sy(1:3*mny) = [AS; APy; AN];
sz(1:3*mnz) = [AB; APz; AF];

iixy(1:4*mny) = repmat(rowx_index,4,1);
jjxy(1:4*mnx) = [reshape(G(3:Nx+2,1:Ny,2:Nz+1),mnx,1); reshape(G(1:Nx,3:Ny+2,2:Nz+1),mnx,1);...
                 reshape(G(3:Nx+2,3:Ny+2,2:Nz+1),mnx,1); reshape(G(1:Nx,1:Ny,2:Nz+1),mnx,1);];
sxy(1:4*mnx) = [ASE; ANW; ANE; ASW];
iixz(1:4*mny) = repmat(rowx_index,4,1);
jjxz(1:4*mnx) = [reshape(G(3:Nx+2,2:Ny+1,1:Nz),mnx,1); reshape(G(1:Nx,2:Ny+1,3:Nz+2),mnx,1);...
                 reshape(G(3:Nx+2,2:Ny+1,3:Nz+2),mnx,1); reshape(G(1:Nx,2:Ny+1,1:Nz),mnx,1);];
sxz(1:4*mnx) = [ABE; AFW; AFE; ABW];
iiyz(1:4*mny) = repmat(rowy_index,4,1);
jjyz(1:4*mnx) = [reshape(G(2:Nx+1,3:Ny+2,1:Nz),mnx,1); reshape(G(2:Nx+1,1:Ny,3:Nz+2),mnx,1);...
                 reshape(G(2:Nx+1,3:Ny+2,3:Nz+2),mnx,1); reshape(G(2:Nx+1,1:Ny,1:Nz),mnx,1);];
syz(1:4*mnx) = [ABN; AFS; AFN; ABS];

% build the sparse matrix
kx = 3*mnx;
ky = 3*mny;
kz = 3*mnz;
kxy = 4*mnx;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nt, Nt);
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), Nt, Nt);
Mz = sparse(iiz(1:kz), jjz(1:kz), sz(1:kz), Nt, Nt);
Mxy = sparse(iixy(1:kxy), jjxy(1:kxy), sxy(1:kxy), Nt, Nt);
Mxz = sparse(iixz(1:kxy), jjxz(1:kxy), sxz(1:kxy), Nt, Nt);
Myz = sparse(iiyz(1:kxy), jjyz(1:kxy), syz(1:kxy), Nt, Nt);
M = Mx + My + Mz + Mxy + Mxz + Myz;
