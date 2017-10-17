function [M, Mx, My, Mz] = diffusionTerm3D(D)
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
Nx = D.domain.dims(1)+2;
Ny = D.domain.dims(2)+2;
Nz = D.domain.dims(3)+2;
Nt = (Nx+2)*(Ny+2)*(Nz+2);
G=reshape((1:Nt), Nx+2, Ny+2, Nz+2);
cellsize_x = [D.domain.cellsize.x(1); D.domain.cellsize.x;...
              D.domain.cellsize.x(end)];
cellsize_y = [D.domain.cellsize.y(1); D.domain.cellsize.y;...
              D.domain.cellsize.y(end)];
cellsize_z = [D.domain.cellsize.z(1); D.domain.cellsize.z;...
              D.domain.cellsize.z(end)];
DX = repmat(cellsize_x, 1, Ny+2, Nz+2);
DY = repmat(cellsize_y', Nx+2, 1, Nz+2);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = cellsize_z;
DZ=repmat(DZ, Nx+2, Ny+2, 1);

% define the vectors to stores the sparse matrix data
a = 5; %a=3
iix = zeros(a*Nt,1);
jjx = zeros(a*Nt,1);
sx = zeros(a*Nt,1);
iiy = zeros(a*Nt,1);
jjy = zeros(a*Nt,1);
sy = zeros(a*Nt,1);
iiz = zeros(a*Nt,1);
jjz = zeros(a*Nt,1);
sz = zeros(a*Nt,1);
mnx = (Nx-2)*(Ny-2)*(Nz-2);
mny = (Nx-2)*(Ny-2)*(Nz-2);
mnz = (Nx-2)*(Ny-2)*(Nz-2);


% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
Dx = D.xvalue;
Dy = D.yvalue;
Dz = D.zvalue;

% reassign the east, west, north, and south velocity vectors for the
% code readability
% IMPORTANT: I need to extend data in extremes twice instead of once
% for x 11 is the diagonal. D1 = (D11, D12, D13) es xvalue
% for y 22 is the diagonal. D2 = (D21, D22, D23) es yvalue
% for z 33 is the diagonal. D3 = (D31, D32, D33) es zvalue

% terms \phi(x,y,y)
Dc_xx = Dx(3:Nx,3:Ny,3:Nz,1)./(DX(3:Nx,3:Ny,3:Nz).^2);
Dc_yy = Dy(3:Nx,3:Ny,3:Nz,2)./(DY(3:Nx,3:Ny,3:Nz).^2);
Dc_zz = Dz(3:Nx,3:Ny,3:Nz,3)./(DZ(3:Nx,3:Ny,3:Nz).^2);
Dc_xy = Dx(3:Nx,3:Ny,3:Nz,2)./(DX(3:Nx,3:Ny,3:Nz).*DY(3:Nx,3:Ny,3:Nz));
%Dc_yx = Dy(3:Nx,3:Ny,3:Nz,1)./(DX(3:Nx,3:Ny,3:Nz).*DY(3:Nx,3:Ny,3:Nz));
Dc_xz = Dx(3:Nx,3:Ny,3:Nz,3)./(DX(3:Nx,3:Ny,3:Nz).*DZ(3:Nx,3:Ny,3:Nz));
%Dc_zx = Dz(3:Nx,3:Ny,3:Nz,1)./(DX(3:Nx,3:Ny,3:Nz).*DZ(3:Nx,3:Ny,3:Nz));
Dc_yz = Dy(3:Nx,3:Ny,3:Nz,3)./(DY(3:Nx,3:Ny,3:Nz).*DZ(3:Nx,3:Ny,3:Nz));
%Dc_zy = Dz(3:Nx,3:Ny,3:Nz,2)./(DY(3:Nx,3:Ny,3:Nz).*DZ(3:Nx,3:Ny,3:Nz));

% terms \phi(x+2dx,y,z), \phi(x-2dx,y,z), \phi(x,y+2dy,z),
%       \phi(x,y-2dy,z), \phi(x,y,z+2dz), \phi(x,y,z-2dz)
De2_xx = Dx(5:Nx+2,3:Ny,3:Nz,1)./(DX(5:Nx+2,3:Ny,3:Nz).^2);
Dw2_xx = Dx(1:Nx-2,3:Ny,3:Nz,1)./(DX(1:Nx-2,3:Ny,3:Nz).^2);
Dn2_yy = Dy(3:Nx,5:Ny+2,3:Nz,2)./(DY(3:Nx,5:Ny+2,3:Nz).^2);
Ds2_yy = Dy(3:Nx,1:Ny-2,3:Nz,2)./(DY(3:Nx,1:Ny-2,3:Nz).^2);
Df2_zz = Dz(3:Nx,3:Ny,5:Nz+2,3)./(DZ(3:Nx,3:Ny,5:Nz+2).^2);
Db2_zz = Dz(3:Nx,3:Ny,1:Nz-2,3)./(DZ(3:Nx,3:Ny,1:Nz-2).^2);

% \phi(x+dx,y,z), \phi(x-dx,y,z) 
De_xx = Dx(4:Nx+1,3:Ny,3:Nz,1)./(DX(4:Nx+1,3:Ny,3:Nz).^2);
De_xy = Dx(4:Nx+1,3:Ny,3:Nz,2)./(DX(4:Nx+1,3:Ny,3:Nz).*DY(4:Nx+1,3:Ny,3:Nz));
De_yx = Dy(4:Nx+1,3:Ny,3:Nz,1)./(DX(4:Nx+1,3:Ny,3:Nz).*DY(4:Nx+1,3:Ny,3:Nz));
De_xz = Dx(4:Nx+1,3:Ny,3:Nz,3)./(DX(4:Nx+1,3:Ny,3:Nz).*DZ(4:Nx+1,3:Ny,3:Nz));
De_zx = Dz(4:Nx+1,3:Ny,3:Nz,1)./(DX(4:Nx+1,3:Ny,3:Nz).*DZ(4:Nx+1,3:Ny,3:Nz));
Dw_xx = Dx(2:Nx-1,3:Ny,3:Nz,1)./(DX(2:Nx-1,3:Ny,3:Nz).^2);
Dw_xy = Dx(2:Nx-1,3:Ny,3:Nz,2)./(DX(2:Nx-1,3:Ny,3:Nz).*DY(2:Nx-1,3:Ny,3:Nz));
Dw_yx = Dy(2:Nx-1,3:Ny,3:Nz,1)./(DX(2:Nx-1,3:Ny,3:Nz).*DY(2:Nx-1,3:Ny,3:Nz));
Dw_xz = Dx(2:Nx-1,3:Ny,3:Nz,3)./(DX(2:Nx-1,3:Ny,3:Nz).*DZ(2:Nx-1,3:Ny,3:Nz));
Dw_zx = Dz(2:Nx-1,3:Ny,3:Nz,2)./(DX(2:Nx-1,3:Ny,3:Nz).*DZ(2:Nx-1,3:Ny,3:Nz));
% \phi(x,y+dy,z), \phi(x,y-dy,z)
Dn_yy = Dy(3:Nx,4:Ny+1,3:Nz,2)./(DY(3:Nx,4:Ny+1,3:Nz).^2);
Dn_yx = Dy(3:Nx,4:Ny+1,3:Nz,1)./(DY(3:Nx,4:Ny+1,3:Nz).*DX(3:Nx,4:Ny+1,3:Nz));
Dn_xy = Dx(3:Nx,4:Ny+1,3:Nz,2)./(DY(3:Nx,4:Ny+1,3:Nz).*DX(3:Nx,4:Ny+1,3:Nz));
Dn_yz = Dy(3:Nx,4:Ny+1,3:Nz,3)./(DY(3:Nx,4:Ny+1,3:Nz).*DZ(3:Nx,4:Ny+1,3:Nz));
Dn_zy = Dz(3:Nx,4:Ny+1,3:Nz,2)./(DY(3:Nx,4:Ny+1,3:Nz).*DZ(3:Nx,4:Ny+1,3:Nz));
Ds_yy = Dy(3:Nx,2:Ny-1,3:Nz,2)./(DY(3:Nx,2:Ny-1,3:Nz).^2);
Ds_yx = Dy(3:Nx,2:Ny-1,3:Nz,1)./(DY(3:Nx,2:Ny-1,3:Nz).*DX(3:Nx,2:Ny-1,3:Nz));
Ds_xy = Dx(3:Nx,2:Ny-1,3:Nz,2)./(DY(3:Nx,2:Ny-1,3:Nz).*DX(3:Nx,2:Ny-1,3:Nz));
Ds_yz = Dy(3:Nx,2:Ny-1,3:Nz,3)./(DY(3:Nx,2:Ny-1,3:Nz).*DZ(3:Nx,2:Ny-1,3:Nz));
Ds_zy = Dz(3:Nx,2:Ny-1,3:Nz,2)./(DY(3:Nx,2:Ny-1,3:Nz).*DZ(3:Nx,2:Ny-1,3:Nz));
% \phi(x,y,z+dz), \phi(x,y,z-dz)
Df_zz = Dz(3:Nx,3:Ny,4:Nz+1,3)./(DZ(3:Nx,3:Ny,4:Nz+1).^2);
Df_zx = Dz(3:Nx,3:Ny,4:Nz+1,1)./(DZ(3:Nx,3:Ny,4:Nz+1).*DX(3:Nx,3:Ny,4:Nz+1));
Df_xz = Dx(3:Nx,3:Ny,4:Nz+1,3)./(DZ(3:Nx,3:Ny,4:Nz+1).*DX(3:Nx,3:Ny,4:Nz+1));
Df_yz = Dy(3:Nx,3:Ny,4:Nz+1,3)./(DZ(3:Nx,3:Ny,4:Nz+1).*DY(3:Nx,3:Ny,4:Nz+1));
Df_zy = Dz(3:Nx,3:Ny,4:Nz+1,2)./(DZ(3:Nx,3:Ny,4:Nz+1).*DY(3:Nx,3:Ny,4:Nz+1));
Db_zz = Dy(3:Nx,3:Ny,2:Nz-1,2)./(DZ(3:Nx,3:Ny,2:Nz-1).^2);
Db_zx = Dy(3:Nx,3:Ny,2:Nz-1,1)./(DZ(3:Nx,3:Ny,2:Nz-1).*DX(3:Nx,3:Ny,2:Nz-1));
Db_xz = Dx(3:Nx,3:Ny,2:Nz-1,2)./(DZ(3:Nx,3:Ny,2:Nz-1).*DX(3:Nx,3:Ny,2:Nz-1));
Db_yz = Dy(3:Nx,3:Ny,2:Nz-1,3)./(DZ(3:Nx,3:Ny,2:Nz-1).*DY(3:Nx,3:Ny,2:Nz-1));
Db_zy = Dz(3:Nx,3:Ny,2:Nz-1,2)./(DZ(3:Nx,3:Ny,2:Nz-1).*DY(3:Nx,3:Ny,2:Nz-1));

% \phi(x+dx,y+dy,z),\phi(x-dx,y+dy,z),\phi(x+dx,y-dy,z),\phi(x-dx,y-dy,z)
Dws_xy = Dx(2:Nx-1,2:Ny-1,3:Nz,2)./(DX(2:Nx-1,2:Ny-1,3:Nz).*DY(2:Nx-1,2:Ny-1,3:Nz));
Dws_yx = Dy(2:Nx-1,2:Ny-1,3:Nz,1)./(DX(2:Nx-1,2:Ny-1,3:Nz).*DY(2:Nx-1,2:Ny-1,3:Nz));
Dwn_xy = Dx(2:Nx-1,4:Ny+1,3:Nz,2)./(DX(2:Nx-1,4:Ny+1,3:Nz).*DY(2:Nx-1,4:Ny+1,3:Nz));
Dwn_yx = Dy(2:Nx-1,4:Ny+1,3:Nz,1)./(DX(2:Nx-1,4:Ny+1,3:Nz).*DY(2:Nx-1,4:Ny+1,3:Nz));
Des_xy = Dx(4:Nx+1,2:Ny-1,3:Nz,2)./(DX(4:Nx+1,2:Ny-1,3:Nz).*DY(4:Nx+1,2:Ny-1,3:Nz));
Des_yx = Dy(4:Nx+1,2:Ny-1,3:Nz,1)./(DX(4:Nx+1,2:Ny-1,3:Nz).*DY(4:Nx+1,2:Ny-1,3:Nz));
Den_xy = Dx(4:Nx+1,4:Ny+1,3:Nz,2)./(DX(4:Nx+1,4:Ny+1,3:Nz).*DY(4:Nx+1,4:Ny+1,3:Nz));
Den_yx = Dy(4:Nx+1,4:Ny+1,3:Nz,1)./(DX(4:Nx+1,4:Ny+1,3:Nz).*DY(4:Nx+1,4:Ny+1,3:Nz));
% \phi(x+dx,y,z+dz),\phi(x-dx,y,z+dz),\phi(x+dx,y,z-dz),\phi(x-dx,y,z-dz)
Dwb_xz = Dx(2:Nx-1,3:Ny,2:Nz-1,3)./(DX(2:Nx-1,3:Ny,2:Nz-1).*DZ(2:Nx-1,3:Ny,2:Nz-1));
Dwb_zx = Dz(2:Nx-1,3:Ny,2:Nz-1,1)./(DX(2:Nx-1,3:Ny,2:Nz-1).*DZ(2:Nx-1,3:Ny,2:Nz-1));
Dwf_xz = Dx(2:Nx-1,3:Ny,4:Nz+1,3)./(DX(2:Nx-1,3:Ny,4:Nz+1).*DZ(2:Nx-1,3:Ny,4:Nz+1));
Dwf_zx = Dz(2:Nx-1,3:Ny,4:Nz+1,1)./(DX(2:Nx-1,3:Ny,4:Nz+1).*DZ(2:Nx-1,3:Ny,4:Nz+1));
Deb_xz = Dx(4:Nx+1,3:Ny,2:Nz-1,3)./(DX(4:Nx+1,3:Ny,2:Nz-1).*DZ(4:Nx+1,3:Ny,2:Nz-1));
Deb_zx = Dz(4:Nx+1,3:Ny,2:Nz-1,1)./(DX(4:Nx+1,3:Ny,2:Nz-1).*DZ(4:Nx+1,3:Ny,2:Nz-1));
Def_xz = Dx(4:Nx+1,3:Ny,4:Nz+1,3)./(DX(4:Nx+1,3:Ny,4:Nz+1).*DZ(4:Nx+1,3:Ny,4:Nz+1));
Def_zx = Dz(4:Nx+1,3:Ny,4:Nz+1,1)./(DX(4:Nx+1,3:Ny,4:Nz+1).*DZ(4:Nx+1,3:Ny,4:Nz+1));
% \phi(x,y+dy,z+dz),\phi(x,y-dy,z+dz),\phi(x,y+dy,z-dz),\phi(x,y-dy,z-dz)
Dsb_yz = Dy(3:Nx,2:Ny-1,2:Nz-1,3)./(DY(3:Nx,2:Ny-1,2:Nz-1).*DZ(3:Nx,2:Ny-1,2:Nz-1));
Dsb_zy = Dz(3:Nx,2:Ny-1,2:Nz-1,2)./(DY(3:Nx,2:Ny-1,2:Nz-1).*DZ(3:Nx,2:Ny-1,2:Nz-1));
Dsf_yz = Dy(3:Nx,2:Ny-1,4:Nz+1,3)./(DY(3:Nx,2:Ny-1,4:Nz+1).*DZ(3:Nx,2:Ny-1,4:Nz+1));
Dsf_zy = Dz(3:Nx,2:Ny-1,4:Nz+1,2)./(DY(3:Nx,2:Ny-1,4:Nz+1).*DZ(3:Nx,2:Ny-1,4:Nz+1));
Dnb_yz = Dy(3:Nx,4:Ny+1,2:Nz-1,3)./(DY(3:Nx,4:Ny+1,2:Nz-1).*DZ(3:Nx,4:Ny+1,2:Nz-1));
Dnb_zy = Dz(3:Nx,4:Ny+1,2:Nz-1,2)./(DY(3:Nx,4:Ny+1,2:Nz-1).*DZ(3:Nx,4:Ny+1,2:Nz-1));
Dnf_yz = Dy(3:Nx,4:Ny+1,4:Nz+1,3)./(DY(3:Nx,4:Ny+1,4:Nz+1).*DZ(3:Nx,4:Ny+1,4:Nz+1));
Dnf_zy = Dz(3:Nx,4:Ny+1,4:Nz+1,2)./(DY(3:Nx,4:Ny+1,4:Nz+1).*DZ(3:Nx,4:Ny+1,4:Nz+1));

% calculate the coefficients for the internal cells
APx = reshape(  6 * Dc_xx + 4 * 2 * Dc_xy, mnx, 1);
APy = reshape(  6 * Dc_yy + 4 * 2 * Dc_xz, mny, 1);
APz = reshape(  6 * Dc_zz + 4 * 2 * Dc_yz, mnz, 1);
AE1 = reshape(- 4 * De_xx - 2 * (De_xy + De_yx + De_xz + De_zx), mnx, 1);
AW1 = reshape(- 4 * Dw_xx - 2 * (Dw_xy + Dw_yx + Dw_xz + Dw_zx), mnx, 1);
AN1 = reshape(- 4 * Dn_yy - 2 * (Dn_xy + Dn_yx + Dn_zy + Dn_yz), mny, 1);
AS1 = reshape(- 4 * Ds_yy - 2 * (Ds_xy + Ds_yx + Ds_zy + Ds_yz), mny, 1);
AF1 = reshape(- 4 * Df_zz - 2 * (Df_xz + Df_zx + Df_yz + Df_zy), mnz, 1);
AB1 = reshape(- 4 * Db_zz - 2 * (Db_xz + Db_zx + Db_yz + Db_zy), mnz, 1);
AE2 = reshape(De2_xx, mnx, 1);
AW2 = reshape(Dw2_xx, mnx, 1);
AN2 = reshape(Dn2_yy, mny, 1);
AS2 = reshape(Ds2_yy, mny, 1);
AF2 = reshape(Df2_zz, mnz, 1);
AB2 = reshape(Db2_zz, mnz, 1);
AWS = reshape(Dws_xy + Dws_yx, mnx, 1);
AWN = reshape(Dwn_xy + Dwn_yx, mnx, 1);
AES = reshape(Des_xy + Des_yx, mnx, 1);
AEN = reshape(Den_xy + Den_yx, mnx, 1);
AWB = reshape(Dwb_xz + Dwb_zx, mny, 1);
AWF = reshape(Dwf_xz + Dwf_zx, mny, 1);
AEB = reshape(Deb_xz + Deb_zx, mny, 1);
AEF = reshape(Def_xz + Def_zx, mny, 1);
ASB = reshape(Dsb_yz + Dsb_zy, mnz, 1);
ASF = reshape(Dsf_yz + Dsf_zy, mnz, 1);
ANB = reshape(Dnb_yz + Dnb_zy, mnz, 1);
ANF = reshape(Dnf_yz + Dnf_zy, mnz, 1);

% build the sparse matrix based on the numbering system
% main terms for x and y dimensions
rowx_index = reshape(G(3:Nx,3:Ny,3:Nz),mnx,1); % main diagonal x
iix(1:5*mnx) = repmat(rowx_index,5,1);
rowy_index = reshape(G(3:Nx,3:Ny,3:Nz),mny,1); % main diagonal y
iiy(1:5*mny) = repmat(rowy_index,5,1);
rowz_index = reshape(G(3:Nx,3:Ny,3:Nz),mnz,1); % main diagonal z
iiz(1:5*mnz) = repmat(rowz_index,5,1);
iixy(1:4*mnx) = repmat(rowx_index,4,1);
iixz(1:4*mnx) = repmat(rowx_index,4,1);
iiyz(1:4*mnx) = repmat(rowy_index,4,1);

jjx(1:5*mnx) = [reshape(G(1:Nx-2,3:Ny,3:Nz),mnx,1); ...
                reshape(G(2:Nx-1,3:Ny,3:Nz),mnx,1); ...
                reshape(G(3:Nx,3:Ny,3:Nz),mnx,1); ...
                reshape(G(4:Nx+1,3:Ny,3:Nz),mnx,1); ...
                reshape(G(5:Nx+2,3:Ny,3:Nz),mnx,1)];
jjy(1:5*mny) = [reshape(G(3:Nx,1:Ny-2,3:Nz),mny,1); ...
                reshape(G(3:Nx,2:Ny-1,3:Nz),mny,1); ...
                reshape(G(3:Nx,3:Ny,3:Nz),mny,1); ...
                reshape(G(3:Nx,4:Ny+1,3:Nz),mny,1); ...
                reshape(G(3:Nx,5:Ny+2,3:Nz),mny,1)];
jjz(1:5*mnz) = [reshape(G(3:Nx,3:Ny,1:Nz-2),mnz,1); ...
                reshape(G(3:Nx,3:Ny,2:Nz-1),mnz,1); ...
                reshape(G(3:Nx,3:Ny,3:Nz),mnz,1); ...
                reshape(G(3:Nx,3:Ny,4:Nz+1),mnz,1); ...
                reshape(G(3:Nx,3:Ny,5:Nz+2),mnz,1)];

jjxy(1:4*mnx) = [reshape(G(2:Nx-1,2:Ny-1,3:Nz),mnx,1); ...
                reshape(G(2:Nx-1,4:Ny+1,3:Nz),mnx,1); ...
                reshape(G(4:Nx+1,2:Ny-1,3:Nz),mnx,1); ...
                reshape(G(4:Nx+1,4:Ny+1,3:Nz),mnx,1)];
jjxz(1:4*mnx) = [reshape(G(2:Nx-1,3:Ny,2:Nz-1),mnx,1); ...
                reshape(G(2:Nx-1,3:Ny,4:Nz+1),mnx,1); ...
                reshape(G(4:Nx+1,3:Ny,2:Nz-1),mnx,1); ...
                reshape(G(4:Nx+1,3:Ny,4:Nz+1),mnx,1)];
jjyz(1:4*mnx) = [reshape(G(3:Nx,2:Ny-1,2:Nz-1),mnz,1); ...
                reshape(G(3:Nx,2:Ny-1,4:Nz+1),mnz,1); ...
                reshape(G(3:Nx,4:Ny+1,2:Nz-1),mnz,1); ...
                reshape(G(3:Nx,4:Ny+1,4:Nz+1),mnz,1)];

sx(1:5*mnx) = [AW2; AW1; APx; AE1; AE2];
sy(1:5*mny) = [AS2; AS1; APy; AN1; AN2];
sz(1:5*mnz) = [AB2; AB1; APz; AF1; AF2];
sxy(1:4*mnx) = [AWS; AWN; AES; AEN];
sxz(1:4*mnx) = [AWB; AWF; AEB; AEF];
syz(1:4*mnz) = [ASB; ASF; ANB; ANF];

% build the sparse matrix
kx = 5*mnx; 
ky = 5*mny;
kz = 5*mnz;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nt, Nt);
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), Nt, Nt);
Mz = sparse(iiz(1:kz), jjz(1:kz), sz(1:kz), Nt, Nt);
kxy = 4*mnx;
kxz = 4*mnx;
kyz = 4*mnx;
Mxy = sparse(iixy(1:kxy), jjxy(1:kxy), sxy(1:kxy), Nt, Nt);
Mxz = sparse(iixz(1:kxz), jjxz(1:kxz), sxz(1:kxz), Nt, Nt);
Myz = sparse(iiyz(1:kyz), jjyz(1:kyz), syz(1:kyz), Nt, Nt);

M_aux = full(Mx + My + Mz + Mxy + Mxz + Myz);
indexes = reshape(G(2:end-1,2:end-1,2:end-1), Nx*Ny*Nz, 1);
M = M_aux(indexes, indexes);

