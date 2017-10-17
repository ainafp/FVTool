function [M, Mx, My] = diffusionTermTensor2D(D)
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
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);

DX = repmat(ones(Nx+2,1), 1, Ny+2); % size voxel?
DY = DX;

% define the vectors to store the sparse matrix data
iix = zeros(5*(Nx+2)*(Ny+2),1);	iiy = zeros(5*(Nx+2)*(Ny+2),1);
jjx = zeros(5*(Nx+2)*(Ny+2),1);	jjy = zeros(5*(Nx+2)*(Ny+2),1);
sx = zeros(5*(Nx+2)*(Ny+2),1);	sy = zeros(5*(Nx+2)*(Ny+2),1);
mnx = (Nx-2)*(Ny-2);	mny = (Nx-2)*(Ny-2);

% reassign the east, west, north, and south velocity vectors
% for x 11 is the diagonal. D1 = (D11, D12) es xvalue
% for y 22 is the diagonal. D2 = (D21, D22) es yvalue
% term \phi(x,y)
Dc_xx = D.xvalue(3:Nx,3:Ny,1)./(DX(3:Nx,3:Ny).^2);
Dc_xy = D.xvalue(3:Nx,3:Ny,2)./(DX(3:Nx,3:Ny).*DY(3:Nx,3:Ny));
Dc_yy = D.yvalue(3:Nx,3:Ny,2)./(DY(3:Nx,3:Ny).^2);
% terms \phi(x+2dx,y), \phi(x-2dx,y), \phi(x,y+2dy), \phi(x,y-2dy)
De2_xx = D.xvalue(5:Nx+2,3:Ny,1)./(DX(5:Nx+2,3:Ny).^2);
Dw2_xx = D.xvalue(1:Nx-2,3:Ny,1)./(DX(1:Nx-2,3:Ny).^2);
Dn2_yy = D.yvalue(3:Nx,5:Ny+2,2)./(DY(3:Nx,5:Ny+2).^2);
Ds2_yy = D.yvalue(3:Nx,1:Ny-2,2)./(DY(3:Nx,1:Ny-2).^2);
% terms \phi(x+dx,y), \phi(x-dx,y), \phi(x,y+dy), \phi(x,y-dy)
De_xx = D.xvalue(4:Nx+1,3:Ny,1)./(DX(4:Nx+1,3:Ny).^2);
De_xy = D.xvalue(4:Nx+1,3:Ny,2)./(DX(4:Nx+1,3:Ny).*DY(4:Nx+1,3:Ny));
Dw_xx = D.xvalue(2:Nx-1,3:Ny,1)./(DX(2:Nx-1,3:Ny).^2);
Dw_xy = D.xvalue(2:Nx-1,3:Ny,2)./(DX(2:Nx-1,3:Ny).*DY(2:Nx-1,3:Ny));
Dn_yy = D.yvalue(3:Nx,4:Ny+1,2)./(DY(3:Nx,4:Ny+1).^2);
Dn_xy = D.xvalue(3:Nx,4:Ny+1,2)./(DY(3:Nx,4:Ny+1).*DX(3:Nx,4:Ny+1));
Ds_yy = D.yvalue(3:Nx,2:Ny-1,2)./(DY(3:Nx,2:Ny-1).^2);
Ds_xy = D.xvalue(3:Nx,2:Ny-1,2)./(DY(3:Nx,2:Ny-1).*DX(3:Nx,2:Ny-1));
% terms \phi(x+dx,y+dy), \phi(x-dx,y+dy), \phi(x+dx,y-dy), \phi(x-dx,y-dy)
Dws_xy = D.xvalue(2:Nx-1,2:Ny-1,2)./(DX(2:Nx-1,2:Ny-1).*DY(2:Nx-1,2:Ny-1));
Dwn_xy = D.xvalue(2:Nx-1,4:Ny+1,2)./(DX(2:Nx-1,4:Ny+1).*DY(2:Nx-1,4:Ny+1));
Des_xy = D.xvalue(4:Nx+1,2:Ny-1,2)./(DX(4:Nx+1,2:Ny-1).*DY(4:Nx+1,2:Ny-1));
Den_xy = D.xvalue(4:Nx+1,4:Ny+1,2)./(DX(4:Nx+1,4:Ny+1).*DY(4:Nx+1,4:Ny+1));

% calculate the coefficients for the internal cells
APx = reshape(  6 * Dc_xx + 4 * Dc_xy, mnx, 1);
APy = reshape(  6 * Dc_yy + 4 * Dc_xy, mny, 1);
AE1 = reshape(- 4 * De_xx - 2 * 2 * De_xy, mnx, 1);
AW1 = reshape(- 4 * Dw_xx - 2 * 2 * Dw_xy, mnx, 1);
AN1 = reshape(- 4 * Dn_yy - 2 * 2 * Dn_xy, mny, 1);
AS1 = reshape(- 4 * Ds_yy - 2 * 2 * Ds_xy, mny, 1);
AE2 = reshape(De2_xx, mnx, 1);
AW2 = reshape(Dw2_xx, mnx, 1);
AN2 = reshape(Dn2_yy, mny, 1);
AS2 = reshape(Ds2_yy, mny, 1);
AWS = reshape(2 * Dws_xy, mnx, 1);
AWN = reshape(2 * Dwn_xy, mnx, 1);
AES = reshape(2 * Des_xy, mny, 1);
AEN = reshape(2 * Den_xy, mny, 1);

% build the sparse matrix based on the numbering system
% main terms for x and y dimensions
rowx_index = reshape(G(3:Nx,3:Ny),mnx,1); % main diagonal x
rowy_index = reshape(G(3:Nx,3:Ny),mny,1); % main diagonal y
iix(1:5*mnx) = repmat(rowx_index,5,1);
iiy(1:5*mny) = repmat(rowy_index,5,1);
iixy(1:4*mnx) = repmat(rowx_index,4,1);
jjx(1:5*mnx) = [reshape(G(1:Nx-2,3:Ny),mnx,1); ...
                reshape(G(2:Nx-1,3:Ny),mnx,1); ...
                reshape(G(3:Nx,3:Ny),mnx,1); ...
                reshape(G(4:Nx+1,3:Ny),mnx,1); ...
                reshape(G(5:Nx+2,3:Ny),mnx,1)];
jjy(1:5*mny) = [reshape(G(3:Nx,1:Ny-2),mny,1); ...
                reshape(G(3:Nx,2:Ny-1),mny,1); ...
                reshape(G(3:Nx,3:Ny),mny,1); ...
                reshape(G(3:Nx,4:Ny+1),mny,1); ...
                reshape(G(3:Nx,5:Ny+2),mny,1)];
jjxy(1:4*mnx) = [reshape(G(2:Nx-1,2:Ny-1),mnx,1); ...
                reshape(G(2:Nx-1,4:Ny+1),mnx,1); ...
                reshape(G(4:Nx+1,2:Ny-1),mnx,1); ...
                reshape(G(4:Nx+1,4:Ny+1),mnx,1)];

%sx(1:5*mnx-2*Nx+4) = [AW2; AW1; APx; AE1; AE2];
sx(1:5*mnx) = [AW2; AW1; APx; AE1; AE2];
sy(1:5*mny) = [AS2; AS1; APy; AN1; AN2];
sxy(1:4*mnx) = [AWS; AWN; AES; AEN];

% build the sparse matrix
kx = 5*mnx;
ky = 5*mny;
Nt = (Nx+2)*(Ny+2);
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nt, Nt);
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), Nt, Nt);
kxy = 4*mnx;
Mxy = sparse(iixy(1:kxy), jjxy(1:kxy), sxy(1:kxy), Nt, Nt);
M_aux = Mx + My + Mxy;
M = M_aux(G(2:end-1,2:end-1),G(2:end-1,2:end-1));
