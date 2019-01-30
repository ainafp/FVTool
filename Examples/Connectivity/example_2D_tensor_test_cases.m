% This example shows different cases of tensor maps forming shapes
% Author: Aina Frau-Pascual

imx = 19;
imy = imx;

%%
% Input image isotropic
im = zeros(imx, imy, 2, 2);
im(:,:,1,1) = 0.2;
im(:,:,2,2) = 0.2;
fa = sum(sum(im,4),3);

% Compute
[meshstruct, D, M, RHSbc] = computeDiffusionMatrix2D(im, [imx, imy]);
[potentials, mask_index] = computePotentials2D(D, meshstruct, M, RHSbc);
interactiveConductance2D(fa, [imx, imy], potentials, mask_index);

%%
% Input image cross
im = zeros(imx, imy, 2, 2);
tensor_base = [1, 0; 0, 0.2];
center = [round(imx / 2), round(imy / 2)];
for x=1:imx
    for y=1:imx
        im(x,y,:,:) = put_tensor(tensor_base,center,[x,y],0);
    end
end
im(10,10,:,:) = [1, 0; 0, 1];
fa = sum(sum(im,4),3);

% Compute
[meshstruct, D, M, RHSbc] = computeDiffusionMatrix2D(im, [imx, imy]);
[potentials, mask_index] = computePotentials2D(D, meshstruct, M, RHSbc);
interactiveConductance2D(fa, [imx, imy], potentials, mask_index);

%%
% Input image star
im = zeros(imx, imy, 2, 2);
tensor_base = [1, 0; 0, 0.2];
center = [round(imx / 2), round(imy / 2)];
for x=1:imx
    for y=1:imx
        im(x,y,:,:) = put_tensor(tensor_base,center,[x,y],1);
    end
end
im(10,10,:,:) = [1, 0; 0, 1];
fa = sum(sum(im,4),3); 

% Compute
[meshstruct, D, M, RHSbc] = computeDiffusionMatrix2D(im, [imx, imy]);
[potentials, mask_index] = computePotentials2D(D, meshstruct, M, RHSbc);
interactiveConductance2D(fa, [imx, imy], potentials, mask_index);

%%
% Input image circle
im = zeros(imx, imy, 2, 2);
tensor_base = [1, 0; 0, 0.2];
center = [round(imx / 2), round(imy / 2)];
for x=1:imx
    for y=1:imx
        im(x,y,:,:) = put_tensor(tensor_base,center,[x,y],2);
    end
end
im(10,10,:,:) = [1, 0; 0, 1];
im = im + eps;
fa = sum(sum(im,4),3); 

% Compute
[meshstruct, D, M, RHSbc] = computeDiffusionMatrix2D(im, [imx, imy]);
[potentials, mask_index] = computePotentials2D(D, meshstruct, M, RHSbc);
interactiveConductance2D(fa, [imx, imy], potentials, mask_index);
