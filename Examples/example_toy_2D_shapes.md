## Example 1: Compute conductance of different shapes, and plot it interactively

Code in [Examples/Connectivity/example_2D_tensor_test_cases.m](Connectivity/example_2D_tensor_test_cases.m)

In this example, we generate data ourselves with 2D tensors making different shapes, and we can check how our condutance work interactively: click on a voxel and an image of the conductance map from that point to the rest will be generated.

For isotropic tensors:
```
im = zeros(imx, imy, 2, 2);
im(:,:,1,1) = 0.2;
im(:,:,2,2) = 0.2;
fa = sum(sum(im,4),3);

[meshstruct, D, M, RHSbc] = computeDiffusionMatrix2D(im, [imx, imy]);
[potentials, mask_index] = computePotentials2D(D, meshstruct, M, RHSbc);
interactiveConductance2D(fa, [imx, imy], potentials, mask_index);
```
<p align="center">
<img src="Images/isotropic_tensors.png" width="150" hspace="40"> <img src="Images/iso.png" width="200"> <img src="Images/iso2.png" width="200">
</p>

For a cross:
```
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

[meshstruct, D, M, RHSbc] = computeDiffusionMatrix2D(im, [imx, imy]);
[potentials, mask_index] = computePotentials2D(D, meshstruct, M, RHSbc);
interactiveConductance2D(fa, [imx, imy], potentials, mask_index);
```
<p align="center">
<img src="Images/cross_tensors.png" width="150" hspace="40"> <img src="Images/cross.png" width="200"> <img src="Images/cross2.png" width="200">
</p>

For a star:
```
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

[meshstruct, D, M, RHSbc] = computeDiffusionMatrix2D(im, [imx, imy]);
[potentials, mask_index] = computePotentials2D(D, meshstruct, M, RHSbc);
interactiveConductance2D(fa, [imx, imy], potentials, mask_index);
```
<p align="center">
<img src="Images/star_tensors.png" width="150" hspace="40"> <img src="Images/star.png" width="200"> <img src="Images/star2.png" width="200">
</p>

For a circle:
```
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

[meshstruct, D, M, RHSbc] = computeDiffusionMatrix2D(im, [imx, imy]);
[potentials, mask_index] = computePotentials2D(D, meshstruct, M, RHSbc);
interactiveConductance2D(fa, [imx, imy], potentials, mask_index);
```
<p align="center">
<img src="Images/circular_tensors.png" width="150" hspace="40"> <img src="Images/circle2.png" width="200"> <img src="Images/circle3.png" width="200">
</p>
