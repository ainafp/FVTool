## Example 2: Compute conductance of Fibercup data, and plot it interactively

Code in [Examples/Connectivity/example_3D_tensor_fibercup.m](Examples/Connectivity/example_3D_tensor_fibercup.m)

Data can be downloaded in [www.tractometer.org](http://www.tractometer.org/original_fibercup/data/).
In this example, we use the file [acq-averaged_b-1500.nii.gz](http://www.tractometer.org/downloads/downloads/fibercup/dwi/acq-averaged_b-1500.nii.gz)

First we need to run DSI Studio to get the tensors, given some inputs: DWI data (data.nii.gz, bvals, bvecs), and a brain mask (nodif_brain_mask.nii.gz)
dsi_studio_run --action=src --source=data.nii.gz --bval=bvals --bvec=bvecs --output=data.src.gz
dsi_studio_run --action=rec --source=data.src.gz --mask=nodif_brain_mask.nii.gz --method=1 --output_tensor=1 --output_dif=1 --output=data.src.gz.012fy.dti.fib.gz
The output of these instructions will be the file data.src.gz.012fy.dti.fib.gz, and contains the tensors resulting from the DWI data:
Alternatively, you could generate your own tensors and load them as a volume of dimensions (im_x, im_y, im_z, tensor_x, tensor_y).
The results of DSI Studio
<p align="center">
<img src="Examples/Images/fibercup.png" width="550">
</p>



