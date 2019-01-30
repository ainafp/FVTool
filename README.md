# FVT4DWI: FVTool for Diffusion-weighted Imaging

## Quantification of Structural Brain Connectivity via a Conductance Model

This is a toolbox to study structural brain connectivity using a combination 
of differential Maxwell’s equations and Kirchhoff’s circuit laws, resulting in 
an equation similar to the heat equation. By solving this partial differential equation
for a certain current configuration between 2 voxels, we find the potential map
for that specific configuration. We further compute the electric conductance 
between each pair of voxels from potential maps, to which all diffusion paths between the pair contribute.  
The same measure can also be computed between a pair of regions of interest (ROIs) 
instead of voxels, by distributing the currents among the ROI voxels.

<p align="center">
<img src="conductance.png" width="550">
</p>

### How do I use it?

You can find a generic example for structural brain connectivity in [Examples/Connectivity/run_conductance_model.m](Examples/Connectivity/run_conductance_model.m)

### Reference

This work is in press in NeuroImage:

[Frau-Pascual, Aina, Morgan Fogarty, Bruce Fischl, Anastasia Yendiki, and Iman Aganj. "Quantification of structural brain connectivity via a conductance model." NeuroImage (2019).](https://www.sciencedirect.com/science/article/pii/S1053811919300333)

### Inspiration

This toolbox is a fork of FVTool, and was extended for dealing with tensors and 
modified for its use in structural connectivity settings.
