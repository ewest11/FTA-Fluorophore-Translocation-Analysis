# FTA-Fluorophore-Translocation-Analysis

This code was developed to analyze colocalization data from Gutnick et al., 2018. The code was developed on 2016b and requires the Matlab [Image Processing Toolbox](https://www.mathworks.com/products/image.html) and the OME's [Bio-Format's library](http://www.openmicroscopy.org/bio-formats/) (developed with version 5.9.1).

The code is developed to segment 2-D fluorescent images with cytoplasmic, nuclear, and organel (mitochondrial) signals and calculate fluorescent colocalization between cytoplasmic and organelle channels. Specifically, the code calculates the ratio of signal colocalization between two channels across an image stack.

Calculating Colocalization Between Fluorescent Channels
------------------------
In this folder you will find all functions and sample images needed to run the main function, FT_colocalization. The function inputs and outputs are as follows:

```matlab
  >> [ColocalizationRatio] = FT_colocalization(mitotiff, sigtiff, nuctiff, mitofudgeFactor, cellfudgeFactor, nucfudgeFactor)

 INPUTS:
%   mitotiff = tiff stack of mitochondrial-localized signal 
%   sigtiff = tiff stack of cytoplasmically localized signal
%   nuctiff = tiff stack of DAPI signal
%   mitofudgeFactor, cellfudgeFactor, nucfudgeFactor = inputs to change the
%   sensitivity of the edge detection algorithm. Reasonable inputs range
%   from 0.3-1.5. Higher values set the segmentation threshold higher.

 OUTPUT:
%   ColocalizationRatio = n x 2 matrix where n is the number of slices in
%   your image stack (i.e. timepoints). The first column gives the stack
%   number and the second gives the colocalization ratio given as
%   (intensity of signal colocalized with mitochondria / total intensity of
%   signal in the cell)
```

To test this function on our sample images, run the following line:
```matlab
>> [ColocalizationRatio] = FT_colocalization('mito.tif', 'YFP.tif', 'nucleus.tif', 0.8, 0.4, 0.95)
```

Questions?
------------
Please reach out to [Emma](mailto:emma_west@g.harvard.edu) with any questions about installing and running the scripts. 
