function [ColocalizationRatio] = FT_colocalization(mitotiff, sigtiff, nuctiff, mitofudgeFactor, cellfudgeFactor, nucfudgeFactor)
%% This function is used for colocalization calculation of two fluorescent channels from the same image.
% The mitotiff is used to determine location of mitochondria by using a
% mitochondrial fluorescent stain. The sigtiff is the TIFF image
% corresponding to the signal of interest (in our example, cytoplasmic YFP). The nuctiff is a TIFF image of a
% DAPI stain, which is used to subtract the nucleus from the analysis.
% These individual TIFF files were extracted from a multi-channel TIFF
% image using ImageJ.

%% To run this code, copy and paste the following line into the command line:

% [ColocalizationRatio] = FT_colocalization('dsRed.tif', 'YFP.tif', 'nucleus.tif', 0.8, 0.4, 0.95)
% INPUTS:
%   mitotiff = tiff stack of mitochondrial-localized signal 
%   sigtiff = tiff stack of cytoplasmically localized signal
%   nuctiff = tiff stack of DAPI signal
%   mitofudgeFactor, cellfudgeFactor, nucfudgeFactor = inputs to change the
%   sensitivity of the edge detection algorithm. Reasonable inputs range
%   from 0.3-1.5. Higher values set the segmentation threshold higher.

% OUTPUT:
%   ColocalizationRatio = n x 2 matrix where n is the number of slices in
%   your image stack (i.e. timepoints). The first column gives the stack
%   number and the second gives the colocalization ratio given as
%   [intensity of signal colocalized with mitochondria / intensity of
%   signal in the cell]
%% Import image stacks as matrices
MITOS=tiff2mat(mitotiff);
SIGNAL=tiff2mat(sigtiff);
NUCLEUS=tiff2mat(nuctiff);

% Convert images stacks to double type
mitoimage= double(MITOS);
cellimage= double(SIGNAL);

% Extract dimensions
N=size(cellimage,3);
M=size(cellimage,1);
L=size(cellimage,2);

%% Initialize matrices for faster calculations
total_cell_intensity=zeros(N,1);
average_intensity_mitos=zeros(N,1);
mitofillpixels=zeros(M,L,N);
mitopixels=zeros(M,L,N);
cell_no_mitos=zeros(M,L,N);
BWs=zeros(M,L,N);
BWsdil=zeros(M,L,N);
BWdfill=zeros(M,L,N);
total_intensity_mitos_new=zeros(1,N);

%% Dynamically define mito outline for each frame using Sobel edge detection and object filling
for i=1:N
image= mitoimage(:,:,i);
[~, threshold] = edge(image, 'sobel');
%mitofudgeFactor = .9;

% Mito Outline
BWs = edge(image,'sobel', threshold * mitofudgeFactor); % fudge factor to fill gaps 
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);

% Mito Partial Fill
BWsdil(:,:,i) = imdilate(BWs, [se90 se0]);

% Create logical matrix for mitochondrial pixels
BWdfill(:,:,i) = imfill(BWsdil(:,:,i), 'holes');

%Create logical indexing array of mito pixels
mitofillpixels(:,:,i) = BWdfill(:,:,i);
mitopixels(:,:,i) = BWsdil(:,:,i);
end

 figure, imshow(BWs(:,:,1)), title('Mito Outline');
 figure, imshow(BWdfill(:,:,1)); title('Mito Fill');
%% Create logical array of cell pixels
%Use first slice to define cell outline - assumes that the signal in
%cellimage fills the cytoplasm in this frame
image= cellimage(:,:,1);

% Apply Sobel edge detection
[~, threshold] = edge(image, 'sobel');
%cellfudgeFactor = .55;
BWs = edge(image,'sobel', threshold * cellfudgeFactor);
%figure, imshow(BWs), title('Cell Outline');

se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);

BWsdil = imdilate(BWs, [se90 se0]);

% Create logical matrix BWdfill with 1 if pixel lies in cell and 0 if pixel does not lie in
% the cell
BWdfill = imfill(BWsdil, 'holes');
figure, imshow(BWdfill);
title('Cell Fill');

%Create logical indexing array of cell pixels
cellpixels = BWdfill;


%% Define logical array for nuclear pixels
%Use frame 90 to define nucleus outline
image= NUCLEUS(:,:,90);

% Apply Sobel edge detection
[~, threshold] = edge(image, 'sobel');
%nucfudgeFactor = .85;
BWs = edge(image,'sobel', threshold * nucfudgeFactor);
%figure, imshow(BWs), title('Nucleus Outline');
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWsdil = imdilate(BWs, [se90 se0]);
%figure, imshow(BWsdil), title('Nucleus Partial Fill');

% Logical nuclear pixel matrix
BWdfill = imfill(BWsdil, 'holes');
figure, imshow(BWdfill);
title('Nucleus Fill');

%Create logical indexing array of nucleus pixels
nucpixels = BWdfill;

%Subtract pixels in mitochondria and nucleus from cytoplasmic cell mask
for i=1:N
cell_no_mitos(:,:,i)=cellpixels-mitofillpixels(:,:,i)-nucpixels;
end 
figure; imshow(cell_no_mitos(:,:,1))
title('Cell Area (not mitos or nucleus)')

figure;imshowpair(mitofillpixels(:,:,1),MITOS(:,:,1))

%Calculate intensity on cell without mitos or nucleus
for i=1:size(cellimage,3)
    image=SIGNAL(:,:,i,1);
   A(:,:,i)=logical(cell_no_mitos(:,:,i));
total_cell_intensity(i)=sum(sum(image(A(:,:,i))));
end

% Calculate average intensity of pixels within mitochondria
for i=1:size(cellimage,3)
    image=cellimage(:,:,i);
    MITOS=logical(mitopixels(:,:,i));
    Mitos=image(MITOS);
    M=sum(sum(mitopixels(:,:,i)));
total_intensity_mitos_new(i)=sum(sum(Mitos));
average_intensity_mitos(i)=sum(sum(Mitos))/M;
end

% Sum over the cell area and mito area for each frame
for i=1:N
 a(i)=sum(sum(cell_no_mitos(:,:,i)));
 b(i)=sum(sum(mitopixels(:,:,i)));
end

inv_tot_int_mitos= total_intensity_mitos_new';

% Calculate ratio of signal colocalized with mitochondria vs. cytoplasm
for i=1:N
mitotoCell(i)=(inv_tot_int_mitos(i)/b(i))./(total_cell_intensity(i)/a(i));
end
%% Generate ratio of GFP intensity colocalized with mitochondria to GFP intensity colocalized with cytoplasm at each time point
x=[1:size(mitoimage,3)];
stackInterval=1; % this can be changed to reflect different times between images in the stack
t=stackInterval.*x;
y=mitotoCell;

ColocalizationRatio = horzcat(t',y');
end

