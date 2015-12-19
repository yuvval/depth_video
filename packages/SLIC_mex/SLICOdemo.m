%======================================================================
%SLICO demo
% Copyright (C) 2015 Ecole Polytechnique Federale de Lausanne
% File created by Radhakrishna Achanta
% Please also read the copyright notice in the file slicmex.c 
%======================================================================
%Input:
%[1] 8 bit images (color or grayscale)
%[2] Number of required superpixels (optional, default is 200)
%
%Ouputs are:
%[1] labels (in raster scan order)
%[2] number of labels in the image (same as the number of returned
%superpixels
%
%NOTES:
%[1] number of returned superpixels may be different from the input
%number of superpixels.
%[2] you must compile the C file using mex slicmex.c before using the code
%below
%----------------------------------------------------------------------
% How is SLICO different from SLIC?
%----------------------------------------------------------------------
% 1. SLICO does not need compactness factor as input. It is calculated
% automatically
% 2. The automatic value adapts to the content of the superpixel. So,
% SLICO is better suited for texture and non-texture regions
% 3. The advantages 1 and 2 come at the cost of slightly poor boundary
% adherences to regions.
% 4. This is also a very small computational overhead (but speed remains
% almost as fast as SLIC.
% 5. There is a small memory overhead too w.r.t. SLIC.
% 6. Overall, the advantages are likely to outweigh the small disadvantages
% for most applications of superpixels.
%======================================================================
%img = imread('someimage.jpg');

rescale_mul = 0.1513;
% rescale_mul = 0.3333;

img = imread('bee.jpg');
img = im;
% img = imresize(img, rescale_mul);
if max(img(:)) <= 1
    img = int8(img*255);
end
[labels, numlabels] = slicomex(img,2000);%numlabels is the same as number of superpixels
figure;
% imagesc(labels);

labels = double(labels);
labels = imresize(labels, rescale_mul, 'nearest');
% labels = round(labels);
% labels = medfilt2(labels);
% labels = medfilt2(labels);
% labels(labels > numlabels) = numlabels;
% labels(labels < 0) = 0;
% imagesc(labels);
[gx,gy] = gradient(labels);
e =1- ((gx.^2+gy.^2)==0);
E = repmat(e, [1,1,3]);
% imagesc(E);
img = imresize(img, rescale_mul);
imagesc(double(img)/255 + E);shg


