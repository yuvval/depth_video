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
clear
close all
%img = imread('someimage.jpg');

% rescale_mul = 0.1513;
rescale_mul = 0.3333;

img = imread('bee.jpg');
% img = im;
img = imresize(img, rescale_mul);
if max(img(:)) <= 1
    img = int8(img*255);
end
[labels, numlabels] = slicomex(img,2000);%numlabels is the same as number of superpixels
figure;
% imagesc(labels);

labels = double(labels);
% labels = imresize(labels, rescale_mul, 'nearest');
% img_edges = double(edge(rgb2gray(img), 'sobel'));
img_edges = double(imgradient(rgb2gray(img), 'sobel'));
img_edges = img_edges/max(img_edges(:));
% imagesc(repmat(img_edges,[1 1 3]) + double(img)/255)
imagesc(img_edges )
title('edge detector ');
hfig = gcf;
save_plot(['fig' num2str(hfig.Number)], gcf, '~/www/figs/', true, false);

[gx,gy] = gradient(labels);
e =1- ((gx.^2+gy.^2)==0);
E = repmat(e, [1,1,3]);

[sp_edges_y, sp_edges_x] = find(e);
sp_edges_yx = [sp_edges_y, sp_edges_x];


edges_grad_yx = nan*sp_edges_yx;
% edges_grad_angle = [];
for k=1:length(edges_grad_yx)
    edges_grad_yx(k,:) = [gy(sp_edges_yx(k,1), sp_edges_yx(k,2)) gx(sp_edges_yx(k,1), sp_edges_yx(k,2))];
    edges_grad_yx(k,:) = edges_grad_yx(k,:)/norm(edges_grad_yx(k,:));
    %     edges_grad_angle(k) = phase( edges_grad_yx(k,2) + 1i*edges_grad_yx(k,1)) ;
end

% neigh_sp_edges_yx_pos = nan*sp_edges_yx;
% neigh_sp_edges_yx_neg = nan*sp_edges_yx;
neigh_sp_pairs = nan*sp_edges_yx;
tmp = 0;

neigh_sp_edges_yx_pos = round(sp_edges_yx + edges_grad_yx);

neigh_sp_edges_yx_pos(neigh_sp_edges_yx_pos(:,1) < 1, 1) = 1;
neigh_sp_edges_yx_pos(neigh_sp_edges_yx_pos(:,2) < 1, 2) = 1;
neigh_sp_edges_yx_pos(neigh_sp_edges_yx_pos(:,1) > size(labels,1), 1) = size(labels,1);
neigh_sp_edges_yx_pos(neigh_sp_edges_yx_pos(:,2) > size(labels,2), 2) = size(labels,2);

neigh_sp_edges_yx_pos2 = round(sp_edges_yx + 2*edges_grad_yx);

neigh_sp_edges_yx_pos2(neigh_sp_edges_yx_pos2(:,1) < 1, 1) = 1;
neigh_sp_edges_yx_pos2(neigh_sp_edges_yx_pos2(:,2) < 1, 2) = 1;
neigh_sp_edges_yx_pos2(neigh_sp_edges_yx_pos2(:,1) > size(labels,1), 1) = size(labels,1);
neigh_sp_edges_yx_pos2(neigh_sp_edges_yx_pos2(:,2) > size(labels,2), 2) = size(labels,2);

neigh_sp_edges_yx_neg = round(sp_edges_yx - edges_grad_yx);

neigh_sp_edges_yx_neg(neigh_sp_edges_yx_neg(:,1) < 1, 1) = 1;
neigh_sp_edges_yx_neg(neigh_sp_edges_yx_neg(:,2) < 1, 2) = 1;
neigh_sp_edges_yx_neg(neigh_sp_edges_yx_neg(:,1) > size(labels,1), 1) = size(labels,1);
neigh_sp_edges_yx_neg(neigh_sp_edges_yx_neg(:,2) > size(labels,2), 2) = size(labels,2);

neigh_sp_edges_yx_neg2 = round(sp_edges_yx - 2*edges_grad_yx);

neigh_sp_edges_yx_neg2(neigh_sp_edges_yx_neg2(:,1) < 1, 1) = 1;
neigh_sp_edges_yx_neg2(neigh_sp_edges_yx_neg2(:,2) < 1, 2) = 1;
neigh_sp_edges_yx_neg2(neigh_sp_edges_yx_neg2(:,1) > size(labels,1), 1) = size(labels,1);
neigh_sp_edges_yx_neg2(neigh_sp_edges_yx_neg2(:,2) > size(labels,2), 2) = size(labels,2);

same_superpxl_ids_pos = neigh_sp_edges_yx_pos(:,1) == neigh_sp_edges_yx_pos(:,2);

sp_edges_yx = [ sp_edges_yx ];
neigh_sp_edges_yx = [neigh_sp_edges_yx_pos];
% neigh_sp_edges_yx = neigh_sp_edges_yx_pos;
% neigh_sp_edges_yx(same_superpxl_ids_pos, :) = neigh_sp_edges_yx_neg(same_superpxl_ids_pos, :);

for k=1:length(neigh_sp_edges_yx)
    E(neigh_sp_edges_yx(k,1), neigh_sp_edges_yx(k,2), :) = [0 0 1];
end

E = 0*e;

sp_edges_ids = sub2ind(size(labels), sp_edges_yx(:,1), sp_edges_yx(:,2));
neigh_sp_edges_pos_ids = sub2ind(size(labels), neigh_sp_edges_yx_pos(:,1), neigh_sp_edges_yx_pos(:,2));
neigh_sp_edges_pos2_ids = sub2ind(size(labels), neigh_sp_edges_yx_pos2(:,1), neigh_sp_edges_yx_pos2(:,2));
neigh_sp_edges_neg_ids = sub2ind(size(labels), neigh_sp_edges_yx_neg(:,1), neigh_sp_edges_yx_neg(:,2));
neigh_sp_edges_neg2_ids = sub2ind(size(labels), neigh_sp_edges_yx_neg2(:,1), neigh_sp_edges_yx_neg2(:,2));

neigh_sp_pairs = [labels(sp_edges_ids) labels(neigh_sp_edges_pos_ids)];
neigh_sp_pairs_edge_energies = nan(length(neigh_sp_pairs),1);

edge_weighs = img_edges(sp_edges_ids) + img_edges(neigh_sp_edges_pos_ids) + ...
    img_edges(neigh_sp_edges_pos2_ids) + img_edges(neigh_sp_edges_neg_ids) +...
    img_edges(neigh_sp_edges_neg2_ids);

E(sp_edges_ids) = edge_weighs;



[u_neigh, i_unique, i_orig ] = unique(neigh_sp_pairs, 'rows');
unique_summed_edge_weighs = 0*i_unique;
marginalized_edge_weighs = 0*edge_weighs;

cnt = 1;
for pair = i_orig.'
    unique_summed_edge_weighs(pair) = unique_summed_edge_weighs(pair) + edge_weighs(cnt);
    cnt = cnt+1;
end
figure;
imagesc(E);shg
title('edges marginalized weights');
hfig = gcf;
save_plot(['fig' num2str(hfig.Number)], gcf, '~/www/figs/', true, false);

E=0*E;

cnt=1;
for pair = i_orig.'
    marginalized_edge_weighs(cnt) = unique_summed_edge_weighs(pair);
    cnt = cnt+1;
end

E(sp_edges_ids) = marginalized_edge_weighs;


figure

imagesc(E);shg
title('edges marginalized weights');
% img = imresize(img, rescale_mul);
% imagesc(double(img)/255 + E);shg
% title('image edges overlapping superpixels edges');

hfig = gcf;
save_plot(['fig' num2str(hfig.Number)], gcf, '~/www/figs/', true, false);


