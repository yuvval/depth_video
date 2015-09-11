function [bbox, mask_pxl_ids, mask_pxl_coords] = mask2bbox(mask_img)
% finds bounding box coordinates of an image mask
% bbox = [x1 x2 y1 y2];
% returns also byproducts that might be useful:
% mask_pxl_ids = a vector of indices of the mask pixels
% mask_pxl_coords = an x,y coordinates of the mask pixels

mask_pxl_ids = find(mask_img(:));
[mask_pxl_coords(:,1), mask_pxl_coords(:,2)] = ind2sub(size(mask_img), mask_pxl_ids);


x1 = min(mask_pxl_coords(:,1));
x2 = max(mask_pxl_coords(:,1));
y1 = min(mask_pxl_coords(:,2));
y2 = max(mask_pxl_coords(:,2));
bbox = [x1 x2 y1 y2];
