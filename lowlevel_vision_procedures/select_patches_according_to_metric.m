function [selected_masks_ids, covered_pixels_by_masks] = ...
    select_patches_according_to_metric(metrics, ...
                                     selection_type, selection_type_params)

switch selection_type
    case 'immobile'
        [selected_masks_ids, covered_pixels_by_masks] = ...
            select_immobile_patches(metrics, selection_type_params);
    % selection_type_params should consist of:
    % params.patch_sizes_range  = [min_size_in_percent, max_size_in_percent] 
    % params.num_pxl_per_mask   = an array that contains the number of pixels per
    %                           mask
    % params.imsize             = [n_rows, n_cols] image size 
    % params.nbest              = number of best masks to select
    % params.masks              = objects proposal masks

    otherwise
        error ('unknown selection_type (%s)', selection_type);
        
end

function [selected_masks_ids, covered_pixels_by_masks] = ...
            select_immobile_patches(metrics, params)

[~, sorted_ids] = sort(metrics,'ascend');
patch_sizes_range = params.patch_sizes_range;% in percent of image size
num_pxl = params.num_pxl_per_mask;
Nimg = prod(params.imsize);

valid_patch_ids = find(patch_sizes_range(1) <= num_pxl/Nimg & num_pxl/Nimg<=patch_sizes_range(2));
sorted_ids(~ismember(sorted_ids, valid_patch_ids)) = []; % trim too small/too big patches

covered_pixels_by_masks = zeros(params.imsize);
nbest = params.nbest;
selected_masks_ids = nan(1, nbest);

cnt = 1;
for c=sorted_ids
    if cnt > nbest
        break;
    end
    [~, mask_pxl_ids, ~] = mask2bbox(params.masks{c});
    
    % if this mask overlap by more than 20% with previous covered areas then we skip this mask
    if sum(covered_pixels_by_masks(mask_pxl_ids))/numel(mask_pxl_ids) > 0.2
        continue 
    else
        selected_masks_ids(cnt) = c;
        cnt = cnt+1;
        covered_pixels_by_masks(mask_pxl_ids) = covered_pixels_by_masks(mask_pxl_ids) + 1; % mark current mask pixels as 'covered'
    end
end
covered_pixels_by_masks = uint8(covered_pixels_by_masks > 0);

