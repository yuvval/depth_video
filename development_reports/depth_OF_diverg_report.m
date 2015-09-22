close all
clear

proj_root_path = '../';
addpath (proj_root_path);
initdirs;
OFcmap = 'jet';


%% Selecting immobile patches of an image
load ([proj_root_path '/preprocessed_videos/approaching_toward_fsmp_15_ppvid.mat'])

t=6;

imrgb = ppvid.frames{t};
imshow(imrgb);
im = double(rgb2gray(ppvid.frames{t}));
Nimg = numel(im);

%% Here is the image of optical flow (OF) magnitudes on a log10 scale
uOF_img = ppvid.uvOFs{t}(:,:,1);
vOF_img = ppvid.uvOFs{t}(:,:,2);
OFmag = sqrt(uOF_img.^2 + vOF_img.^2);
figure
imshow(log10(OFmag));
caxis('auto');
colormap(OFcmap);
h = colorbar;
ylabel(h, 'log10(Optical flow magnitude)');
title('optical flow (OF) magnitudes on a log10 scale')
shg
%% Here is the image of optical flow divergence 
%%
% 
%  It looks like the movevemnt of the person towards the camera induces a lot 
%  of variance around the movement area (divergence in both polarities).
% 

OFdiv = divergence(uOF_img, vOF_img);
figure;
imshow(OFdiv);
caxis([-1 1]);
%  caxis('auto');
colormap(OFcmap);
h = colorbar;
ylabel(h, 'Optical flow divergence');
title('optical flow divergence')

%% We extract candidates of image patches using MCG 
%%
% 
%  (http://arxiv.org/abs/1503.00848)

%%  Let's look on the OF and its statistics of the some candidates 

candidates = im2mcg(imrgb, 'fast');
% evaluate each of the masks
N_cand = size(candidates.scores,1);
mask_img = cell(N_cand,1);
k = 1;
for id = 1:N_cand
    mask_img{id} = get_mcg_mask_from_candidates(candidates, id);
    k=k+1;
end


for id = 1:5:30
    [bbox, mask_pxl_ids, mask_pxl_coords] = mask2bbox(mask_img{id});
    mask_crop = (mask_img{id}(bbox(1):bbox(2), bbox(3):bbox(4)));
    OFdiv_mask = (OFdiv(bbox(1):bbox(2), bbox(3):bbox(4))).*mask_crop;
    figure
    imshow(OFdiv_mask); caxis([-1 1]); 
    mean_OFdiv =  mean(OFdiv_mask(mask_crop));
    std_OFdiv =  std(OFdiv_mask(mask_crop));
    title(sprintf('mean = %2.3f, std = %2.3f', mean_OFdiv ,std_OFdiv)) ;
    caxis(1*[-1 1]);
    colormap(OFcmap);
end
%%
% 
%  We see that the OF divergence mean can be anything, and 
%  what discriminates motion is a large variance with the OF divergence
% 

%% Let's look at the statitics of the OF div. of all of the patches 
%%
% 
%  Each point on the scatter plots signify a patch candidate. 
%  The color coding for the scatter plot means the relative area of the
%  patch out of the image.
% 

metric_OFdiv = [];
for c = 1:N_cand
    [bbox, mask_pxl_ids, mask_pxl_coords] = mask2bbox(mask_img{c});
    mask_crop = (mask_img{c}(bbox(1):bbox(2), bbox(3):bbox(4)));
    OFdiv_mask = (OFdiv(bbox(1):bbox(2), bbox(3):bbox(4))).*mask_crop;
    mean_OFdiv(c) =  mean(OFdiv_mask(mask_crop));
    std_OFdiv(c) =  std(OFdiv_mask(mask_crop));
    metric_OFdiv(c) =  mean_OFdiv(c).^2 + (std_OFdiv(c)).^2;
    num_pxl(c) = sum(mask_crop(:));
end
% figure; hist(mean_OFdiv, -0.5:0.05:1.5); title('histogram of mean of OF divergence');
% figure; hist(std_OFdiv, 0:0.025:1.5); title('histogram of std of OF divergence');
figure; scatter(std_OFdiv.', mean_OFdiv.', [],log10(num_pxl/Nimg)); xlabel('std'); ylabel('mean'); 
title('scatter of mean vs std of OF divergence');
colormap('jet');
h = colorbar;
ylabel(h, 'log10(% area of image)');
figure; scatter(std_OFdiv.', mean_OFdiv.', [],log10(num_pxl/Nimg)); xlabel('std'); ylabel('mean'); 
title('(ZOOM IN) scatter of mean vs std of OF divergence');
colormap('jet');
h = colorbar;
ylabel(h, 'log10(% area of image)');
xlim([0, 0.1])
ylim([-.05 .05])

%% Selection of immobile patches
%%
% 
%  For adding a constraint of immobility to patches of an image, we would like to constrain areas
%  of the image that have both low abs mean and variance of OF divergence. 
%  Therefore, we can choose a metric that selects the patches of image with 
%  the lowest 2nd moment (mean^2 + var).
%  Since most of the patches with the lowest values are tiny, we only
%  select from patches that cover 1%-10% of the image, and we allow up to 20%
%  overlap with previously covered pixels. We select 20 such minimal OF
%  diversity patches.
% 

[~, sorted_minOF_ids] = sort(metric_OFdiv,'ascend');
patch_size_range = [0.01 0.1] ;% in percent of image size
valid_patch_ids = find(patch_size_range(1) <= num_pxl/Nimg & num_pxl/Nimg<=patch_size_range(2));
sorted_minOF_ids(~ismember(sorted_minOF_ids, valid_patch_ids)) = []; % trim too small/too big patches

covered_pixels_by_masks = zeros(size(im));
nbest = 20;
cnt = 1;
figure
for c=sorted_minOF_ids
    if cnt > nbest
        break;
    end
    [bbox, mask_pxl_ids, mask_pxl_coords] = mask2bbox(mask_img{c});
    
    % if this mask overlap by more than 20% with previous covered areas then we skip this mask
    if sum(covered_pixels_by_masks(mask_pxl_ids))/numel(mask_pxl_ids) > 0.2
        continue 
    else
        cnt = cnt+1;
    end
    covered_pixels_by_masks(mask_pxl_ids) = covered_pixels_by_masks(mask_pxl_ids) + 1; % mark current mask pixels as 'covered'
    
    mask_crop = (mask_img{c}(bbox(1):bbox(2), bbox(3):bbox(4)));
    OFdiv_mask = (OFdiv(bbox(1):bbox(2), bbox(3):bbox(4))).*mask_crop;
    
    subplot(5, ceil(nbest/5), cnt-1);
    imshow(OFdiv_mask); caxis([-1 1]);
    colormap(OFcmap);
    if cnt == ceil(nbest/5)
        title(sprintf('Optical flow divergence of each of the selected patches'));
    end
end
figure
covered_pixels_by_masks = covered_pixels_by_masks > 0;
covered_pixels_by_masks = uint8(covered_pixels_by_masks);
masked_im = imrgb.*repmat(covered_pixels_by_masks,[1,1,3]);
% masked_im = masked_im/255;
% masked_im= masked_im.^0.5; %gamma correction
imshow(masked_im);
title('masks (in white) of selected minimal OF div. patches');
% % caxis('auto');
return


