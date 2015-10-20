function vis_preproc(ppvid)
close all
N = size(ppvid.depth_frames,3);


gt_depth_frames1 = ppvid.gt_depth_frames(:,:,1);
davg_gt = mean(gt_depth_frames1(:));
dstd_gt = std(gt_depth_frames1(:));
dmn_gt = max(davg_gt-dstd_gt, 0);
dmx_gt = davg_gt+dstd_gt;

rawest_depth_frames1 = ppvid.depth_frames(:,:,1);
davg = mean(rawest_depth_frames1(:));
dstd = std(rawest_depth_frames1(:));
dmn = max(davg-dstd, 0);
dmx = davg+dstd;

tmp = ppvid.rgb_frames(:,:,:,1);
n255 = 1; 
if max(tmp(:)) > 1
    n255 = 255;
end
for n=1:N
    subplot(2,2,1)
    imshow(ppvid.rgb_frames(:,:,:,n)/n255);

    if ~isempty(ppvid.gt_depth_frames)
        subplot(2,2,2)
        imshow(ppvid.gt_depth_frames(:,:,n));
        caxis([dmn_gt dmx_gt])
        title('ground truth depth');
    end
    
    subplot(2,2,4)
    
    tmp = ppvid.depth_frames(:,:,n);
%     tmp(tmp > dmx) = 0;
    imshow(tmp);
%     caxis([dmn dmx])
    if n == 1
        caxis('auto')
        ax = gca;
        my_clim = ax.CLim;
    end
    caxis(my_clim);
    title('raw depth estimation');

    pause(0.1)
end