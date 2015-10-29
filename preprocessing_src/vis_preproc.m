function vis_preproc(ppvid, video_name, fname_gif)
if nargin<3
    fname_gif = []; % indicates that we don't save this visualization
end

N = size(ppvid.depth_frames,3);


%%
gt_depth_frames1 = ppvid.gt_depth_frames(:,:,1);
gt_depth_frames1(gt_depth_frames1 == 0 ) = []; % trim non estimated values.
davg_gt = mean(gt_depth_frames1(:));
dstd_gt = std(gt_depth_frames1(:));
% dmn_gt = max(davg_gt-2*dstd_gt, 0);
% dmx_gt = davg_gt+2*dstd_gt;

dmn_gt = min(gt_depth_frames1(:));
dmx_gt = max(gt_depth_frames1(:));

rawest_depth_frames1 = ppvid.depth_frames(:,:,1);
davg = mean(rawest_depth_frames1(:));
dstd = std(rawest_depth_frames1(:));
% dmn = max(davg-2*dstd, 0);
% dmx = davg+2*dstd;

dmn = min(rawest_depth_frames1(:));
dmx = max(rawest_depth_frames1(:));

tmp = ppvid.rgb_frames(:,:,:,1);
n255 = 1; 
if max(tmp(:)) > 1
    n255 = 255;
end

%%
imsize = size(tmp(:,:,1));
[OFdiv_frames] = deal(zeros([imsize N]));
for n=1:N
    if n<N
        uOF = ppvid.opflow_frames{n}(:,:,1);
        vOF = ppvid.opflow_frames{n}(:,:,2);
        OFdiv_frames(:,:,n) = divergence(uOF, vOF);    
    end
end
mOFdiv = mean(OFdiv_frames(:));
sOFdiv = std(OFdiv_frames(:));
sigmoid = @(x, a, b) (1./(1+exp(-a*(x-b))));

%%
for n=1:N
    subplot(2,9,1:4)
    rgb_frame = ppvid.rgb_frames(:,:,:,n)/n255;
    imshow(rgb_frame);
    title(video_name,'interpreter','none');

    if ~isempty(ppvid.gt_depth_frames)
        subplot(2,9,5:8)
        imshow(ppvid.gt_depth_frames(:,:,n));
        caxis([dmn_gt dmx_gt])
        colormap((jet));
        title('ground truth depth');
        subplot(2,9,9)
        imshow(zeros(20,1));
        caxis([dmn_gt dmx_gt])
        colormap((jet));
        colorbar;
    end
    
    subplot(2,9,14:17)
    
    tmp = ppvid.depth_frames(:,:,n);
%     tmp(tmp > dmx) = 0;
    imshow(tmp);
    caxis([dmn dmx])
    colormap((jet));
    title('raw depth estimation');
    subplot(2,9,18)
    imshow(zeros(20,1));
    caxis([dmn dmx])
    colormap((jet));
    colorbar;

    % OF divergence
    if n<N

    subplot(2,9,10:13)

    % color render the OF divergence
    gim = double(rgb2gray(rgb_frame));
    OFdivRender = repmat(gim, 1,1,3);
    normOFdiv = (OFdiv_frames(:,:,n) - mOFdiv)/sOFdiv;
    sig_a = 10;
    normOFdiv_pos = normOFdiv.*(normOFdiv > 0);
    normOFdiv_neg = normOFdiv.*(normOFdiv <= 0);
    color_gamma_pos = 2*(1-sigmoid(normOFdiv_pos,sig_a,0));
    color_gamma_neg = 2*(sigmoid(normOFdiv_neg,sig_a,0));
    ch = 1;
    OFdivRender(:,:,ch) = OFdivRender(:,:,ch).^color_gamma_pos;
    ch = 3;
    OFdivRender(:,:,ch) = OFdivRender(:,:,ch).^color_gamma_neg;
    imshow(OFdivRender);
    title('Divergence of optical flow')
    
    if ~isempty(fname_gif)
        save_animated_gif_frame(fname_gif, n==1);        
    end
    
    end
    pause(0.001)
end
if ~isempty(fname_gif)
    save_animated_gif_frame(fname_gif, false);
    save_animated_gif_frame(fname_gif, false);
end
