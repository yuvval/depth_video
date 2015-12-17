function vis_dep_est(ppvid, video_name, depth_estimate, mnFrameID, ...
                     mxFrameID, fname_gif, visibility_flag, ptmp)


    if nargin<6
        fname_gif = []; % indicates that we don't save this visualization
    end
    
    if nargin < 7
        visibility_flag = false;
    end
    
    % an API for temporary experimental stuff
    if nargin < 8
        ptmp = [];
    else
        boxes = ptmp;
    end
    
    
    if visibility_flag 
        visibility_flag = 'on';
    else
        visibility_flag = 'off';
    end

    N = size(ppvid.depth_frames,3);
    if mxFrameID == -1
        mxFrameID = N;
    end


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

    davg_qp = mean(depth_estimate(:));
    dstd_qp = std(depth_estimate(:));
    dmn_qp = min(depth_estimate(:));
    dmx_qp = max(depth_estimate(:));

    tmp = ppvid.rgb_frames(:,:,:,1);
    n255 = 1; 
    if max(tmp(:)) > 1
        n255 = 255;
    end

    %%
    %imsize = size(tmp(:,:,1));

    %%
    h = figure('Visible',visibility_flag);
    a = axes('Visible',visibility_flag);

    t=1;
    for n=mnFrameID:mxFrameID
        subplot(2,9,1:4)
        rgb_frame = ppvid.rgb_frames(:,:,:,n)/n255;
        if is_in_str('eigen', fname_gif)
            rgb_frame = imresize(rgb_frame, [109, 147]);
        end % if is_in_str('eigen'..

        imshow(rgb_frame);
        title(video_name,'interpreter','none');

        if exist('boxes')
            colors = {'r', 'y'};
            for trkr = 1:length(boxes)
                
                
                d = t;
                x1 = boxes{trkr}(d,1);
                x2 = boxes{trkr}(d,2);
                y1 = boxes{trkr}(d,3);
                y2 = boxes{trkr}(d,4);
                %             label = sprintf('%s, %2.3f', label, ppvid.scores{t}(d));
                %             feat_name = 'velocity_abs';
                %             feat_id = find(ismember(tracker_feats.names, feat_name));
                %             feat_val = ppvid.scores{t}(d);
                %             feat_val = tracker_feats.values{t}(d_prev, d, feat_id);
                line([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]',   'color', colors{trkr},'linewidth', 1.5, 'linestyle', '-');
                
            end
            
            
        end
%         disp('press any key');pause
        drawnow;
        
        
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
        
        subplot(2,9,10:13)
        
        tmp = ppvid.depth_frames(:,:,n);
        %     tmp(tmp > dmx) = 0;
        imshow(tmp);
        caxis([dmn dmx])
        colormap((jet));
        title('raw depth estimation');

        % QP estimated depth
        subplot(2,9,14:17)
        
        tmp = depth_estimate(:,:,t);
        imshow(tmp);
        caxis([dmn_qp dmx_qp])
        colormap((jet));
        title('QP estimated depth');
        subplot(2,9,18)
        imshow(zeros(20,1));
        caxis([dmn dmx])
        colormap((jet));
        colorbar;
        
%         disp('press key');pause
        save_animated_gif_frame(fname_gif, t==1, h);

        t=t+1;
    end
    if ~isempty(fname_gif)
        save_animated_gif_frame(fname_gif, t==1, h)
        save_animated_gif_frame(fname_gif, t==1, h)
    end
end
