function [ppvid, res_fname, fname_OF, depth_im ] = preprocess_video_3d(vid_fname, detection_thresh, frame_sample_interval, take_top_n_detections )
%% init
if nargin<1
    vid_fname = '../videos/after2sec_diagonal.avi'; % Person approaches a chair.
end

if nargin < 2
    detection_thresh = -1.5;
end
if nargin < 3
    frame_sample_interval = 15; % Sample a frame from video once every X frames.
end
if nargin < 4
    take_top_n_detections = 3;
end

trim_first_seconds = 3;

addpath ../
addpath ../optical_flow
addpath ../optical_flow/algorithms/CLG-TV/

init_obj_detect;
fname_split = regexp(vid_fname, '[\./]', 'split');
vid_name = fname_split{end-1};

obj = VideoReader(vid_fname);
video = obj.read();

%% Iterate on frames.
[boxes, classes, scores, centers, projected_centers ] = deal({});

Nframes = size(video,4);
t = 1; % Sampled frames counter.
tic 
depth_im{1} = getDepth_Fayao(video(:,:,:,1)); % eval the depth for the 1st frame
toc
for k=1 + (trim_first_seconds*30):frame_sample_interval:Nframes
    
    % Evaluate detections.
    im=video(:,:,:,k);
    tic
    [boxes{t}, classes{t}, scores{t}, classes_names]  = obj_detect_frame(im, detection_thresh);
    toc
    boxes{t} = ceil(boxes{t}); % convert to integer
    % take only top-n detections
    if take_top_n_detections ~= inf
        % sort detection scores
        [~, sorted_ids] = sort(scores{t}, 'descend');
        n_det = min(numel(sorted_ids), take_top_n_detections);
        scores{t} = scores{t}(sorted_ids(1:n_det));
        boxes{t} = boxes{t}(sorted_ids(1:n_det), :);
        classes{t} = classes{t}(sorted_ids(1:n_det));
    end
    
    if k <= (Nframes-frame_sample_interval) % we can't eval OF for the last frame, so we skip it.
    %% Evaluate depth
        im2 = video (:,:,:,k+frame_sample_interval);
        tic
        depth_im{t+1} = getDepth_Fayao(im2);
        toc
    
    %% Evaluate optical flow. + depth flow
        
        [im1gray, im2gray] = acquistionSeq(im, im2);
        
        uvOF = OpticalFlowCLG_TV(im1gray, im2gray);
        fname_OF = ['../OF_frames/' vid_name '_OF_' num2str(t) '_' num2str(frame_sample_interval)];
        save (fname_OF, 'uvOF');
        % Evaluate detection boxes mean optical flow.
        
        n_detections = size(boxes{t},1);
        [centers{t}, projected_centers{t}] = deal(nan(n_detections ,3));
        
        max_x_dim = size(im,2);
        max_y_dim = size(im,1);
        for d=1:n_detections
            x1 = boxes{t}(d,1);
            x2 = boxes{t}(d,2);
            y1 = boxes{t}(d,3);
            y2 = boxes{t}(d,4);
            
            u = uvOF(:,:,1).';
            v = uvOF(:,:,2).';
            u_box = u(x1:x2, y1:y2);
            u_avg_box = mean(u_box(:));
            
            v_box = v(x1:x2, y1:y2);
            v_avg_box = mean(v_box(:));
            
            center_x = (x1+x2)/2;
            center_y = (y1+y2)/2;

            len_x = x2-x1;
            len_y = y2-y1;
            % taking the only the depth measurement around the center, so
            % we can avoid the background depth.
            trim_len_x = ceil(len_x/4);
            trim_len_y = ceil(len_y/4);
            
            range_depth_box_x = ceil((center_x-trim_len_x):(center_x+trim_len_x));
            range_depth_box_y = ceil((center_y-trim_len_y):(center_y+trim_len_y));
            depth_box = depth_im{t}(range_depth_box_x, range_depth_box_y);
            depth_box_avg = mean(depth_box(:));

            center_z = depth_box_avg;
            centers{t}(d,:) = [center_x, center_y, center_z];
            
            
            projected_x_range = round((x1:x2)+u_avg_box);
            projected_x_range(projected_x_range > max_x_dim) = [];
            projected_x_range(projected_x_range < 1) = [];
            
            projected_y_range = round((y1:y2)+v_avg_box);
            projected_y_range(projected_y_range > max_y_dim) = [];
            projected_y_range(projected_y_range < 1) = [];

            projected_depth_box = depth_im{t+1}(projected_x_range, projected_y_range);
            projected_center_z = mean(projected_depth_box(:));

            projected_center_x = center_x + u_avg_box;
            projected_center_y = center_y + v_avg_box;
%             projected_center_z = projected_depth_box_avg;
            projected_centers{t}(d,:) = [projected_center_x, projected_center_y, projected_center_z];
            
        end
    else % eval centers for last frame (no OF there)
        n_detections = size(boxes{t},1);
        centers{t} = nan(n_detections ,3);
        for d=1:n_detections
            x1 = boxes{t}(d,1);
            x2 = boxes{t}(d,2);
            y1 = boxes{t}(d,3);
            y2 = boxes{t}(d,4);

            depth_box = depth_im{t}(x1:x2, y1:y2);
            depth_box_avg = mean(depth_box(:));

            center_x = (x1+x2)/2;
            center_y = (y1+y2)/2;
            center_z = depth_box_avg;

            centers{t}(d,:) = [center_x, center_y, center_z];
        end
    end
    
    if false % (don't) visualize detections
        imshow(im);
        %     line([boxes{k}(:,1) boxes{k}(:,1) boxes{k}(:,2) boxes{k}(:,2) boxes{k}(:,1)]', [boxes{k}(:,3) boxes{k}(:,4) boxes{k}(:,4) boxes{k}(:,3) boxes{k}(:,3)]', 'color', 'r', 'linewidth', 3, 'linestyle', '-');
        for d=1:size(boxes{t},1)
            x1 = boxes{t}(d,1);
            x2 = boxes{t}(d,2);
            y1 = boxes{t}(d,3);
            y2 = boxes{t}(d,4);
            label = classes_names{classes{t}(d)};
            line([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]', 'color', 'r', 'linewidth', 3, 'linestyle', '-');
            text(x1, y1, label, 'Color', 'white');
        end
        drawnow
        shg
    end
    
    t=t+1; % Increase sampled frames counter.
end

%% save results
if take_top_n_detections ~= inf
    top_str = ['_top_' num2str(take_top_n_detections)];
else
    top_str = '';
end
th_str = regexprep(num2str(detection_thresh), '-', 'm');
th_str = regexprep(th_str, '\.', '_');

res_fname = ['../preprocessed_videos/' vid_name '_3D_detections_th' th_str top_str];
save(res_fname, 'vid_fname', 'boxes', 'classes', 'scores', 'classes_names', 'centers', 'projected_centers' ,'detection_thresh', 'frame_sample_interval', 'take_top_n_detections');
ppvid = load(res_fname);

function uvOF = OpticalFlowCLG_TV(im1gray, im2gray)
%%settings
settings.lambda = 2200; % the weighting of the data term
settings.pyramid_factor = 0.5;
settings.resampling_method = 'bicubic'; % the resampling method used to build pyramids and upsample the flow
settings.warps = 5; % the number of warps per level
settings.interpolation_method = 'cubic'; % the interpolation method used for warping
settings.its = 10; % the number of iterations used for minimization
settings.use_diffusion = 1; % apply a weighting factor to the regularization term (the diffusion coefficient)
settings.use_bilateral = 1; % the data term weighting: bilateral or gaussian
settings.wSize = 5; % the window's size for the data fidelity term (Lukas-Kanade)
settings.sigma_d = settings.wSize/6; % sigma for the distance gaussian of the bilateral filter
settings.sigma_r = 0.1; % sigma for the range gaussian of the bilateral filter
settings.use_ROF_texture = 0; % apply ROF texture to the images (1 yes, 0 no)
settings.ROF_texture_factor = 0.95; % ROF texture; I = I - factor*ROF(I);
show_flow = 0; % display the flow during computation

tic
if show_flow
    fig_hnd = figure('Name', 'Optical flow');
    [u, v] = coarse_to_fine(im1gray, im2gray, settings, show_flow, fig_hnd);
    close all
else
    fig_hnd = nan;
    [u, v] = coarse_to_fine(im1gray, im2gray, settings, show_flow, fig_hnd);
end
toc

uvOF(:, :, 1) = u;
uvOF(:, :, 2) = v;

function [im1gray, im2gray] = acquistionSeq(im1, im2)
im1gray=mat2gray(im1);
im2gray=mat2gray(im2);

im1gray=rgb2gray(im1gray);
im2gray=rgb2gray(im2gray);
