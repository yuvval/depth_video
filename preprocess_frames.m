function [ppvid, res_fname ] = preprocess_video(vid_fname, frame_sample_interval, trim_first_seconds)
%% init
if nargin<1
    vid_fname = 'videos/approaching_toward.avi'; 
end

if nargin < 2
    frame_sample_interval = 15; % Sample a frame from video once every X frames.
end

if nargin < 3
    trim_first_seconds = 3;
end

initdirs

fname_split = regexp(vid_fname, '[\./]', 'split');
vid_name = fname_split{end-1};

obj = VideoReader(vid_fname);
video = obj.read();

frames_per_sec = 30;
%% Iterate on frames.
[frames, depths_superpxl, depths_pxl, edges, uvOFs] = deal({});

Nframes = size(video,4);

t = 1; % Sampled frames counter.
[depths_superpxl{1}, depths_pxl{1}] = getDepth_Fayao(video(:,:,:,1));
for k=1 + (trim_first_seconds*frames_per_sec):frame_sample_interval:Nframes
    im=video(:,:,:,k);
    frames{t} = im;



    %% Optical Flow
    if k <= (Nframes-frame_sample_interval) % we can't eval OF for the last frame, so we skip it.
        %% Evaluate depth estimation
        im2 = video (:,:,:,k+frame_sample_interval);
        tic
        [depths_superpxl{t+1}, depths_pxl{t+1}] = getDepth_Fayao(im2);
        toc

        [im1gray, im2gray] = OFacquistionSeq(im, im2);
        
        uvOFs{t} = OpticalFlowCLG_TV(im1gray, im2gray);
        
    end
    
    t=t+1;
end
%% save results
fsmp_str = [ '_fsmp_' num2str(frame_sample_interval)];
res_fname = ['preprocessed_videos/' vid_name fsmp_str '_ppvid'];
ppvid = v2struct(vid_fname, frames, depths_superpxl, depths_pxl, edges, uvOFs, frame_sample_interval, trim_first_seconds);
fname_depth_gif = visualize_preproc(vid_name, ppvid);
save(res_fname, 'ppvid', '-v7.3');


function [fname_depth_gif, fname_frames_gif] = visualize_preproc(vid_name, ppvid)
N = length(ppvid.depths_pxl);
fsmp_str = [ '_fsmp_' num2str(ppvid.frame_sample_interval)];

fname_depth_gif = sprintf('depth_%s%s.gif', vid_name, fsmp_str);
fname_frames_gif = sprintf('frames_%s%s.gif', vid_name, fsmp_str);

max_depth = max(max(cell2mat(ppvid.depths_pxl)));
for t=1:N
    figure(1)
    imshow(ppvid.depths_pxl{t}/max_depth);
    colormap(flipud(parula));shg
    colorbar;shg
    drawnow
    
    save_animated_gif_frame(fname_depth_gif, t==1);
    if t==N % saving last frame again
        save_animated_gif_frame(fname_depth_gif, t==1);
        save_animated_gif_frame(fname_depth_gif, t==1);
    end

    figure(2)
    imshow(ppvid.frames{t});
    drawnow
    save_animated_gif_frame(fname_frames_gif, t==1);
    if t==N % saving last frame again
        save_animated_gif_frame(fname_frames_gif, t==1);
        save_animated_gif_frame(fname_frames_gif, t==1);
    end

end







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

function [im1gray, im2gray] = OFacquistionSeq(im1, im2)
im1gray=mat2gray(im1);
im2gray=mat2gray(im2);

im1gray=rgb2gray(im1gray);
im2gray=rgb2gray(im2gray);
