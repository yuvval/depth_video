function [depth_frames, opflow_uvframes, superpxl_frames] = ...
    preprocess_frames(rgb_frames, prepr_params)

% frame_sample_interval = prepr_params.sample_interval;
video = rgb_frames;

%% Iterate on frames.
% [frames, depths_superpxl, depths_pxl, superpxl_frames, edges, opflow_uvframes, Gx, Gy] = deal({});

[opflow_uvframes, superpxl_frames, depths_superpxl, depths_pxl] = deal({});
Nframes = size(video,4);

t = 1; % Sampled frames counter.

im=video(:,:,:,1);

switch prepr_params.depth_method
    case 'fayao'
        [depths_superpxl{1}, depths_pxl{1}, superpxl_frames{1}] ...
            = getDepth_Fayao(im);
    otherwise
        error('Unknown depth method');
end

% imgray = double(rgb2gray(im));
imgray=mat2gray(im);  			
imgray=rgb2gray(imgray);  		

imsize = size(imgray);
depth_frames = nan(imsize(1), imsize(2), Nframes);
depth_frames(:,:,1) = depths_pxl{1};

for k=1:Nframes
    frames{t} = im;
    
    

    if k <= (Nframes-1) % we can't eval OF for the last frame, so we skip it.
        %% Evaluate depth estimation
        im2 = video (:,:,:,k+1);
%         im2gray = double(rgb2gray(im2));
        im2gray=mat2gray(im2);
        im2gray=rgb2gray(im2gray);
        
        switch prepr_params.depth_method
            case 'fayao'
                [depths_superpxl{t+1}, depths_pxl{t+1}, ...
                    superpxl_frames{t+1}] = getDepth_Fayao(im2);
            otherwise
                error('Unknown depth method');
        end
        depth_frames(:,:,t+1) = depths_pxl{t+1};
        
        %% Optical Flow
        switch prepr_params.opflow_method
            case 'CLGTV'
                opflow_uvframes{t} = OpticalFlowCLG_TV(imgray, im2gray);
            otherwise
                error('Unknown optical flow method.');
        end
        
        
        
        % preparing for next iteration
        im = im2;
        imgray = im2gray;
    end
    
    t=t+1;
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

