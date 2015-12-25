function [depth_frames, opflow_uvframes, superpxl_frames, depths_superpxl, n_superpxl] = ...
    preprocess_frames(rgb_frames, prepr_params)

initdirs
% frame_sample_interval = prepr_params.sample_interval;
video = rgb_frames;

%% Iterate on frames.
% [frames, depths_superpxl, depths_pxl, superpxl_frames, edges, opflow_uvframes, Gx, Gy] = deal({});

[opflow_uvframes, superpxl_frames, n_superpxl, depths_superpxl, depths_pxl] = deal({});
Nframes = size(video,4);

t = 1; % Sampled frames counter.

im=video(:,:,:,1);

n255 = 1; 
if max(im(:)) > 1
    n255 = 255;
end
im = im/n255;

frames_tmp = {};
if strcmp(prepr_params.depth_method, 'eigen') && ~isfield(prepr_params, 'superpixels')
    frames_tmp{1} = im;
    im = imresize(im, [109 147]);
else
    frames_tmp{1} = im;
    [superpxl_frames{1}, n_superpxl{1}] = slicomex(int8(im*255),prepr_params.n_superpixels);
    superpxl_frames{1} = superpxl_frames{1} +1; % start labeling sp from 1 (instead of 0)
end % if strcmp


% imgray = double(rgb2gray(im));
imgray=mat2gray(im);  			
imgray=rgb2gray(imgray);  		

imsize = size(imgray);
depth_frames = nan(imsize(1), imsize(2), Nframes);

for k=1:Nframes
    frames{t} = im;
    

    if k <= (Nframes-1) % we can't eval OF for the last frame, so we skip it.
        %% Evaluate depth estimation
        im2 = video (:,:,:,k+1)/n255;
        if strcmp(prepr_params.depth_method, 'eigen')  && ~isfield(prepr_params, 'superpixels')
            frames_tmp{t+1} = im2;
            im2 = imresize(im2, [109 147]);
        else
            frames_tmp{t+1} = im2;
            [superpxl_frames{t+1}, n_superpxl{t+1}] = slicomex(int8(im2*255),prepr_params.n_superpixels);
            superpxl_frames{t+1} = superpxl_frames{t+1} +1; % start labeling sp from 1 (instead of 0)
        end % if strcmp

%         im2gray = double(rgb2gray(im2));
        im2gray=mat2gray(im2);
        im2gray=rgb2gray(im2gray);
        
        %        switch prepr_params.depth_method
        %            case 'fayao'
        %                [depths_superpxl{t+1}, depths_pxl{t+1}, ...
        %   superpxl_frames{t+1}] = getDepth_Fayao(im2);
        %            otherwise
        %error('Unknown depth method');
        %        end
        %        depth_frames(:,:,t+1) = depths_pxl{t+1}; 
        
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

% Evaluate depths

switch prepr_params.depth_method
  case 'fayao'
    for t = 1:length(frames)
        [depths_superpxl{t}, depths_pxl{t}, superpxl_frames{t}] ...
            = getDepth_Fayao(frames{t});
        depth_frames(:,:,t) = depths_pxl{t}; 
    end % for k =1:length(frames)
  case 'eigen'
    system(['source ' proj_root_path 'theano-env/bin/activate'])
    gpunum = randi(4) -1;
    pythoncmd = sprintf(['THEANO_FLAGS=floatX=float32,device=gpu%d ' ...
                        'python'], gpunum);
    funcargs.frames = frames_tmp;
    res = matpyfs('infer_depth_and_normals_frames_seq', funcargs, ...
                  'eigen_depth', [proj_root_path 'preprocessing_src'] ...
                  , pythoncmd);

    depth_frames = shiftdim(res.arg0, 2);
    
    if isfield(prepr_params, 'superpixels')
        % Iterate on all the frames and assign mean depth to each superpxl.
        for t = 1:length(frames) 
            depth_frame = depth_frames(:,:,t);
            superpxl_frame = double(superpxl_frames{t});
            depth_frame = imresize(depth_frame, size(superpxl_frame));            
            depths_superpxl{t} = nan(n_superpxl{t},1);            
            for k = 1:n_superpxl{t}
                depths_superpxl{t}(k) = mean(depth_frame(superpxl_frame==k));
            end
            assert(all(~isnan(depths_superpxl{t}(:))));
        end        
    end

  otherwise
    error('Unknown depth method');
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

