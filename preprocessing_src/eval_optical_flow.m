function [opflow_uvframes] = eval_optical_flow(rgb_frames, ...
                                               prepr_params)

initdirs

dataset = prepr_params.dataset;
ofmethod = prepr_params.opflow.opflow_method;


fname = generate_preproc_fname('opflow', prepr_params);
fullpath = fullfile(preproc_base, 'OpF', ofmethod, dataset);
system(['mkdir -p ' fullpath]); % creates the dir if it doesn't exists

fullfname = fullfile(fullpath, fname);

do_force = take_from_struct(prepr_params, 'do_force', 0);

preproc_vars = {'opflow_uvframes'};
[do_stage, opflow_uvframes] ...
    = cond_load([fullfname '.mat'], do_force, preproc_vars{1:end});

if do_stage
    video = rgb_frames;

    [opflow_uvframes] = deal({});
    Nframes = size(video,4);

    im=video(:,:,:,1);

    n255 = 1; 
    if max(im(:)) > 1
        n255 = 255;
    end

    im = im/n255;
    imgray=mat2gray(im);  			
    imgray=rgb2gray(imgray);  		

    for k=1:(Nframes-1) % we can't eval OF for the last frame, so we skip it.
        %% Evaluate depth estimation
        im2 = video (:,:,:,k+1)/n255;

        im2gray=mat2gray(im2);
        im2gray=rgb2gray(im2gray);
        
        %% Optical Flow
        switch ofmethod
          case 'CLGTV'
            opflow_uvframes{k} = OpticalFlowCLG_TV(imgray, im2gray);
          otherwise
            error('Unknown optical flow method.');
        end
        
        % preparing for next iteration
        im = im2;
        imgray = im2gray;
        
    end

    preproc_vars{end+1} = 'fieldNames';
    ofvid = v2struct(preproc_vars);
    save(fullfname, '-struct', 'ofvid', '-v7.3');
end % if do_stage




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




