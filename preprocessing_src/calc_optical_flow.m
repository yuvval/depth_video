function [opflow_uvframes] = calc_optical_flow(video, prepr_params)
%
% [opflow_uvframes] = calc_optical_flow(video, prepr_params)
%

fullfname = generate_preproc_fname('opflow', prepr_params);
do_force = take_from_struct(prepr_params, 'do_force', 0);
preproc_vars = {'opflow_uvframes'};
[do_stage, opflow_uvframes] = cond_load([fullfname '.mat'], do_force, ...
                                        preproc_vars{1:end});
if ~do_stage, fprintf('Loaded OF from %s\n', fullfname);return;end


ofmethod = prepr_params.opflow.opflow_method;

[opflow_uvframes] = deal({});
Nframes = size(video,4);

im=video(:,:,:,1);
imgray=mat2gray(im);  			
imgray=rgb2gray(imgray);  		

for k=1:(Nframes-1) % skip OF for the last frame 
    im2 = video (:,:,:,k+1);
    im2gray = mat2gray(im2);
    im2gray = rgb2gray(im2gray);
    
    %% Optical Flow
    switch ofmethod
      case 'CLGTV'
        opflow_uvframes{k} = OpticalFlowCLG_TV(imgray, im2gray, prepr_params.opflow);
      otherwise
        error('Unknown optical flow method [%s]\n', ofmethod);
    end
    
    % preparing for next iteration
    imgray = im2gray;
end
save(fullfname, preproc_vars{1:end}, '-v7.3');

function uvOF = OpticalFlowCLG_TV(im1gray, im2gray, settings)
[~, settings] = take_from_struct(settings, 'lambda', 2200); % the
                                                            % weighting
                                                            % of
                                                            % the
                                                            % data
                                                            % term
[~, settings] = take_from_struct(settings, 'pyramid_factor', 0.5);
[~, settings] = take_from_struct(settings, 'resampling_method', ...
                                 'bicubic'); % the resampling
                                             % method
                                             % used to build
                                             % pyramids
                                             % and upsample the
                                             % flow
[~, settings] = take_from_struct(settings, 'warps', 5); % the
                                                        % number of
                                                        % warps per
                                                        % level
[~, settings] = take_from_struct(settings, 'interpolation_method', ...
                                 'cubic'); % the interpolation
                                           % method
                                           % used for warping
[~, settings] = take_from_struct(settings, 'its', 10); % the number
                                                       % of
                                                       % iterations
                                                       % used for
                                                       % minimization
[~, settings] = take_from_struct(settings, 'use_diffusion', 1); % apply a weighting factor to the
                                                                % regularization
                                                                % term
                                                                % (the
                                                                % diffusion
                                                                % coefficient)
[~, settings] = take_from_struct(settings, 'use_bilateral', 1); % the data term weighting: bilateral or
                                                                % gaussian
[~, settings] = take_from_struct(settings, 'wSize', 5); % the
                                                        % window's
                                                        % size for
                                                        % the data
                                                        % fidelity
                                                        % term
                                                        % (Lukas-Kanade)
[~, settings] = take_from_struct(settings, 'sigma_d', settings.wSize/6); ...
% sigma for the distance
% gaussian of the bilateral
% filter
[~, settings] = take_from_struct(settings, 'sigma_r', 0.1); % sigma
                                                            % for
                                                            % the
                                                            % range
                                                            % gaussian
                                                            % of
                                                            % the
                                                                    
% bilateral filter
[~, settings] = take_from_struct(settings, 'use_ROF_texture', 0); % apply ROF texture to the images (1
                                                                  % yes,
                                                                  % 0
                                                                  % no)
[~, settings] = take_from_struct(settings, 'ROF_texture_factor', ...
                                 0.95); % ROF texture; I', I -
                                        % factor*ROF(I);


show_flow = 0; % don't display the flow during computation
fig_hnd = nan;
[u, v] = coarse_to_fine(im1gray, im2gray, settings, show_flow, ...
                            fig_hnd);


uvOF = cat(3, u, v);
