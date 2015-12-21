clear
addpath ../
initdirs


vidname = 'front_pick_bin1'
interval = 1;
ppvid = preprocess_wrapper('mat', vidname, struct('depth_method', ...
    'eigen', ...
    'opflow_method', ...
    'CLGTV', ...
    'sample_interval', ...
    interval, ...
    'superpixels', ...
    'slico', ...
    'n_superpixels', ...
    2000));

depth_est_params = struct('rho_OF', 1, 'C_OF', 1, 'rho_div', 1, 'dist_around_OF', 1);

% t = 5
% imagesc(superpxl_intensities_to_img(ppvid.superpxl_frames{t}, ppvid.depths_superpxl{t}, ppvid.n_superpxl{t}))


mnFrameID = 1;
mxFrameID = -1; % -1 means last frame

% vidname = 'br_occ_turn0'
% interval = 3;
%     ppvid = preprocess_wrapper('princeton', vidname, struct('depth_method', ...
%                                                       'eigen', ...
%                                                       'opflow_method', ...
%                                                       'CLGTV', ...
%                                                       'sample_interval', ...
%                                                       interval, ...
%                                                       'scale_to_resolution', ...
%                                                       [240 320]));
% mnFrameID = 1;
% mxFrameID = 30; % -1 means till last frame
% depth_est_params = struct('rho_OF', 1, 'C_OF', 1, 'rho_div', 1, 'dist_around_OF', 1);

%%
boxes{1} = [ 152 186 51 147 ;
    152 190 48 165
    152 190 60 175
    152 190 59 165
    142 190 79 165
    152 190 55 200
    152 190 55 200
    ];
boxes{2} = [ 140 181 165 207 ;
    140 181 165 207 ;
    140 181 165 207 ;
    140 181 165 207 ;
    140 181 165 207 ;
    140 186 93 130 ;
    140 186 93 130 ;
    ];
N = length(boxes);
tactile_inds = zeros(2,2,N);
tactile_inds(:,:,1) = eye(2); % no tactile cue
tactile_inds(:,:,2) = eye(2); % no tactile cue
tactile_inds(:,:,3) = eye(2); % no tactile cue
tactile_inds(:,:,4) = eye(2); % no tactile cue
tactile_inds(1,2,5) = 1; % tactile indication masks 1:2
tactile_inds(1,2,6) = 1; % tactile indication masks 1:2
tactile_inds(2,1,7) = 1; % tactile indication masks 1:2

% generate masks from boxes
masks = {};
for t = 1:size(boxes{1},1)
    for b = 1:length(boxes)
        masks{t}{b} = zeros(240,320);
        masks{t}{b}(boxes{b}(t, 3):boxes{b}(t, 4), boxes{b}(t, 1):boxes{b}(t, 2)) = 1;
        masks{t}{b} = sparse(imresize(masks{t}{b}, [109, 147])>=0.5); % resizing according to eigen depth estimation imsize
    end
end

%%
depth_est_params.C_intra = 1;
depth_est_params.rho_intra = 1;

depth_est_params.C_OF = 1;

% depth_est_params.detect_masks = masks;
depth_est_params.tactile_inds = tactile_inds;
depth_est_params.C_tact = 5;
depth_est_params.rho_tact = 10;
dvid_sp = est_depth_video_sp(ppvid, mnFrameID, mxFrameID, ...
    depth_est_params);


dvid = nan(240,320,length(dvid_sp));
for t=1:length(dvid_sp)
    dvid(:,:,t) = superpxl_intensities_to_img(ppvid.superpxl_frames{t}, dvid_sp{t}, ppvid.n_superpxl{t});
end
save(sprintf('~/tmp/%s_est_sp_d.mat', vidname), 'dvid_sp', 'mnFrameID',...
    'mxFrameID', 'interval', 'vidname', 'depth_est_params')

fname_gif = sprintf('~/www/figs/%s_est_d_intra_C_%f_std_%f.gif', vidname, ...
                    depth_est_params.C_intra, depth_est_params.rho_intra);

%%
vis_dep_est(ppvid, vidname, dvid, mnFrameID, ...
    mxFrameID, fname_gif, false)
%     mxFrameID, fname_gif, true, boxes)
scp_file(fname_gif);
