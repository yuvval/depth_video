addpath ../
initdirs


vidname = 'front_pick_bin1'
interval = 1;
    ppvid = preprocess_wrapper('mat', vidname, struct('depth_method', ...
                                                      'eigen', ...
                                                      'opflow_method', ...
                                                      'CLGTV', ...
                                                      'sample_interval', ...
                                                      interval));

depth_est_params = struct('rho_OF', 1, 'C_OF', 1, 'rho_div', 1, 'dist_around_OF', 1);
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
tacktile_ind = zero(2,2,N);
tacktile_ind(:,:,1) = eye(2);
tacktile_ind(:,:,2) = eye(2);
tacktile_ind(:,:,3) = eye(2);
tacktile_ind(:,:,4) = eye(2);
tacktile_ind(1,2,5) = 1;
tacktile_ind(1,2,6) = 1;
tacktile_ind(2,1,7) = 1;

%%
depth_est_params.boxes = boxes;
depth_est_params.tacktile_ind = tacktile_ind;
dvid = est_depth_video(ppvid, mnFrameID, mxFrameID, ...
                                depth_est_params);

save(sprintf('~/tmp/%s_est_d.mat', vidname), 'dvid', 'mnFrameID',...
'mxFrameID', 'interval', 'vidname', 'depth_est_params')

fname_gif = sprintf('~/tmp/%s_est_d.gif', vidname);

%%
vis_dep_est(ppvid, vidname, dvid, mnFrameID, ...
                     mxFrameID, fname_gif, true, boxes)
