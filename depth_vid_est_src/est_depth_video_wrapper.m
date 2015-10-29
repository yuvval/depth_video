function est_depth_vid = est_depth_video(ppvid, mnFrameID, mxFrameID, depth_est_params)


ppvid = preprocess_wrapper('internal', 'approaching_toward.avi', struct('depth_method', 'fayao', 'opflow_method', 'CLGTV', 'sample_interval', 1, 'scale_to_resolution', [240 320]));

depth_est_params = struct('rho_OF', 1, 'C_OF', 1, 'rho_div', 1, 'dist_around_OF', 1);
