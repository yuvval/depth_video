function ppvid = preprocess_wrapper(prepr_params)
% function ppvid = preprocess_wrapper(prepr_params)
% join preprocessed (mostly precomputed) data to a single struct

[ppvid.rgb_frames, ppvid.gt_depth_frames, ...
 ppvid.depth_holes_mask_frames, ppvid.vid_info] = ...
    load_video(prepr_params);

[ppvid.depth_frames, ppvid.depths_superpxl, ...
 ppvid.superpxl_frames, ppvid.n_superpxl] = ...
    calc_depth(ppvid.rgb_frames, prepr_params);

ppvid.opflow_frames = calc_optical_flow(ppvid.rgb_frames, prepr_params);
ppvid.prepr_params = prepr_params;
