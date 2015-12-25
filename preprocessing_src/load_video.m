function [rgb_frames, gt_depth_frames, depth_holes_mask_frames, vid_info] = load_video(prepr_params)
initdirs


load_vid_params = prepr_params.load_vid;

switch prepr_params.dataset
  case 'princeton'

    scaled_resln = take_from_struct(load_vid_params, 'scale_to_resolution', []);
    if isempty(scaled_resln)
        [rgb_frames, gt_depth_frames, ...
         ~, ~, vid_info] = ...
            get_frames_princeton(prepr_params.video_name, ...
                                 load_vid_params.sample_interval, scaled_resln);
    else
        [~, ~, rgb_frames, gt_depth_frames, ...
         vid_info] = ...
            get_frames_princeton(prepr_params.video_name, ...
                                 load_vid_params.sample_interval, scaled_resln);
    end
    depth_holes_mask_frames = [];

  case 'msrv3d'
    dataDir = prepr_params.video_name.dir;
    clipIndex = prepr_params.video_name.idx;
    [rgb_frames, gt_depth_frames, ...
     depth_holes_mask_frames, vid_info] = ...
        get_frames_msrv3d(dataDir, clipIndex, load_vid_params.sample_interval);

  case 'internal'
    vid_fname = fullfile(proj_root_path, 'videos', prepr_params.video_name);
    obj = VideoReader(vid_fname);
    video = obj.read();
    rgb_frames = video(:,:,:, 1:load_vid_params.sample_interval:end);
    gt_depth_frames = [];
    vid_info = [];
    depth_holes_mask_frames = [];

  case 'mat'
    vid_fname = fullfile(proj_root_path, 'videos', prepr_params.video_name);
    load(vid_fname)
    if ~exist('video');
        video = frames;
    end
    rgb_frames = video(:,:,:, 1:load_vid_params.sample_interval:end);
    gt_depth_frames = [];
    depth_holes_mask_frames = [];
    vid_info = [];
    
  otherwise
    error('Unknown dataset')
end
