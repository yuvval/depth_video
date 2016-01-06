function [rgb_frames, gt_depth_frames, depth_holes_mask_frames, ...
          vid_info] = load_video(prepr_params)
%
% [rgb_frames, gt_depth_frames, depth_holes_mask_frames, ...
%              vid_info] = load_video(prepr_params)
%


fullfname = generate_preproc_fname('load_vid', prepr_params);
do_force = take_from_struct(prepr_params, 'do_force', 0);
preproc_vars = split_csv(['rgb_frames, gt_depth_frames, ' ...
                    'depth_holes_mask_frames, vid_info']);
[do_stage, rgb_frames, gt_depth_frames, depth_holes_mask_frames, ...
 vid_info] = cond_load(fullfname, do_force, preproc_vars{1:end});

if ~do_stage, fprintf('Loaded video from %s\n', fullfname);return;end


if ~isfield (prepr_params, 'load_vid')
    prepr_params.load_vid = struct();
end

load_vid_params = prepr_params.load_vid;

switch prepr_params.dataset
  case 'princeton'

    scaled_resln = take_from_struct(load_vid_params, 'scale_to_resolution', []);
    if isempty(scaled_resln)
        [rgb_frames, gt_depth_frames, ...
         ~, ~, vid_info] = ...
            get_frames_princeton(prepr_params.video_name, ...
                                 1, scaled_resln);
    else
        [~, ~, rgb_frames, gt_depth_frames, ...
         vid_info] = ...
            get_frames_princeton(prepr_params.video_name, ...
                                 1, scaled_resln);
    end
    depth_holes_mask_frames = [];

  case 'msrv3d'
    dataDir = prepr_params.video_name.dir;
    clipIndex = prepr_params.video_name.idx;
    [rgb_frames, gt_depth_frames, ...
     depth_holes_mask_frames, vid_info] = ...
        get_frames_msrv3d(dataDir, clipIndex, 1);

  case 'internal'
    vid_fname = fullfile(proj_root_path, 'videos', prepr_params.video_name);
    obj = VideoReader(vid_fname);
    video = obj.read();
    rgb_frames = video;
    gt_depth_frames = [];
    vid_info = [];
    depth_holes_mask_frames = [];

  case 'mat'
    vid_fname = fullfile(proj_root_path, 'videos', prepr_params.video_name);
    load(vid_fname)
    if ~exist('video', 'var');
        video = frames;
    end
    rgb_frames = video;
    gt_depth_frames = [];
    depth_holes_mask_frames = [];
    vid_info = [];
    
  otherwise
    error('Unknown dataset')
end

rgb_frames = cast_video_to_double(rgb_frames);

save(fullfname, preproc_vars{1:end}, '-v7.3');


function rgb_frames = cast_video_to_double(rgb_frames)
% 
if ~isa(rgb_frames, 'uint8')
    return
end
rgb_frames = rgb_frames/255;
