function [rgb_frames, depth_frames, depth_holes_frames, clip_info] = ...
    get_frames_msrv3d(dataDir, clipIndex, sample_interval,  base_path)


if nargin <3
    sample_interval = 1;
end

if nargin <4
    base_path = [get_dirs('datasets_base') '/MSRV3D/'];
end

addpath(base_path);


[rgb_frames_all, ~, ~, ~, holesMask_all, depth_frames_all, clip_info] = ...
    loadStereoRGBDdata([base_path dataDir], clipIndex, [431 320]);

max_frames = size(rgb_frames_all, 4);
frames_range = 1:sample_interval:max_frames;

rgb_frames = (rgb_frames_all(95:334,1:320,:, frames_range));
depth_frames = (depth_frames_all(95:334,1:320,:, frames_range));
depth_holes_frames = (holesMask_all(95:334,1:320,:, frames_range));
