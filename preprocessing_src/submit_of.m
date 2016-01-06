clear
addpath ../
initdirs

[vid_ids, ~, basedir]  = get_msrv3d_vid_ids();

baseparams.opflow.opflow_method = 'CLGTV';
baseparams.dataset = 'msrv3d';
baseparams.mnFrameID = 1;
baseparams.mxFrameID = -1; % -1 means last frame

%start_parallel_pool(45, '/tmp/yuvval/');

baseparams.sample_interval = 1; % 30 fps
% parfor k=1:length(vid_ids)
for k=1:length(vid_ids)
    params = baseparams;
    params.video_name = vid_ids{k};
    video = load_video(params);
    calc_optical_flow(video, params);
end



baseparams.sample_interval = 6; % 5fps
for k=1:length(vid_ids)
% parfor k=1:length(vid_ids)
    params = baseparams;
    params.video_name = vid_ids{k};
    video = load_video(params);
    calc_optical_flow(video, params);
end


