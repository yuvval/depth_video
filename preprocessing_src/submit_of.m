addpath ../
initdirs

[vid_ids, ~, basedir]  = get_msrv3d_vid_ids();

baseparams.opflow.opflow_method = 'CLGTV';
baseparams.dataset = 'msrv3d';
baseparams.load_vid.sample_interval = 1;

start_parallel_pool(45, '/tmp/yuvval/');
parfor k=1:length(vid_ids)
    params = baseparams;
    params.video_name = vid_ids{k};
    rgb_frames = load_video(params);
    calc_optical_flow(rgb_frames, params);
end



