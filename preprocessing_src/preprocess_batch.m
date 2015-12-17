
initdirs
base_path = [datasets_base 'princeton_tracking_RGBD'];
[~, lsres] = system(['ls ' base_path '/EvaluationSet'])
videos = regexp(lsres, '\s+', 'split')

depth_method = 'eigen';
interval = 3;

for k=1:length(videos)
    vidname = videos{k};
    if strcmp(vidname,'')
        continue;
    end
    fname_gif = sprintf('%s__%s_smpinterv_%d.gif', vidname, depth_method,interval)
    ppvid = preprocess_wrapper('princeton', vidname, struct('depth_method', ...
                                                      'eigen', ...
                                                      'opflow_method', ...
                                                      'CLGTV', ...
                                                      'sample_interval', ...
                                                      interval, ...
                                                      'scale_to_resolution', ...
                                                      [240 320]));
    vis_preproc(ppvid, vidname, fname_gif);
end % for vidname = videos


