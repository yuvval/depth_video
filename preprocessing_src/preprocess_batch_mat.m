
initdirs
base_path = ['../videos/'];
[~, lsres] = system(['ls ' base_path '*.mat'])
videos = regexp(lsres, '\s+', 'split')

depth_method = 'eigen';
interval = 1;

for k=1:length(videos)
    vidname = regexprep(videos{k}, '.mat', '');
    vidname = regexprep(vidname, '../videos/', '');
    if strcmp(vidname,'')
        continue;
    end
    fname_gif = sprintf('%s__%s_smpinterv_%d.gif', vidname, depth_method,interval)
    ppvid = preprocess_wrapper('mat', vidname, struct('depth_method', ...
                                                      'eigen', ...
                                                      'opflow_method', ...
                                                      'CLGTV', ...
                                                      'sample_interval', ...
                                                      interval));
    vis_preproc(ppvid, vidname, fname_gif);
end % for vidname = videos


