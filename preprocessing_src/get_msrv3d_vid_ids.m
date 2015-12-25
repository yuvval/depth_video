function [vid_ids, dir_names, basedir] = get_msrv3d_vid_ids(basedir)
if nargin<1
    basedir = '/storage/yuvval/data/video/';
end

lsres = ls([basedir 'MSRV3D/']);
ls_raw_list = regexp(lsres, '\s', 'split');

dir_names = {};
cnt = 1;
for k = 1:length(ls_raw_list)
    if is_in_str('Indoor', ls_raw_list{k});
        dir_names{cnt} = ls_raw_list{k};
        cnt = cnt+1;
    end
end


cnt = 1;
for k = 1:length(dir_names)
    dir = dir_names{k};
    info = load([basedir 'MSRV3D/' dir '/info.mat']);
    num_vids_in_dir = length(info.clips);
    
    for m = 1:num_vids_in_dir
        vid_id.dir = dir;
        vid_id.idx = m;
        vid_ids{cnt} = vid_id;
        cnt = cnt+1;
    end
end

