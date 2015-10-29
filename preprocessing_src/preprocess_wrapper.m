function ppvid = preprocess_wrapper(dataset, video_name, prepr_params)
initdirs


fname = generate_preproc_fname(dataset, video_name, prepr_params);
fullpath = fullfile(preproc_base,...
    prepr_params.depth_method, prepr_params.opflow_method, dataset);
system(['mkdir -p ' fullpath]); % creates the dir if it doesn't exists

fullfname = fullfile(fullpath, fname);

do_force = 0;
preproc_vars = ...
    {'rgb_frames', 'depth_frames', 'opflow_frames' 'superpxl_frames', ....
    'gt_depth_frames', 'camera_info'};
[do_preproc, ...
    rgb_frames, depth_frames, opflow_frames, superpxl_frames, ...
    gt_depth_frames, camera_info] ...
    = cond_load([fullfname '.mat'], do_force, preproc_vars{1:end});

if do_preproc
    
    switch dataset
        case 'princeton'

            scaled_resln = take_from_struct(prepr_params, 'scale_to_resolution', []);
            if isempty(scaled_resln)
                [rgb_frames, gt_depth_frames, ...
                    ~, ~, camera_info] = ...
                    get_frames_princeton(video_name, ...
                    prepr_params.sample_interval, scaled_resln);
            else
                [~, ~, rgb_frames, gt_depth_frames, ...
                    camera_info] = ...
                    get_frames_princeton(video_name, ...
                    prepr_params.sample_interval, scaled_resln);
            end
            
            
            
        otherwise
            error('Unknown dataset')
    end
    [depth_frames, opflow_frames, superpxl_frames] = ...
        preprocess_frames(rgb_frames, prepr_params);
    
    preproc_vars{end+1} = 'fieldNames';
    ppvid = v2struct(preproc_vars);
    save(fullfname, '-struct', 'ppvid', '-v7.3');
else
    preproc_vars{end+1} = 'fieldNames';
    ppvid = v2struct(preproc_vars);
end



% def preprocess_wrapper(dataset, video_name, pp_params)
% 
% 
% ppfname = generate_pp_fname()
% condload(ppfname)
% 
% if do_..
% 
%   switch dataset
%     case ...
%        reso = take_from_struct(pp..res, [])
%        get video + depthGT(reso)
%   
%   save (ppfname, ppvid)
