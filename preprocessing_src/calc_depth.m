function [depth_frames, depth_sp] = calc_depth(video, prepr_params)
%
%

vars = {'depth_frames'};
fullname = generate_preproc_fname('depth', prepr_params);

[do_stage, depth_frames] = cond_load(fullname, 0, vars{1: end});
if ~do_stage, return, end;

num_frames  = size(video,4);
proj_root_path = get_dirs('proj_root_path');
[depth_sp, superpxl_frames] = deal( cell(num_frames,1));


switch prepr_params.depth.depth_method
  case 'fayao'
    depth_frames = nan(size(video,1), size(video,2), num_frames);
    for t = 1:num_frames
        [depth_sp{t}, depth_frames(:,:,t), superpxl_frames{t}] = ...
            getDepth_Fayao(video(:,:,:,t)); 
    end

  case 'eigen'
    frames_cell = cell(1, num_frames);
    for t = 1:num_frames
        frames_cell{t} = video(:,:,:,t);
    end
    system(['source ' proj_root_path 'theano-env/bin/activate']);
    
    gpunum = take_from_struct(prepr_params, 'gpunum', 0);
    pythoncmd = sprintf('THEANO_FLAGS=floatX=float32,device=gpu%d python', ...
        gpunum);
    funcargs.frames = frames_cell;
    res = matpyfs('infer_depth_and_normals_frames_seq', funcargs, ...
                  'eigen_depth', [proj_root_path 'preprocessing_src'], ...
                  pythoncmd);
    depth_frames = shiftdim(res.arg0, 2);
    
    % Load precalculated super pixels 
    [sp_frames, num_sps] = calc_sp(video, prepr_params);

    % For each frame, assign mean depth to each superpxl.
    for t = 1:num_frames
        depth_frame = depth_frames(:,:,t);
        sp_frame = double(sp_frames{t});
        depth_frame = imresize(depth_frame, size(sp_frame));            
        depth_sp{t} = nan(num_sps{t},1);            

        % eval mean depth per superpixel - TEMPORARLY COMMENTED OUT
        % for k = 1:num_sps{t} % NOT EFFICIENT TODO - IMPLEMENT FASTER
        %     depth_sp{t}(k) = mean(depth_frame(sp_frame==k));
        % end
        % assert(all(~isnan(depth_sp{t}(:))));

        % Started to code faster implementations. TODO!!
        % integ_depths = cumsum(depth_frame(:));
        % sample_integ_ids = diff(sp_frame(:), numel(sp_frame));
        % sampled_integ_depths = integ_depths(sample_integ_ids);
        % assert(all(~isnan(depth_sp{t}(:))));
            

    end

    otherwise, 
      error('Invalid depth method [%s]\n', prepr_params.depth_method);

 end


save(fullname, vars{1:end});