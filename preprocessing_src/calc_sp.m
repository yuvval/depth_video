function [sp_frames, num_sps] = calc_sp(video, prepr_params)
%
vars = {'sp_frames', 'num_sps'};
fullname = generate_preproc_fname('superpixels', prepr_params);

[do_stage, sp_frames, num_sps] = cond_load(fullname, 0, vars{1: end});
if ~do_stage, return, end;


num_frames = size(video, 4);
switch prepr_params.superpixels.sp_method
    case 'slico',
      for t=1:num_frames
          im = video(:,:,:,t);
          [sp_frames{t}, num_sps{t}] = slicomex(int8(im*255), prepr_params.superpixels.n_superpixels);
          sp_frames{t} = sp_frames{t}+1; % start labeling sp from 1
                                         % (instead of 0)
      end
  otherwise, 
    error('Invalid sp_method = [%s]\n', prepr_params.sp_method);
end



save(fullname, vars{1:end});