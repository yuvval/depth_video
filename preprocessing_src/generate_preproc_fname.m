function [fullname, fname] = generate_preproc_fname(ppstage, prepr_params)
% function [fullname, fname] = generate_preproc_fname(ppstage, prepr_params)


if ~isfield (prepr_params, 'load_vid')
    prepr_params.('load_vid') = struct();
end


% params_for_fname struct holds the parameters to derive the filename
if strcmp(ppstage, 'sp_opflow_overlap')
    params_for_fname = catstruct(prepr_params.opflow, prepr_params.superpixels);
else 
    params_for_fname = prepr_params.(ppstage); 
end

if isfield(prepr_params.load_vid, 'scale_to_resolution')
    params_for_fname.scale_to_resolution = ...
        prepr_params.load_vid.scale_to_resolution;

    %shorten field name 
    params_for_fname = rename_struct_field(params_for_fname, 'scale_to_resolution', 'res');
end

% optical flow depends on the video sampling interval (between frames)
if strcmp(ppstage, 'opflow') || strcmp(ppstage, 'sp_opflow_overlap')
    if isfield(prepr_params, 'sample_interval')
        params_for_fname.smp_ivl = prepr_params.sample_interval;
        params_for_fname.mnF = prepr_params.mnFrameID;
        params_for_fname.mxF = prepr_params.mxFrameID;
    end
end


dataset = prepr_params.dataset;
tmp_name = prepr_params.video_name;

if strcmp(dataset, 'msrv3d') % msrv3d has a unique way of video naming
    video_name = sprintf('%s_%d', tmp_name.dir, tmp_name.idx);
else
    video_name = tmp_name;
end % if strcmp

fname = regexprep(sprintf('%s_%s_%s', dataset, video_name, ...
buildStringFromStruct(params_for_fname, '__')), '\.', '_' );
fname = [fname '.mat']; % add extension

fullpath = fullfile(get_dirs('preproc_base'), dataset, ppstage);
system(['mkdir -p ' fullpath]); % creates the dir 
%TODO Check that system command worked
fullname = fullfile(fullpath, fname);


