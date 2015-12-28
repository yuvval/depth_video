function [fullname, fname] = generate_preproc_fname(ppstage, prepr_params)
%
params_for_fname = prepr_params.(ppstage);

if isfield(prepr_params.load_vid, 'scale_to_resolution')
    params_for_fname.res = prepr_params.load_vid.scale_to_resolution;
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

fullpath = fullfile(get_dirs('preproc_base'), dataset, ppstage);
system(['mkdir -p ' fullpath]); % creates the dir 
%TODO Check that system command worked
fullname = fullfile(fullpath, fname);