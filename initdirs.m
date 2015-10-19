
switch get_domain()
    case 'cortex'
        homedir = '~';
        datasets_base = '/cortex/data/video';
        preproc_base = '/cortex/users/yuvval/depth_vl/preproc/';
    case 'csail'
        homedir = '~';
        datasets_base = '~/data/video';
        preproc_base = '~/depth_vl/preproc/';
end

sep = filesep; % equals '/' in linux

if ~exist('proj_root_path', 'var')
    proj_root_path = [homedir '/depth_video/'];
end

addpath([proj_root_path sep 'lowlevel_vision_procedures/']);
addpath([proj_root_path sep 'packages/']);
addpath([proj_root_path sep 'optical_flow/']);
addpath([proj_root_path sep 'optical_flow/algorithms/CLG-TV//']);
addpath([proj_root_path sep 'preprocessing_src//']);

%%
mcg_root = [homedir '/externals/mcg-2.0/pre-trained'];
addpath(mcg_root);
addpath(fullfile(mcg_root,'lib'));
addpath(fullfile(mcg_root,'scripts'));
addpath(fullfile(mcg_root,'datasets'));
addpath(genpath(fullfile(mcg_root,'src')));
addpath(genpath(fullfile(mcg_root,'src', 'bboxes')));
