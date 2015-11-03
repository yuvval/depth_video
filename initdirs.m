switch get_domain()
    case 'cortex'
        homedir = '/home/lab/yuvval';
        datasets_base = '/cortex/data/video/';
        preproc_base = '/cortex/users/yuvval/depth_vl/preproc/';
    case 'csail'
        homedir = '/afs/csail.mit.edu/u/y/yuvval';
        datasets_base = [homedir '/data/video/'];
        preproc_base = [homedir '/depth_vl/preproc/'];
end

sep = filesep; % equals '/' in linux

if ~exist('proj_root_path', 'var')
    proj_root_path = [homedir '/depth_video/'];
end

addpath([proj_root_path sep 'lowlevel_vision_procedures/']);
addpath([proj_root_path sep 'packages/']);
addpath([proj_root_path sep 'packages/matpyfs']);
addpath([proj_root_path sep 'optical_flow/']);
addpath([proj_root_path sep 'optical_flow/algorithms/CLG-TV/']);
addpath([proj_root_path sep 'preprocessing_src/']);
addpath([proj_root_path sep 'depth_vid_est_src/']);

%%
addpath([homedir sep 'externals/matpyfs']);

%%
mcg_root = [homedir '/externals/mcg-2.0/pre-trained'];
addpath(mcg_root);
addpath(fullfile(mcg_root,'lib'));
addpath(fullfile(mcg_root,'scripts'));
addpath(fullfile(mcg_root,'datasets'));
addpath(genpath(fullfile(mcg_root,'src')));
addpath(genpath(fullfile(mcg_root,'src', 'bboxes')));
