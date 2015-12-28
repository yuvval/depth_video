
homedir = get_dirs('homedir');
if ~exist('proj_root_path', 'var')
    proj_root_path = [homedir '/depth_video/'];
end

addpath([proj_root_path '/lowlevel_vision_procedures/']);
addpath([proj_root_path '/packages/']);
addpath([proj_root_path '/packages/SLIC_mex/']);
addpath([proj_root_path '/packages/matpyfs']);
addpath([proj_root_path '/optical_flow/']);
addpath([proj_root_path '/optical_flow/algorithms/CLG-TV/']);
addpath([proj_root_path '/preprocessing_src/']);
addpath([proj_root_path '/depth_vid_est_src/']);

%%
addpath([homedir '/externals/matpyfs']);

%%
mcg_root = [homedir '/externals/mcg-2.0/pre-trained'];
addpath(mcg_root);
addpath(fullfile(mcg_root,'lib'));
addpath(fullfile(mcg_root,'scripts'));
addpath(fullfile(mcg_root,'datasets'));
addpath(genpath(fullfile(mcg_root,'src')));
addpath(genpath(fullfile(mcg_root,'src', 'bboxes')));
