
csail_home = getenv('CSHOME');
ctx_home = getenv('CTXHOME');

homedir = csail_home;
if strcmp(homedir, '' )
    homedir = '~';
end
sep = filesep; % equals '/' in linux

addpath packages/
addpath optical_flow
addpath optical_flow/algorithms/CLG-TV/ 


%%
mcg_root = [homedir '/externals/mcg-2.0/pre-trained'];
addpath(mcg_root);
addpath(fullfile(mcg_root,'lib'));
addpath(fullfile(mcg_root,'scripts'));
addpath(fullfile(mcg_root,'datasets'));
addpath(genpath(fullfile(mcg_root,'src')));
addpath(genpath(fullfile(mcg_root,'src', 'bboxes')));
