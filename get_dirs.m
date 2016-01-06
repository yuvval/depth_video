function dir = get_dirs(dir_type)
% 
domain = get_domain();
switch domain 
    case 'cortex'
      homedir = '/home/lab/yuvval';
      switch dir_type
        case 'homedir', dir = homedir;
        case 'datasets_base', dir = '/cortex/data/video/';
        case 'preproc_base', dir = ['/cortex/users/yuvval/depth_vl/' ...
                              'preproc/'];
        case 'proj_root_path', dir = [homedir '/depth_video/'];
      end
    case 'csail'
      homedir = '/afs/csail.mit.edu/u/y/yuvval';
      switch dir_type
        case 'homedir', dir = homedir;
        case 'datasets_base', dir  = ['/storage/yuvval/data/video/'];
        case 'preproc_base', dir  = ['/storage/yuvval/depth_vl/' ...
                              'preproc/'];
        case 'proj_root_path', dir = [homedir '/depth_video/'];
      end
    otherwise
      error('Unknown domain: %s', domain);
end

% TODO where is proj_root_path used? 
%    
