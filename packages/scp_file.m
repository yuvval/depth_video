function scp_file(fullfname, remotepath)
% function scp_file(fullfname, remotepath)
% default remote path is www/figs on cortex
if nargin < 2
    % if not on cortex then scp to cortex by default
    if isempty(regexp(get_hostname(), 'ctx', 'match'))
        remotepath = 'yuvval@ctx01.lnx.biu.ac.il:~/www/figs/';
    else % otherwise, return (do nothing)
        return
    end

end

[~, fname, ext] = fileparts(fullfname);

system(['scp -rp ' fullfname ' ' remotepath fname ext]);
