function scp_file(fname, remotepath)
% default remote path is www/figs on cortex
if nargin < 2
    % if not on cortex then scp to cortex by default
    if isempty(regexp(get_hostname(), 'ctx', 'match'))
        remotepath = 'yuvval@ctx01.lnx.biu.ac.il:';
    else % otherwise, return (do nothing)
        return
    end

end

system(['scp ' fname ' ' remotepath fname]);
