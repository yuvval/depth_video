function start_parallel_pool(max_workers, parpool_base_path)
% function start_parallel_pool(max_workers, temp_path)
% start a parpool
% max_workers = max number of workers. if 0, then we don't create the pool
% temp_path for the parpool

% check if parpool is available
isOpen = ~isempty(gcp('nocreate')); % If no pool, do not create new one.
if isOpen
    delete(gcp) % delete current pool if available
end

if max_workers >0
    pc = parcluster;
    pc.NumWorkers = max_workers;
    
    parpool_path=tempname(parpool_base_path);
    mkdir(parpool_path);
    pc.JobStorageLocation=parpool_path;
    fprintf('parpool path = %s\n', parpool_path)
    parpool(pc) %init a new pool
end

end
