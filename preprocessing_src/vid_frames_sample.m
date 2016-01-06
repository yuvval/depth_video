function samples_vec = vid_frames_sample(video, prepr_params)
% function samples_vec = vid_frames_sample(video, prepr_params)

if prepr_params.mxFrameID == -1
    if iscell(video)
        prepr_params.mxFrameID = length(video);
    else 
        prepr_params.mxFrameID = size(video,4);
    end
end

samples_vec = ...
    prepr_params.mnFrameID:prepr_params.sample_interval: ...
    prepr_params.mxFrameID;


