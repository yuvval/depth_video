%nel command window settare i=:
clear
clc
close all
i=0;
%for i=42:47
    switch i      
        case 0
            seq='me';
            alg='CLG-TV';
            name=strcat(seq,'_',alg);
            %cd videos
            obj = VideoReader('../videos/2chairs_approach_behind.avi');
            video = obj.read();
            %cd ./..
            frame_sample_interval=10;
            for k=1:frame_sample_interval:size(video,4)
            	acquisitionSequence
                main2
            end
                        
    end
    clear
    close all
    
%end