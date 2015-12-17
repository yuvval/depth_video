function [rgb_frames, depth_frames, ...
    scaled_rgb_frames, scaled_depth_frames, camera_info] = ...
    get_frames_mat(video_name, sample_interval, scale_to_resolution, base_path)

initdirs

if nargin <2
    sample_interval = 1;
end

if nargin <3
    scale_to_resolution = []; % Ys, Xs
end

if nargin <4
    base_path = ['../videos/'];
end

load([base_path  video_name '.mat']);

camera_info = frames;

max_frames = length(frames);
frames_range = 1:sample_interval:max_frames;
Nframes = length(frames_range);

cnt = 1;
for frameId = frames_range
    im = frames(:,:,:,frameId);
    if cnt == 1
        imsize = size(im);         
        rgb_frames = nan(imsize(1), imsize(2), 3, Nframes);
        depth_frames = nan(imsize(1), imsize(2), Nframes);    
        
        if ~isempty(scale_to_resolution)
            scaled_rgb_frames = nan(scale_to_resolution(1), scale_to_resolution(2), 3, Nframes);
            scaled_depth_frames = nan(scale_to_resolution(1), scale_to_resolution(2), Nframes);    
        else
            scaled_rgb_frames = [];
            scaled_depth_frames = [];
        end            
    end
    rgb_frames(:,:,:, cnt) = im;
    
    depth_frames(:,:, cnt) = zeros(size(im,1), size(im, 2));
    
    if ~isempty(scale_to_resolution)
        scaled_rgb_frames(:,:,:, cnt)  = imresize(rgb_frames(:,:, :, cnt), scale_to_resolution);
        scaled_depth_frames(:,:, cnt)  = imresize(depth_frames(:,:, cnt), scale_to_resolution);
    end
    
    cnt = cnt +1;

end
%% show the 2D image
%     subplot(1,5,1:2); imshow(rgb);
%     dmn = min(depth{1}(:));
%     dmx = max(depth{1}(:));
%     dmx = (dmx-dmn)/2 + dmn;
%     davg = mean(depth{1}(:));
%     dstd = std(depth{1}(:));
%     dmn = max(davg-dstd, 0);
%     dmx = davg+dstd;
% %     subplot(1,2,2); imshow(log2(depth)); caxis([8 13]); colormap(flipud(bone));
%     subplot(1,5,3:4); imshow(depth{frameId}); caxis([dmn dmx]); colormap(cmap);
%     title('Kinect depth');
% %     subplot(1,2,2); imshow(depth); caxis('auto'); colormap(flipud(bone));
% % colorbar
%     subplot(1,5,5); 
%     imshow(ones(20,1));
%     caxis([dmn dmx]/1000);  colormap(cmap); h = colorbar;
%     ylabel(h, 'meters');
% pause(0.01)
    
%% 3D point for the frame

% K = frames.K; % K is [fx 0 cx; 0 fy cy; 0 0 1];
% cx = K(1,3);cy = K(2,3);
% fx = K(1,1);fy = K(2,2);

%     depthInpaint = depth/1000;
%     [x,y] = meshgrid(1:640, 1:480); 
%     Xworld = (x-cx).*depthInpaint*1/fx;
%     Yworld = (y-cy).*depthInpaint*1/fy;
%     Zworld = depthInpaint;
%     validM = depth~=0;
%     XYZworldframe = [Xworld(:)'; Yworld(:)'; Zworld(:)'];
%     valid = validM(:)';   
%     
%     % XYZworldframe 3xn and RGB 3xn
%     RGB = [reshape(rgb(:,:,1),1,[]);reshape(rgb(:,:,2),1,[]);reshape(rgb(:,:,3),1,[])];
%     XYZpoints = XYZworldframe(:,valid);
%     RGBpoints = RGB(:,valid);
%     
%     % display in 3D: subsample to avoid too much to display.
%     XYZpoints = XYZpoints(:,1:20:end);
%     RGBpoints = RGBpoints(:,1:20:end);
%     figure, scatter3(XYZpoints(1,:),XYZpoints(2,:),XYZpoints(3,:),ones(1,size(XYZpoints,2)),double(RGBpoints)'/255,'filled');
%     axis equal; view(0,-90);
%     pause;
