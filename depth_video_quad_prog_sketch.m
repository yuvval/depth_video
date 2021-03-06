close all
clear
initdirs

%%
load preprocessed_videos/approaching_toward_fsmp_15_ppvid.mat

imsize = size(ppvid.frames{1}(:,:,1));
Nframes = length(ppvid.frames);
Nimg = prod(imsize);
N = Nimg*Nframes;
%% convert depth, gradients and optical flow from cell array to a 3D array (2D sptial, time).

[Gmag_frames, depth_frames] = deal(zeros([imsize Nframes]));
[uOFframes, vOFframes] =  deal(zeros([imsize (Nframes-1)]));
for t=1:Nframes
    im = double(rgb2gray(ppvid.frames{t}));
%     [ppvid.Gx{t}, ppvid.Gy{t}] = imgradientxy(im, 'IntermediateDifference');
%     Gmag_frames(:,:,t) = imgradient(ppvid.Gx{t}, ppvid.Gy{t});
    Gmag_frames(:,:,t) = double(edge(im, 'canny'));

    depth_frames(:,:,t) = ppvid.depths_pxl{t};
    if t<Nframes
        uOFframes(:,:,t) = ppvid.uvOFs{t}(:,:,1);
        vOFframes(:,:,t) = ppvid.uvOFs{t}(:,:,2);
    end
end
% depth_frames = depth_frames/mean(depth_frames(:));

%%
[XX0, YY0] = meshgrid(1:imsize(2), 1:imsize(1));
dist_around_OF = 2; % L1 radius around an OF pixel
N_neigh_per_pxl = (1 + 2*dist_around_OF)^2; % number of neighbours around each OF pixel
[orig_ids_map, OF_ids_map, OFmag] = deal(nan(Nimg*N_neigh_per_pxl, (Nframes -1)));

for t=1:(Nframes-1)
    XpostOF = round(XX0 + uOFframes(:,:,t));
    
    YpostOF = round(YY0 + vOFframes(:,:,t));

    [all_neigh_pairs_inds, all_neigh_L2_dist] = get_inds_of_all_pixels_neighbours_following_OF(imsize , XpostOF, YpostOF, dist_around_OF);
    n = length(all_neigh_L2_dist);
    orig_ids_map(1:n,t) = all_neigh_pairs_inds(:,1) + (t-1)*Nimg;
    OF_ids_map(1:n,t) = t*Nimg + all_neigh_pairs_inds(:,1);
    OFmag(1:n,t) = all_neigh_L2_dist;
end
orig_ids_map = orig_ids_map(:);
OF_ids_map = OF_ids_map(:);
OFmag = OFmag(:);

orig_ids_map(isnan(orig_ids_map)) = [];
OF_ids_map(isnan(OF_ids_map)) = [];
OFmag(isnan(OFmag)) = [];

%%
rhoOF =1e-1;
C_OF = 1/sqrt(2*pi*rhoOF);
OFweights = C_OF*exp((-1/rhoOF)*OFmag(:).^2);
Winter = sparse(orig_ids_map, OF_ids_map(:), OFweights, N, N);
%%
d = double(depth_frames(:));

Gmag2 = Gmag_frames.^2;
Gmag2vec = Gmag2(:);

all_pairs = get_inds_of_all_pixels_neighbours(size(depth_frames));

%%
rhoE = 0.01;
C_E = 0.1;
% imtmp = exp((-1/rhoE)*Gmag2);
% imshow(imtmp/max(imtmp(:)))
all_pairs_edge_weights = C_E*exp((-1/rhoE)*Gmag2vec(all_pairs(:,2)));
Wintra = sparse(all_pairs(:,1), all_pairs(:,2), all_pairs_edge_weights, N, N);
W = Winter;
Wc = sparse(1:N, 1:N, sum(W,1));
Wr = sparse(1:N, 1:N, sum(W,2));
Wq = sparse(Wr + Wc -2*W +speye(size(W)));

%%
tic
X = quadprog(Wq, -2*d);
toc

%% Visualization and post processing

%% 
qp_depth_vid = reshape(X, size(depth_frames,1), size(depth_frames,2), size(depth_frames,3));
fname_gif = sprintf('depth_video_quad_prog_rhoE_%2.1f_rhoOF_%2.3f.gif', rhoE, rhoOF);
himg = figure;
for t=1:Nframes
    imshow(qp_depth_vid(:,:,t));shg
    xlabel(sprintf('\\rhoE = %2.3f, \\rhoOF = %2.3f\nafter quad prog with edges constraints', rhoE, rhoOF),'Interpreter','Tex');    
    caxis([min(X) max(X)]);
    colormap(flipud(parula));shg
    drawnow
    pause(.5)
    save_animated_gif_frame(fname_gif, t==1);
end
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);
% pause(5)
% %% eval Z movement of a stationary area
% x1x2 = 220:310;
% y1y2 = 140:180;
% x1 = min(x1x2); x2 = max(x1x2); 
% y1 = min(y1y2); y2 = max(y1y2); 
% figure(himg)
% line([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]', 'color', 'b', 'linewidth', 3, 'linestyle', '-');
% 
% stat_area_depths = qp_depth_vid(y1y2, x1x2,:);
% stat_area_depths_diff = 0*stat_area_depths;
% stat_area_depths_diff(:,:,end) = [];
% 
% for t = 1:7
%     stat_area_depths_diff(:,:,t) = stat_area_depths(:,:,t+1) - stat_area_depths(:,:,t);
% end
% figure
% hist(stat_area_depths_diff(:));
% title('histogram of Z axis movement of a stationary area');
% shg
% 
% 
% %% eval Z movement of a moving area
% x1x2 = 136:160;
% y1y2 = 127:154;
% x1 = min(x1x2); x2 = max(x1x2); 
% y1 = min(y1y2); y2 = max(y1y2); 
% figure(himg)
% line([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]', 'color', 'y', 'linewidth', 3, 'linestyle', '-');
% 
% 
% moving_area_depths = qp_depth_vid(y1y2, x1x2,:);
% moving_area_depths_diff = 0*moving_area_depths;
% moving_area_depths_diff(:,:,end) = [];
% 
% for t = 1:7
%     moving_area_depths_diff(:,:,t) = moving_area_depths(:,:,t+1) - moving_area_depths(:,:,t);
% end
% figure
% hist(moving_area_depths_diff(:));
% title('histogram of Z axis movement of a moving area');
% shg

% figure
% hist(stat_area_depths_diff(:));



%% Save RAW depth video to gif
Gmag_mask_frames = (1-Gmag_frames);
mn_depth = min(depth_frames (:));
depth_vid_masked = depth_frames.*Gmag_mask_frames;
depth_vid_masked(depth_vid_masked==0) = depth_vid_masked(depth_vid_masked==0) + mn_depth;

fname_gif = sprintf('depth_video_raw.gif');
himg = figure;
for t=1:Nframes
    imshow(depth_vid_masked(:,:,t));shg
    xlabel(sprintf('depth video, raw depth estimation'),'Interpreter','Tex');    
    caxis([min(depth_frames(:)) max(depth_frames(:))]);
    colormap(flipud(parula));shg
    drawnow
    save_animated_gif_frame(fname_gif, t==1);
end
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);
