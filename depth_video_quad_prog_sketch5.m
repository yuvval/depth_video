% function qp_depth_vid = depth_video_quad_prog_sketch3(rhoOF, C_OF)

rhoOF =1;
% C_OF = 1/sqrt(2*pi*rhoOF);
C_OF = 1;
rho_div = 1;
dist_around_OF = 1; % L1 radius around an OF pixel

close all
proj_root_path = './';
initdirs

%%
load preprocessed_videos/approaching_toward_fsmp_15_ppvid.mat

imsize = size(ppvid.frames{1}(:,:,1));
Nframes = length(ppvid.frames);
Nimg = prod(imsize);
N = Nimg*Nframes;
%% convert depth, gradients and optical flow from cell array to a 3D array (2D sptial, time).

[Gmag_frames, depth_frames, OFdiv_frames] = deal(zeros([imsize Nframes]));
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
        OFdiv_frames(:,:,t) = divergence(uOFframes(:,:,t), vOFframes(:,:,t));
    end
end
% depth_frames = depth_frames/mean(depth_frames(:));
%%
d = double(depth_frames(:));
Npxl = length(d);

%%
[XX0, YY0] = meshgrid(1:imsize(2), 1:imsize(1));
N_neigh_per_pxl = (1 + 2*dist_around_OF)^2; % number of neighbours around each OF pixel
[orig_ids_map, OF_ids_map, OFmag] = deal(nan(Nimg*N_neigh_per_pxl, (Nframes -1)));

for t=1:(Nframes-1)
    XpostOF = round(XX0 + uOFframes(:,:,t));
    
    YpostOF = round(YY0 + vOFframes(:,:,t));

    [all_neigh_pairs_inds, all_neigh_L2_dist] = get_inds_of_all_pixels_neighbours_following_OF(imsize , XpostOF, YpostOF, dist_around_OF);
    n = length(all_neigh_L2_dist);
    orig_ids_map(1:n,t) = all_neigh_pairs_inds(:,1) + (t-1)*Nimg;
    OF_ids_map(1:n,t) = t*Nimg + all_neigh_pairs_inds(:,2);
    OFmag(1:n,t) = all_neigh_L2_dist;
end
orig_ids_map = orig_ids_map(:);
OF_ids_map = OF_ids_map(:);
OFmag = OFmag(:);

orig_ids_map(isnan(orig_ids_map)) = [];
OF_ids_map(isnan(OF_ids_map)) = [];
OFmag(isnan(OFmag)) = [];

%%
OFweights = C_OF*exp((-1/rhoOF)*OFmag(:).^2);
% Winter = sparse(orig_ids_map, OF_ids_map(:), OFweights, N, N);

%% evaluate the divergence of superpixels / candidates masks

%% find pixel ids of immobile patches

[masks, num_pxl, mean_OFdiv, std_OFdiv] = deal(cell(Nframes-1, 1));
Nmasks = 0;
 
for t=1:(Nframes-1)
    
    % get objects proposals per frame
    [masks{t}, num_pxl{t}, ~] = ...
        get_obj_det_candidates(ppvid.frames{t}, 'fayao_orig_superpix', ppvid.sp_info{t});
    Nmasks = Nmasks+length(masks{t});

    % get OF divegrence statistics per mask (per frame)
    [mean_OFdiv{t}, std_OFdiv{t}, ~] = ...
        get_featuremap_statistics_per_mask(masks{t}, OFdiv_frames(:,:,t));

end

%%
mean_OFdiv_vec = nan(Nmasks, 1);
Dij_div = sparse(Nmasks, Npxl);
cnt =1;
for t=1:(Nframes-1)
    curr_frm_Nmasks = length(masks{t});
    XpostOF = round(XX0 + uOFframes(:,:,t));
    YpostOF = round(YY0 + vOFframes(:,:,t));

    [all_neigh_r0_pairs_inds, ~] = get_inds_of_all_pixels_neighbours_following_OF(imsize , XpostOF, YpostOF, 0);

    for c = 1:curr_frm_Nmasks
        mean_OFdiv_currmask = mean_OFdiv{t}(c);
        orig_ids_map_curr_mask = ismember(orig_ids_map(:,1), find(masks{t}{c}));
        OFdiv_weights_curr_mask = exp(-(mean_OFdiv_currmask)^2/rho_div^2);
        OFweights(orig_ids_map_curr_mask) = OFdiv_weights_curr_mask*OFweights(orig_ids_map_curr_mask);
%         Dij_div(c, (t-1)*Nimg + find(masks{t}{c})) = -1;
%         orig_ids_curr_mask = ismember(all_neigh_r0_pairs_inds(:,1), find(masks{t}{c}));
%         OF_ids_curr_mask = all_neigh_r0_pairs_inds(orig_ids_curr_mask, 2) + t*Nimg;
%         Dij_div(c, OF_ids_curr_mask) = 1;
    end
    num_pxl_per_mask_curr_frame = arrayfun(@(k) nnz(masks{t}{k}), 1:curr_frm_Nmasks);
    mean_OFdiv_vec(cnt:(cnt+curr_frm_Nmasks-1)) = mean_OFdiv{t}.*num_pxl_per_mask_curr_frame;
    cnt = cnt + curr_frm_Nmasks;
end
Sdiv=[]; Nmasks=0;
%Sdiv = sparse([mean_OFdiv_vec -eye(Nmasks)/400000]);
% Sdiv = sparse([mean_OFdiv_vec]);
%%
% Gmag2 = Gmag_frames.^2;
% Gmag2vec = Gmag2(:);

% all_pairs_intra = get_inds_of_all_pixels_neighbours(size(depth_frames));

%%
% rhoE = 0.0001;
% C_E = 1;
% imtmp = exp((-1/rhoE)*Gmag2);
% imshow(imtmp/max(imtmp(:)))
% all_pairs_edge_weights = C_E*exp((-1/rhoE)*Gmag2vec(all_pairs_intra(:,2)));
% Wintra = sparse(all_pairs_intra(:,1), all_pairs_intra(:,2), all_pairs_edge_weights, N, N);
%%

% OFweights(OFweights<(mean(OFweights)*1e-3)) = 0;


Gmagvec = Gmag_frames(:);

all_pairs_intra = get_inds_of_all_pixels_neighbours(size(depth_frames));
Npairs_intra = length(all_pairs_intra(:,1));
Npairs_inter = length(OFweights);
% Winter = sparse(orig_ids_map, OF_ids_map(:), OFweights, N, N);

Npxl = length(d);
Hdiag = [ones(Npxl,1) ; (1-Gmagvec(all_pairs_intra(:,2))); OFweights(:); zeros(Npxl,1);0;zeros(size(Sdiv,2),1)];
% Hdiag = [ones(Npxl,1) ; (1-Gmagvec(all_pairs_intra(:,2)));  zeros(Npxl,1)];
Nunknowns = length(Hdiag);
H = sparse(1:Nunknowns, 1:Nunknowns, Hdiag);
f = zeros(Nunknowns,1);

DiDj = sparse(1:Npairs_intra, all_pairs_intra(:,1), ones(Npairs_intra, 1)) + sparse(1:Npairs_intra, all_pairs_intra(:,2), -ones(Npairs_intra, 1));
Sij = sparse(1:Npairs_intra, 1:Npairs_intra, ones(Npairs_intra,1));

DmDn = sparse(1:Npairs_inter, orig_ids_map(:), ones(Npairs_inter, 1), Npairs_inter, Npxl) + sparse(1:Npairs_inter, OF_ids_map(:), -ones(Npairs_inter, 1), Npairs_inter, Npxl);
Smn = sparse(1:Npairs_inter, 1:Npairs_inter, ones(Npairs_inter,1));

Si = sparse(1:Npxl, 1:1:Npxl, -ones(Npxl,1));
Di = sparse(1:Npxl, 1:1:Npxl, ones(Npxl,1));

Dix = sparse(kron(eye(Nframes), ones(1, Nimg)));
Dx = sparse([sparse(Npairs_intra + Npairs_inter+Npxl,1); ones(Nframes,1) ; sparse(Nmasks,1)]);
% Aeq = [ [sparse(Npairs_intra + Npairs_inter, Npxl);Si] [Sij;sparse(Npairs_inter+Npxl, Npairs_intra)] [sparse(Npairs_intra, Npairs_inter); Smn ;sparse(Npxl, Npairs_inter)] [ DiDj ; DmDn; Di] ] ;
Aeq = [ [sparse(Npairs_intra + Npairs_inter, Npxl);Si;sparse(Nframes, Npxl)] [Sij;sparse(Npairs_inter+Npxl, Npairs_intra);sparse(Nframes, Npairs_intra)] [sparse(Npairs_intra, Npairs_inter); Smn ;sparse(Npxl, Npairs_inter); sparse(Nframes, Npairs_inter)] [ DiDj ; DmDn; Di; Dix] Dx] ;
% Aeq = [ [sparse(Npairs_intra + Npairs_inter, Npxl);Si;sparse(Nframes+Nmasks, Npxl)] ...
%     [Sij;sparse(Npairs_inter+Npxl, Npairs_intra);sparse(Nframes+Nmasks, Npairs_intra)] ...
%     [sparse(Npairs_intra, Npairs_inter); Smn ;sparse(Npxl, Npairs_inter); sparse(Nframes+Nmasks, Npairs_inter)] ...
%     [ DiDj ; DmDn; Di; Dix ; Dij_div ] Dx ...
%     [sparse(Npairs_intra + Npairs_inter + Npxl + Nframes, size(Sdiv,2)) ; Sdiv] ] ;

beq = [sparse(Npairs_intra+Npairs_inter,1); d;sparse(Nframes+Nmasks,1)];


%%
disp('zz');
tic
X = solveQuadWithEq(H,f, 0, Aeq, beq);
% X = quadprog(H, f, [], [], Aeq, beq);
toc

%% Visualization and post processing
%% 
% qp_depth_vid = reshape(X, size(depth_frames,1), size(depth_frames,2), size(depth_frames,3));
X = full(X);
qp_depth_vid = reshape(X((end-Npxl-size(Sdiv,2)):(end-1-size(Sdiv,2))), size(depth_frames,1), size(depth_frames,2), size(depth_frames,3));
fname_gif = sprintf('depth_video_quad_prog_method5_rhoOF_%2.3f_COF_%2.3f_rhodiv=%2.3f_rOF_%d.gif', rhoOF, C_OF, rho_div, dist_around_OF);
himg = figure;
Gmag_mask_frames = (1-Gmag_frames);
mn_depth = min(qp_depth_vid (:));
depth_vid_masked = qp_depth_vid.*Gmag_mask_frames;
depth_vid_masked(depth_vid_masked==0) = depth_vid_masked(depth_vid_masked==0) + mn_depth;

for t=1:Nframes
    imshow(depth_vid_masked(:,:,t));shg
    xlabel(sprintf('\\rho_{OF} = %2.3f, C_{OF} = %2.3f \\rho_{div} = %2.3f, radius around OF = %d\nafter quad prog method5', rhoOF, C_OF, rho_div, dist_around_OF),'Interpreter','Tex');    
    caxis([min(qp_depth_vid (:)) max(qp_depth_vid (:))]);
    colormap(flipud(parula));shg
    drawnow
    pause(1.5)
    save_animated_gif_frame(fname_gif, t==1);
end
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);
% return
% pause(5)
%% eval Z movement of a stationary area
x1x2 = 220:310;
y1y2 = 140:180;
x1 = min(x1x2); x2 = max(x1x2); 
y1 = min(y1y2); y2 = max(y1y2); 
figure(himg)
line([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]', 'color', 'b', 'linewidth', 3, 'linestyle', '-');

stat_area_depths = qp_depth_vid(y1y2, x1x2,:);
% stat_area_depths_diff = 0*stat_area_depths;
% stat_area_depths_diff(:,:,end) = [];
% 
% for t = 1:7
%     c(:,:,t) = stat_area_depths(:,:,t+1) - stat_area_depths(:,:,t);
% end
% figure
% hist(stat_area_depths_diff(:));
% title('histogram of Z axis movement of a stationary area');
% shg


%% eval Z movement of a moving area
x1x2 = 136:160;
y1y2 = 127:154;
x1 = min(x1x2); x2 = max(x1x2); 
y1 = min(y1y2); y2 = max(y1y2); 
figure(himg)
line([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]', 'color', 'y', 'linewidth', 3, 'linestyle', '-');


moving_area_depths = qp_depth_vid(y1y2, x1x2,:);
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

%% 3D side view - inspired by http://stackoverflow.com/questions/22301571/plot-depth-image-in-matlab


figure;
for t = 1:size(depth_frames,3);

depth = (qp_depth_vid(:,:,t));
img = double((ppvid.frames{t}))/256;
[y x] = ndgrid( 1:imsize(1), 1:imsize(2) );
sel = depth > 0 ; % which points to plot
% "flatten" the matrices for scatter plot
x = x(:);
y = y(:);
img = reshape( img, [], size(img,3) );
depth = -depth(:);
% depth(:) = 1;
scatter3( depth(sel), y(sel), x(sel), 20, img( sel, : ), 'filled' );
colormap(gray);
xlim([-4.5 -2.5]);
fname_gif = sprintf('side_view_method5_rhoOF_%2.3f_COF_%2.3f_rhodiv=%2.3f_rOF_%d.gif', rhoOF, C_OF, rho_div, dist_around_OF);
xlabel(sprintf('depth\n\\rho_{OF} = %2.3f, C_{OF} = %2.3f \\rho_{div} = %2.3f, radius around OF = %d\nside view method5', rhoOF, C_OF, rho_div, dist_around_OF),'Interpreter','Tex');
ylabel('y');
view(0, -90);shg
drawnow
pause(.5);
save_animated_gif_frame(fname_gif, t==1);
end
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);
