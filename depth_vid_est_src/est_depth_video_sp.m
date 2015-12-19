function est_depths_sp_vid = est_depth_video_sp(ppvid, mnFrameID, mxFrameID, depth_est_params)
if mxFrameID == -1
    mxFrameID = size(ppvid.depth_frames,3);
end % if mxFrameID == -1

C_tact = depth_est_params.C_tact;
rho_tact = depth_est_params.rho_tact;
rho_OF = depth_est_params.rho_OF;
C_OF = depth_est_params.C_OF;
rho_div = depth_est_params.rho_div;
dist_around_OF = depth_est_params.dist_around_OF; % L1 radius around an OF pixel

%%

tmp = ppvid.depth_frames(:,:,1);
imsize_depth = size(tmp);
imsize_OF = size(ppvid.opflow_frames{1}(:,:,1));
Nframes = 1+mxFrameID - mnFrameID;%size(ppvid.depth_frames,3);
Nimg = prod(imsize_depth);
N = Nimg*Nframes;
%% convert depth, gradients and optical flow from cell array to a 3D array (2D sptial, time).
[superpxl_frames, depths_superpxl, n_superpxl] = deal({});
depth_frames = ppvid.depth_frames(:,:, mnFrameID:mxFrameID);
[uOFframes, vOFframes] =  deal(zeros([imsize_OF (Nframes-1)]));
t = 1;

d = [];
for n=mnFrameID:mxFrameID
    superpxl_frames{t} = ppvid.superpxl_frames{n};
    depths_superpxl{t} = ppvid.depths_superpxl{n};
    d = [d ; depths_superpxl{t}];
    n_superpxl{t} = ppvid.n_superpxl{n};
    
    if t<Nframes
        uOFframes(:,:,t) = ppvid.opflow_frames{n}(:,:,1);
        vOFframes(:,:,t) = ppvid.opflow_frames{n}(:,:,2);
    end
    t=t+1;
end
% depth_frames = depth_frames/mean(depth_frames(:));
%%
% d = double(depth_frames(:));
Npxl = length(d);

%%
[XX0, YY0] = meshgrid(1:imsize_OF(2), 1:imsize_OF(1));
N_neigh_per_pxl = (1 + 2*dist_around_OF)^2; % number of neighbours around each OF pixel
[orig_ids_map, OF_ids_map, OF_const_mag] = deal([]);
sp_pairs_abs_OF = {};
sp_pairs_overlap_prop_cell = {};
orig_ids_map_cell = {};
n = 0;
for t=1:(Nframes-1)
    XpostOF = round(XX0 + uOFframes(:,:,t));
    
    YpostOF = round(YY0 + vOFframes(:,:,t));
    
    [sp_pairs_inds, sp_pairs_overlap_prop, sp_pairs_abs_OF_curr] = ...
        get_superpxls_relative_overlaps_following_OF...
        (ppvid.superpxl_frames{t}, ppvid.n_superpxl{t}, ...
        ppvid.superpxl_frames{t+1}, ppvid.n_superpxl{t+1}, ...
        XpostOF, YpostOF);
    
    orig_ids_map = [orig_ids_map ; sp_pairs_inds(:,1) + n];
    OF_ids_map   = [OF_ids_map   ; n + double(ppvid.n_superpxl{t}) + sp_pairs_inds(:,2)]; % -1 because superpxl labels count starts from 1
    OF_const_mag = [OF_const_mag ; sp_pairs_overlap_prop];
    orig_ids_map_cell{end+1} = orig_ids_map;
    sp_pairs_overlap_prop_cell{end+1} = sp_pairs_overlap_prop;
    sp_pairs_abs_OF{end+1} = sp_pairs_abs_OF_curr;
    n = n + double(ppvid.n_superpxl{t});
    
    %     n = length(sp_pairs_inds);
    %     orig_ids_map(1:n,t) = sp_pairs_inds(:,1) + (t-1)*double(ppvid.n_superpxl{t});
    %     OF_ids_map(1:n,t) = t*double(ppvid.n_superpxl{t}) + sp_pairs_inds(:,2)-1; % -1 because superpxl labels count starts from 1
    %     OFmag(1:n,t) = sp_pairs_overlap_prop;
end
orig_ids_map = orig_ids_map(:);
OF_ids_map = OF_ids_map(:);
OF_const_mag = OF_const_mag(:);

orig_ids_map(isnan(orig_ids_map)) = [];
OF_ids_map(isnan(OF_ids_map)) = [];
OF_const_mag(isnan(OF_const_mag)) = [];

%%
OFweights = C_OF*OF_const_mag.^(1);
% OFweights = C_OF*exp((-1/rho_OF)*OFmag(:).^2);
% Winter = sparse(orig_ids_map, OF_ids_map(:), OFweights, N, N);

%% Adding a costraint that mean depth of super pixels with 50% lowest OF
% should be equal (up to a slack)

Dmean_d = spalloc(Nframes, Npxl, Npxl);

n=0;
for t=1:(Nframes)
    %%% commented out code doesn't work
%     stable_overlap_sp_ids = find(sp_pairs_overlap_prop_cell{t} > 0.8);  % get SP with >80% overlap
%     [~, SPs_IDs_abs_OF_asc]  = sort(sp_pairs_abs_OF{t});
%     stable_OF_SPs_ids = SPs_IDs_abs_OF_asc(1:floor(length(SPs_IDs_abs_OF_asc)/2));
%     stable_SPs_ids = intersect(stable_overlap_sp_ids, stable_OF_SPs_ids);
%     stable_SPs = orig_ids_map_cell{t}(stable_SPs_ids);
%     n_stable = numel(stable_SPs);
%     Dmean_d(t, n+stable_SPs) = 1/n_stable;

    
    stable_SPs = 1:ppvid.n_superpxl{t};
    Dmean_d(t, n+stable_SPs) = 1;
    
    n = n + double(ppvid.n_superpxl{t});
end


%% Adding a prior that object who touch should have an avg same
% %% depth (up to a slack that we minimize)
%
%
%
Dtact = [];
Stact = [];
Htact = [];
btact = [];

% %% IMPORTANT!! Currectly only relevant to tactile pairs!!
% % To upgrade to support more complex tactile relationship we need
% % to cluster and traverse the tactile relationships graphs with BFS/DFS
% if isfield(depth_est_params, 'detect_masks')
%     detect_masks = depth_est_params.detect_masks;
%     tactile_inds = depth_est_params.tactile_inds;
%
%     for t=1:(Nframes) % TODO: ?? should it be (Nframes - 1)??
%         [rx, ry, cx, cy, umask, stdmask]  = deal({});
%         dpframe = depth_frames(:,:,t);
%
%         Nmasks = length(detect_masks{t});
%         % tactile_neigh_inds = get_tactile_neighbours(detect_masks{t}, )
%
%         [tact1, tact2] = find(tactile_inds(:,:,t));
%         Ntact = length(tact1);
%         for m = 1:Nmasks
%             % prepare tie avg depths per mask
%             s = sum(detect_masks{t}{m},1);
%             rx{m} = nnz(s/max(s)>0.5)/2; %approx X radius of a box
%             s = sum(detect_masks{t}{m},2);
%             ry{m} = nnz(s/max(s)>0.5)/2; %approx Y radius of a box
%
%             [r,c] = find(detect_masks{t}{m}==1);
%             cx{m} = mean(c);
%             cy{m} = mean(r);
%
%             umask{m} = mean(dpframe(detect_masks{t}{m}==1));
%             stdmask{m} = std(dpframe(detect_masks{t}{m}==1));
%
%
% %             Dtact(end+1,1:Npxl) = sparse(1,Npxl);
% %             mask_inds = find(detect_masks{t}{m}(:)) + (t-1)*Nimg;
% %             Dtact(end,mask_inds) = 1;
%         end
%
%         for nt = 1:Ntact
%             % check if there is a tactile interaction between 2 detection masks
%             if tact1(nt)~=tact2(nt) % if yes
%                 ujoint = mean([umask{tact1(nt)} umask{tact2(nt)}]);
%                 % for each pixel k in those masks we add a constraint:
%                 % dk = ujoint + sk, where we add f(std{m})*sk^2 to the objective
%
%                 m = tact1(nt);
%                 stdx = stdmask{m}*rx{m}/2; % setting box boundary to 2*stddev
%                 stdy = stdmask{m}*ry{m}/2; % setting box boundary to 2*stddev
%                 msk = detect_masks{t}{m}==1;
%                 [r,c] = find(msk);
%
%                 mask_inds = find(msk(:)) + (t-1)*Nimg;
%                 Ntact = length(mask_inds);
%
% %                 Htact1 = (r.^2/stdy^2 + c.^2/stdx^2); % TODO: double check if not need to commute r-c
%                 Htact1 = exp((-1)*((r-cy{m}).^2/stdy^2 + (c-cx{m}).^2/stdx^2)/rho_tact); % TODO: double check if not need to commute r-c
%                 Dtact1 = sparse(1:Ntact, mask_inds, ones(Ntact,1), Ntact,Npxl);
%                 btact1 = ones(Ntact,1)*ujoint;
%
%
%                 m = tact2(nt);
%                 stdx = stdmask{m}*rx{m}/2; % setting box boundary to 2*stddev
%                 stdy = stdmask{m}*ry{m}/2; % setting box boundary to 2*stddev
%                 msk = detect_masks{t}{m}==1;
%                 [r,c] = find(msk);
%
%                 mask_inds = find(msk(:)) + (t-1)*Nimg;
%                 Ntact = length(mask_inds);
%
% %                 Htact2 = (r.^2/stdy^2 + c.^2/stdx^2); % TODO: double check if not need to commute r-c
%                 Htact2 = exp((-1)*((r-cy{m}).^2/stdy^2 + (c-cx{m}).^2/stdx^2)/rho_tact); % TODO: double check if not need to commute r-c
%                 Dtact2 = sparse(1:Ntact, mask_inds, ones(Ntact,1), Ntact,Npxl);
%                 btact2 = ones(Ntact,1)*ujoint;
%
%                 Htact = [Htact ; Htact1; Htact2];
%                 Dtact = [Dtact; Dtact1; Dtact2];
%                 btact= [btact; btact1; btact2];
%             end
%         end
% %         Stact_curr = sparse(Stact_curr);
% %         Stact = blkdiag(Stact,Stact_curr);
%
%     end
% end
% % Nmasks = size(Stact, 1);
Ntact = length(btact);
% Stact = speye(Ntact);
% Htact = double(Htact)*C_tact;
% btact = double(btact);

%%
% Gmagvec = Gmag_frames(:);

% all_pairs_intra = get_inds_of_all_pixels_neighbours(size(depth_frames));
% Npairs_intra = length(all_pairs_intra(:,1));
Npairs_inter = length(OFweights);
% Winter = sparse(orig_ids_map, OF_ids_map(:), OFweights, N, N);
Smean_d = sparse([sparse(Npairs_inter+Ntact,1); ones(Nframes,1); sparse(Npxl,1) ]);

% Smean_d = [];
% Dmean_d = [];
% Nframes = 0;

% Hdiag = [ones(Npxl,1) ; OFweights(:); sparse(Htact) ; zeros(Npxl,1);0];
Hdiag = [ones(Npxl,1) ; OFweights(:); sparse(Htact) ; 0; zeros(Npxl,1)];
% Hdiag = [ones(Npxl,1) ; OFweights(:); sparse(Htact) ;  zeros(Npxl,1)];
%Hdiag = [ones(Npxl,1) ; (1-Gmagvec(all_pairs_intra(:,2))); ...
%OFweights(:); zeros(Npxl,1);0;zeros(size(Sdiv,2),1)]; % with OF divergence
% Hdiag = [ones(Npxl,1) ; (1-Gmagvec(all_pairs_intra(:,2)));  zeros(Npxl,1)];
Nunknowns = length(Hdiag);
H = sparse(1:Nunknowns, 1:Nunknowns, Hdiag);
f = zeros(Nunknowns,1);

% DiDj = sparse(1:Npairs_intra, all_pairs_intra(:,1), ones(Npairs_intra, 1)) + sparse(1:Npairs_intra, all_pairs_intra(:,2), -ones(Npairs_intra, 1));
% Sij = sparse(1:Npairs_intra, 1:Npairs_intra, ones(Npairs_intra,1));

DmDn = sparse(1:Npairs_inter, orig_ids_map(:), ones(Npairs_inter, 1), Npairs_inter, Npxl) + sparse(1:Npairs_inter, OF_ids_map(:), -ones(Npairs_inter, 1), Npairs_inter, Npxl);
Smn = sparse(1:Npairs_inter, 1:Npairs_inter, ones(Npairs_inter,1));

Si = sparse(1:Npxl, 1:1:Npxl, -ones(Npxl,1));
Di = sparse(1:Npxl, 1:1:Npxl, ones(Npxl,1));

% % Dix = sparse(kron(eye(Nframes), ones(1, Nimg)));
% % Dx = sparse([sparse(Npairs_intra + Npairs_inter+Ntact+Npxl,1); ones(Nframes,1) ]);
% Dix = [];
% Dx = []
%Dx = sparse([sparse(Npairs_intra + Npairs_inter+Npxl,1); ...
%ones(Nframes,1) ; sparse(Nmasks,1)]); % with OF divergence
% Aeq = [ [sparse(Npairs_intra + Npairs_inter, Npxl);Si] [Sij;sparse(Npairs_inter+Npxl, Npairs_intra)] [sparse(Npairs_intra, Npairs_inter); Smn ;sparse(Npxl, Npairs_inter)] [ DiDj ; DmDn; Di] ] ;
Aeq = [ [sparse(Npairs_inter+Ntact, Npxl);sparse(Nframes, Npxl);Si] [Smn ;sparse(Ntact+Nframes, Npairs_inter); sparse(Npxl, Npairs_inter)] [sparse(Npairs_inter, Ntact);Stact;sparse(Npxl+Nframes,Ntact)] Smean_d [ DmDn; Dtact;Dmean_d;Di]] ;
% Aeq = [ [sparse(Npairs_intra + Npairs_inter, Npxl);Si;sparse(Nframes+Nmasks, Npxl)] ...
%     [Sij;sparse(Npairs_inter+Npxl, Npairs_intra);sparse(Nframes+Nmasks, Npairs_intra)] ...
%     [sparse(Npairs_intra, Npairs_inter); Smn ;sparse(Npxl, Npairs_inter); sparse(Nframes+Nmasks, Npairs_inter)] ...
%     [ DiDj ; DmDn; Di; Dix ; Dij_div ] Dx ...
%     [sparse(Npairs_intra + Npairs_inter + Npxl + Nframes, size(Sdiv,2)) ; Sdiv] ] ;

% beq = [sparse(Npairs_inter,1); btact; d;sparse(Nframes,1)];
beq = [sparse(Npairs_inter,1); btact; sparse(Nframes,1); d];
%beq = [sparse(Npairs_intra+Npairs_inter,1); ...
%d;sparse(Nframes+Nmasks,1)]; % with OF divergence

if strcmp(get_hostname(), 'yuval')
    % running quadprog on my desktop machine is too slow, so exiting
    return
end

%%
disp('zz');
tic
X = solveQuadWithEq(H,f, 0, Aeq, beq);
% X = quadprog(H, f, [], [], Aeq, beq);
toc

qp_res_d = X((end-Npxl+1):(end));
n=1;
for t = 1:(1+mxFrameID - mnFrameID)
    est_depths_sp_vid{t} = qp_res_d(n:(n+ppvid.n_superpxl{t}-1));
    n = n+ppvid.n_superpxl{t};
end

% est_depth_vid = reshape(full(X((end-Npxl):(end-1))), size(depth_frames,1), size(depth_frames,2), size(depth_frames,3));
%est_depth_vid =
%reshape(full(X((end-Npxl-size(Sdiv,2)):(end-1-size(Sdiv,2)))),...
%size(depth_frames,1), size(depth_frames,2), size(depth_frames,3));% with OF divergence
