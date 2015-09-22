function [all_neigh_pairs_inds, all_neigh_L2_dist, all_pairs_OF_mag] = get_inds_of_all_pixels_neighbours_following_OF(siz, XpostOF, YpostOF, L1rad)
% Returns all combinations pixel ids of the -L1rad:L1rad (x,y) neightbouring pixels following optical flow
% ids, s.t. the image is stacked as a column vector. Pixels outside the image boundaries are discarded.
% siz = [Nx, Ny, Nframes]
%
% if siz is 3D, then ... (TBD)
%
% e.g. The pixels of a 2x3 image are indexed in a stacked vector as:
%    1  3  5
%    2  4  6
% 
% The given cooridunates post optical flow are (note that pixel #2 has moved out of the image boundaries):
% Xpost =
% 
%      2     2     3
%      0     2     3
% 
% Ypost =
% 
%      1     1     1
%      2     2     2
%
% We return the following combinations (notice how pixels #1 & #2 map):
% get_inds_of_all_pixels_neighbours_following_OF([2 3 1], Xpost, Ypost, 1)
% 
% ans =
%      1     1
%      1     3
%      1     5
%      1     2
%      1     4
%      1     6
%      2     1
%      2     2
%      3     1
%      3     3
%      3     5
%      3     2
%      3     4
%      3     6
%      4     1
%      4     3
%      4     5
%      4     2
%      4     4
%      4     6
%      5     3
%      5     5
%      5     4
%      5     6
%      6     3
%      6     5
%      6     4
%      6     6

if nargin <2
    L1rad = 1;
end

Nimg = prod(siz(1:2));
arr_ids = (1:Nimg).';

I = YpostOF(:);
J = XpostOF(:);

all_neighbours_ind_combs = allcomb(arr_ids, -L1rad:L1rad, -L1rad:L1rad);


all_neighbours_sub_combs = [zeros(size(all_neighbours_ind_combs,1),1) all_neighbours_ind_combs];
all_neighbours_sub_combs(:,1:2) = [I(all_neighbours_ind_combs(:,1)) J(all_neighbours_ind_combs(:,1))];
all_neigh_L2_dist = sqrt(sum(all_neighbours_sub_combs(:,3:4).^2,2));
all_neighbours_sub_with_margins = [all_neighbours_ind_combs(:,1), (all_neighbours_sub_combs(:,1:2) + all_neighbours_sub_combs(:,3:4)), all_neigh_L2_dist];
all_neighbours_sub = all_neighbours_sub_with_margins;
all_neighbours_sub(all_neighbours_sub(:,2) < 1, :) = [];
all_neighbours_sub(all_neighbours_sub(:,3) < 1, :) = [];
all_neighbours_sub(all_neighbours_sub(:,2) > siz(1), :) = [];
all_neighbours_sub(all_neighbours_sub(:,3) > siz(2), :) = [];

inds = sub2ind(siz, all_neighbours_sub(:,2), all_neighbours_sub(:,3));
all_neigh_pairs_inds = [all_neighbours_sub(:,1) , inds];
all_neigh_L2_dist = all_neighbours_sub(:,4);

% % repeat for all frames, by duplicating and adding a respective offset  
% % (of number of pixels per image) per duplication
% N_neighbours_per_frame = size(all_neigh_pairs_inds,1);
% Nframes = siz(3);
% all_neigh_pairs_inds = repmat(all_neigh_pairs_inds, Nframes,1);
% all_neigh_L2_dist = repmat(all_neigh_L2_dist, Nframes,1);
% for t=2:Nframes
% all_neigh_pairs_inds(((t-1)*N_neighbours_per_frame + 1):(t*N_neighbours_per_frame),:) ...
%     = (t-1)*Nimg + all_neigh_pairs_inds(((t-1)*N_neighbours_per_frame + 1):(t*N_neighbours_per_frame),:);
% end
