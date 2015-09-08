function [all_neigh_pairs_inds, all_neigh_L2_dist] = get_inds_of_all_pixels_neighbours(siz, L1rad)
% Returns all combinations of the -L1rad:L1rad (x,y) proximal pixels ids, s.t. the image is stacked as a
% column vector. Pixels outside the image boundaries are discarded.
% siz = [Nx, Ny, Nframes]
%
% if siz is 3D, then we stacked them to a single vector and repeat the process for each frame
%
% e.g. The pixels of a 2x3 image are indexed in a stacked vector as:
%    1  3  5
%    2  4  6
% 
%       7  9 11
%       8 10 12
%
% We return the following combinations:
% get_inds_of_all_pixels_neighbours([2 3 2], 1)
% 
% ans =
%      1     1
%      1     3
%      1     2
%      1     4
%      2     1
%      2     3
%      2     2
%      2     4
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
%      7     7
%      7     9
%      7     8
%      7    10
%      8     7
%      8     9
%      8     8
%      8    10
%      9     7
%      9     9
%      9    11
%      9     8
%      9    10
%      9    12
%     10     7
%     10     9
%     10    11
%     10     8
%     10    10
%     10    12
%     11     9
%     11    11
%     11    10
%     11    12
%     12     9
%     12    11
%     12    10
%     12    12

if nargin <2
    L1rad = 1;
end

Nimg = prod(siz(1:2));
arr_ids = (1:Nimg).';
[I,J] = ind2sub(siz, arr_ids);

all_neighbours_ind_combs = allcomb(arr_ids, -L1rad:L1rad, -L1rad:L1rad);
all_neighbours_sub_combs = [zeros(size(all_neighbours_ind_combs,1),1) all_neighbours_ind_combs];
all_neighbours_sub_combs(:,1:2) = [I(all_neighbours_ind_combs(:,1)) J(all_neighbours_ind_combs(:,1))];
all_neigh_L2_dist = sqrt(sum(all_neighbours_sub_combs(:,3:4).^2,2));
all_neighbours_sub_with_margins = [all_neighbours_ind_combs(:,1), (all_neighbours_sub_combs(:,1:2) + all_neighbours_sub_combs(:,3:4)), all_neigh_L2_dist];
all_neighbours_sub = all_neighbours_sub_with_margins;
all_neighbours_sub(all_neighbours_sub(:,2) < 1, :) = [];
all_neighbours_sub(all_neighbours_sub(:,3) < 1, :) = [];
all_neighbours_sub(all_neighbours_sub(:,2) > max(I), :) = [];
all_neighbours_sub(all_neighbours_sub(:,3) > max(J), :) = [];

inds = sub2ind(siz, all_neighbours_sub(:,2), all_neighbours_sub(:,3));
all_neigh_pairs_inds = [all_neighbours_sub(:,1) , inds];
all_neigh_L2_dist = all_neighbours_sub(:,4);

% repeat for all frames, by duplicating and adding a respective offset  
% (of number of pixels per image) per duplication
N_neighbours_per_frame = size(all_neigh_pairs_inds,1);
Nframes = siz(3);
all_neigh_pairs_inds = repmat(all_neigh_pairs_inds, Nframes,1);
all_neigh_L2_dist = repmat(all_neigh_L2_dist, Nframes,1);
for t=2:Nframes
all_neigh_pairs_inds(((t-1)*N_neighbours_per_frame + 1):(t*N_neighbours_per_frame),:) ...
    = (t-1)*Nimg + all_neigh_pairs_inds(((t-1)*N_neighbours_per_frame + 1):(t*N_neighbours_per_frame),:);
end
