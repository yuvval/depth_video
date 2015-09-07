function all_neigh_pairs_inds = get_inds_of_all_pixels_neighbours(siz)
% Returns all combinations of the -1,0,+1 (x,y) proximal pixels ids, s.t. the image is stacked as a
% column vector. Pixels outside the image boundaries are discarded.
%
% e.g. The pixels of a 2x3 image are indexed in a stacked vector as:
%    1  3  5
%    2  4  6
% 
% We return the following combinations:
% get_inds_of_all_pixels_neighbours([2 3])
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


N = prod(siz);
arr_ids = (1:N).';
[I,J] = ind2sub(siz, arr_ids);

all_neighbours_ind_combs = allcomb(arr_ids, -1:1, -1:1);
all_neighbours_sub_combs = [zeros(size(all_neighbours_ind_combs,1),1) all_neighbours_ind_combs];
all_neighbours_sub_combs(:,1:2) = [I(all_neighbours_ind_combs(:,1)) J(all_neighbours_ind_combs(:,1))];

all_neighbours_sub_with_margins = [all_neighbours_ind_combs(:,1) (all_neighbours_sub_combs(:,1:2) + all_neighbours_sub_combs(:,3:4)) ];
all_neighbours_sub = all_neighbours_sub_with_margins;
all_neighbours_sub(all_neighbours_sub(:,2) == 0, :) = [];
all_neighbours_sub(all_neighbours_sub(:,3) == 0, :) = [];
all_neighbours_sub(all_neighbours_sub(:,2) > max(I), :) = [];
all_neighbours_sub(all_neighbours_sub(:,3) > max(J), :) = [];

inds = sub2ind(siz, all_neighbours_sub(:,2), all_neighbours_sub(:,3));
all_neigh_pairs_inds = [all_neighbours_sub(:,1) , inds];


