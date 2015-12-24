function [all_sp_pairs_intra, weighs_sp_intra] = ...
    get_super_pixels_intra_weighs(img, spLblImg, all_neigh_pairs_inds)

if nargin < 3 
    % preferably evaluate this once externally and reuse in each call.
    % than use default assignment
    [all_neigh_pairs_inds, ~] = ...
        get_inds_of_all_pixels_neighbours(size(spLblImg), 1);

end
spLblImg = double(spLblImg);

img_edges = double(imgradient(rgb2gray(img), 'sobel'));
%normalizing edge to a unit variance around its mean 
img_edges = img_edges - mean(img_edges(:));
img_edges = img_edges/std(img_edges(:));
img_edges = img_edges + min(img_edges(:));


% find the super pixels labels for all pixels pairs in distance 1
sp_pairs = [ spLblImg(all_neigh_pairs_inds(:,1)) ...
             spLblImg(all_neigh_pairs_inds(:,2)) ];
% trim the cases that two pixels belong to the same super pixel

same_lbl_ids = sp_pairs(:,1) == sp_pairs(:,2);
edges_neigh_pairs_inds = all_neigh_pairs_inds;
edges_neigh_pairs_inds(same_lbl_ids, :) = [];
sp_pairs(same_lbl_ids, :) = [];
unique_edge_pixels  = unique(edges_neigh_pairs_inds(:,1));
edge_pixels_neighbours = ...
    all_neigh_pairs_inds(ismember(all_neigh_pairs_inds(:,1), ...
                                  unique_edge_pixels),:);

% So now we know for each pixel which super pixels edge is relates
% to and its radius 1 neighbours

% We take the edge detector weighs per super pixel edge pixel and
% its neighbours
edge_detector_weighs_per_sp_edge_pixel = ...
    img_edges(unique_edge_pixels);
edge_detector_weighs_per_sp_edge_pixel_neigh = ...
    img_edges(edge_pixels_neighbours(:,2));


ew = edge_detector_weighs_per_sp_edge_pixel;
integ_ew_neigh = ...
    cumsum(edge_detector_weighs_per_sp_edge_pixel_neigh);
sample_integ_ids = (diff([edge_pixels_neighbours(:,1);length(edge_pixels_neighbours)]) ~= 0);
smp_integ_ew_neigh = integ_ew_neigh(sample_integ_ids);
ew_neigh_weighs = diff([0; smp_integ_ew_neigh]);
u_sp_pxls = edge_pixels_neighbours(sample_integ_ids, 1);
[~,u_sp_pxls_permute] = sort(u_sp_pxls); 
edge_detector_weighs_per_sp_edge_pixel = ew+ew_neigh_weighs(u_sp_pxls_permute);

% Now we have the weigh per pixel on a super pixel edge, next we
% will sum all the weighs per edge

u_pairs_weigh_mat = zeros(max(sp_pairs(:)), max(sp_pairs(:)));

E = zeros(size(img_edges));
E(unique_edge_pixels) = edge_detector_weighs_per_sp_edge_pixel;

sp_pairs_ind = sub2ind(size(u_pairs_weigh_mat), sp_pairs(:,1), sp_pairs(:,2));
sumE =  E(edges_neigh_pairs_inds(:,1)) + E(edges_neigh_pairs_inds(:,2));
weighs_vec = accumarray(sp_pairs_ind, sumE(:));
u_pairs_weigh_mat(1:numel(weighs_vec)) = weighs_vec(:);

% represent as vectors of indices and values
[i,j,weighs_sp_intra] = find(u_pairs_weigh_mat);
all_sp_pairs_intra = [i j];
