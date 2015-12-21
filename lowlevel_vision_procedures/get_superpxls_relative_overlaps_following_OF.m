function [sp_pairs_inds, sp_pairs_overlap_prop, sp_pairs_abs_OF, sp_pairs_OF] = get_superpxls_relative_overlaps_following_OF(currSP, currSPn, nextSP, nextSPn, XpostOF, YpostOF)





srcSPvec = currSP(:);
siz = size(currSP);
% tgtSPvec = nextSP(:);


I = YpostOF(:);
J = XpostOF(:);

[XX0, YY0] = meshgrid(1:size(currSP,2), 1:size(currSP,1));
XX0 = XX0(:);
YY0 = YY0(:);
OFdiff_vec = [J-XX0, I-YY0];
OFdiff_abs = sqrt(sum(OFdiff_vec.^2,2));


srcSPvec(I<1) = [];
J (I<1) = []; 
OFdiff_vec (I<1, :) = []; 
I (I<1) = []; 

srcSPvec(I>siz(1)) = [];
J(I>siz(1)) = [];
OFdiff_vec(I>siz(1), :) = [];
I(I>siz(1)) = [];

srcSPvec(J<1) = [];
I (J<1) = []; 
OFdiff_vec (J<1, :) = []; 
J (J<1) = []; 
srcSPvec(J>siz(2)) = [];
I(J>siz(2)) = [];
OFdiff_vec(J>siz(2), :) = [];
J(J>siz(2)) = [];

% srcProjectOnTgt = full(sparse(I, J, double(srcSPvec)));
srcProjectOnTgt = nan(size(currSP));
for k = 1:length(I)
    srcProjectOnTgt(I(k), J(k)) = double(srcSPvec(k));
end

sp_src_ids = [];
sp_tgt_ids = [];
sp_pairs_overlap_prop = [];
sp_pairs_OF = [];
sp_pairs_abs_OF = [];
for k=double(1:currSPn)
    k_tgt_SPids = double(1:nextSPn).';
    kSPOnTgt_vec = double(nextSP(srcProjectOnTgt==k));
    OF_diffs_vec = OFdiff_abs(srcProjectOnTgt==k);
    if ~isempty(kSPOnTgt_vec)

        k_overlap_prop = hist(kSPOnTgt_vec, k_tgt_SPids);
        k_overlap_prop = k_overlap_prop/sum(k_overlap_prop);
        
        % taking only non zeros
        k_tgt_SPids(k_overlap_prop==0) = [];
        k_overlap_prop = k_overlap_prop(k_tgt_SPids);

        % evaluate mean optical flow from source SP to each of target SP
        mean_OF_to_tgt_SPs = zeros(numel(k_tgt_SPids),2);
        mean_abs_OF_to_tgt_SPs = 0*k_tgt_SPids;
        c = 1;
        for ktgt = k_tgt_SPids.'
            mean_OF_to_tgt_SPs(c, :) = mean(OFdiff_vec(kSPOnTgt_vec == ktgt, :), 1);
            mean_abs_OF_to_tgt_SPs(c) = norm(mean_OF_to_tgt_SPs(c,:));       
            c=c +1;
        end

        %%%%% visualization - for debugging
%         im = zeros(size(currSP));
%         im(currSP == k) = 1;
%         [xsrc, ysrc] = find(currSP == k);
%         xsrc = mean(xsrc);
%         ysrc = mean(ysrc);
%         
%         
%         c = 1;
%         for kd = k_tgt_SPids.'
%             im(nextSP == kd) = k_overlap_prop(c);
%             [xtgt, ytgt] = find(nextSP == kd);
%             xtgt = mean(xtgt);
%             ytgt = mean(ytgt);
%             if norm([xtgt ytgt] - [xsrc ysrc]) >40
%                 keyboard
%             end
% 
%             c=c+1;
%         end
%         
%         imshow(imresize(im, [240 320]*2));shg
% %         pause(0.2);
        
        % concat with previous results
        sp_src_ids = [sp_src_ids ; k*ones(numel(k_tgt_SPids),1)];
        sp_tgt_ids = [sp_tgt_ids ; k_tgt_SPids];
        sp_pairs_overlap_prop = [sp_pairs_overlap_prop ; k_overlap_prop.'];
        sp_pairs_OF = [sp_pairs_OF;mean_OF_to_tgt_SPs];
        sp_pairs_abs_OF = [sp_pairs_abs_OF;mean_abs_OF_to_tgt_SPs];
    end
end
assert(min(sp_src_ids) >=1);
assert(min(sp_tgt_ids) >=1);
assert(max(sp_src_ids) <= currSPn);
assert(max(sp_tgt_ids) <= nextSPn);

sp_pairs_inds = [sp_src_ids, sp_tgt_ids];
