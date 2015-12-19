function img = superpxl_intensities_to_img(superpxls_frame, superpxls_I, n_superpxls)
img = nan(size(superpxls_frame));
for k=1:n_superpxls
    img(superpxls_frame==k) = superpxls_I(k);
end
