function [feat_mean_per_mask, feat_std_per_mask, num_pxl_per_mask, feat_metric_per_mask] = get_featuremap_statistics_per_mask(masks, featuremap)
[feat_mean_per_mask, feat_std_per_mask, num_pxl_per_mask, feat_metric_per_mask] = deal(nan(1, length(masks)));
for c = 1:length(masks)
    feat_mean_per_mask(c) =  mean(featuremap(masks{c}));
    feat_std_per_mask(c) =  std(featuremap(masks{c}));
    feat_metric_per_mask(c) =  feat_mean_per_mask(c).^2 + (feat_std_per_mask(c)).^2;
    num_pxl_per_mask(c) = sum(masks{c}(:));
end
