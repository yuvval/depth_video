close all
clear

initdirs;

%% This solver is for a single frame with a QP with Objective & constraints

%%
load preprocessed_videos/approaching_toward_fsmp_15_ppvid.mat

t=4;

imrgb = ppvid.frames{t};
im = double(rgb2gray(ppvid.frames{t}));
Nimg = numel(im);

candidates = im2mcg(imrgb, 'fast');
[~, cand_sort_ids] = sort(candidates.scores, 'descend');
%%
masks_heat_map = zeros(size(im));
N_cand = size(candidates.scores,1);
mask_img = cell(N_cand,1);
heatmap_history = {};
k = 1;
for id = cand_sort_ids.'
mask_img{id} = get_mcg_mask_from_candidates(candidates, id);
masks_heat_map = masks_heat_map + mask_img{id};
if mod(k,100) == 0
    heatmap_history{end+1} = masks_heat_map;
end
k=k+1;
end
imshow(masks_heat_map/N_cand);
caxis([0 max(masks_heat_map(:)/N_cand)]);
colorbar; shg;
%%
figure;
for k=1:length(heatmap_history)
    imshow(heatmap_history{k}/N_cand);
    caxis([0 max(masks_heat_map(:)/N_cand)]);
    colorbar; shg;
    pause(.5)
end
return
%%
metric_OFdiv = [];
for c = 1:N_cand
[bbox, mask_pxl_ids, mask_pxl_coords] = mask2bbox(mask_img{c});
uOF = ppvid.uvOFs{t}(bbox(1):bbox(2), bbox(3):bbox(4),1);
vOF = ppvid.uvOFs{t}(bbox(1):bbox(2), bbox(3):bbox(4),2);



OFdiv = divergence(uOF, vOF);
% imshow(abs(OFdiv));shg
mask_crop = (mask_img{c}(bbox(1):bbox(2), bbox(3):bbox(4)));
zz = (OFdiv).*mask_crop;
% zz = abs(OFdiv);
% figure(1);
% imshow(mask_img{c}); 
% figure(2);
% imshow(zz); caxis([-1 1]); shg
norm_coef = sqrt(sum(mask_crop(:)));
norm_coef = 1;
mean_OFdiv(c) =  mean(OFdiv(mask_crop))/norm_coef;
std_OFdiv(c) =  std(OFdiv(mask_crop))/norm_coef;
metric_OFdiv(c) =  sqrt(mean_OFdiv(c).^2 + (std_OFdiv(c)).^2)/norm_coef;
num_pxl(c) = sum(mask_crop(:));
% title(sprintf('OF divegence for the detection crop\nlog10(SNR) = %2.1f, mean = %2.3f, std = %2.3f', log10(snr_OFdiv(c)), mean_OFdiv(c) ,std_OFdiv(c))) ;
% figure(3)


% subplot(1,2,1)
% hist(OFdiv(:));shg
% title('histogram of OFdiv for bbox')
% [mean(OFdiv(:)), std(OFdiv(:))]
% subplot(1,2,2);
% hist(OFdiv(mask_crop));shg
% title('histogram of OFdiv for mask crop')
% [mean(OFdiv(mask_crop)), std(OFdiv(mask_crop))]
% tilefigs();
% pause(20);
end
%%
close all
bestk = 20;
[~, sorted_minOF_ids] = sort(metric_OFdiv,'ascend');
for c=sorted_minOF_ids(1:bestk)
[bbox, mask_pxl_ids, mask_pxl_coords] = mask2bbox(mask_img{c});
uOF = ppvid.uvOFs{t}(bbox(1):bbox(2), bbox(3):bbox(4),1);
vOF = ppvid.uvOFs{t}(bbox(1):bbox(2), bbox(3):bbox(4),2);

OFdiv = divergence(uOF, vOF);
imshow(abs(OFdiv));shg
mask_crop = (mask_img{c}(bbox(1):bbox(2), bbox(3):bbox(4)));
zz = (OFdiv).*mask_crop;
% zz = abs(OFdiv);
figure(1);
imshow(mask_img{c}); 
figure(2);
imshow(zz); caxis([-1 1]); 
title(sprintf('OF divegence for the detection crop\nmetric_OFdiv = %1.4f, mean = %1.3f, std = %1.3f', metric_OFdiv(c), mean_OFdiv(c) ,std_OFdiv(c))) ;
shg
pause(3);
end
return

%%
close all
figure; hist(mean_OFdiv, -0.5:0.05:1.5); title('histogram of mean of OF divergence');
figure; hist(std_OFdiv, 0:0.025:1.5); title('histogram of std of OF divergence');
figure; plot(std_OFdiv, mean_OFdiv, '.'); xlabel('std'); ylabel('mean'); title('scatter of mean vs std of OF divergence');
% figure; transparentScatter(std_OFdiv.', mean_OFdiv.', 0.005, 0.5); xlabel('std'); ylabel('mean'); title('scatter of mean vs std of OF divergence');
% figure; transparentScatter(std_OFdiv.', mean_OFdiv.', 0.0005, 0.1); xlabel('std'); ylabel('mean'); title('scatter of mean vs std of OF divergence');
% xlim([0, 0.1])
% ylim([-.05 .05])
figure; scatter(std_OFdiv.', mean_OFdiv.', [],log10(num_pxl)); xlabel('std'); ylabel('mean'); title('scatter of mean vs std of OF divergence');
colormap(hsv);
h = colorbar;
ylabel(h, 'log10(num of pixels in region)');
figure; scatter(std_OFdiv.', mean_OFdiv.', [],log10(num_pxl)); xlabel('std'); ylabel('mean'); title('scatter of mean vs std of OF divergence');
colormap(hsv);
h = colorbar;
ylabel(h, 'log10(num of pixels in region)');
xlim([0, 0.1])
ylim([-.05 .05])
return

% mask_uv(:,1) = uOF(mask_pxl_ids);
% mask_uv(:,2) = vOF(mask_pxl_ids);


%%

%% image gradients (edges)
Gmag = double(edge(im, 'canny'));

imshow(Gmag/max(Gmag(:)));shg
d = double(ppvid.depths_pxl{t}(:));

Gmagvec = Gmag(:);

all_pairs = get_inds_of_all_pixels_neighbours([size(im) 1]);
Npairs = length(all_pairs(:,1));
Hdiag = [ones(Nimg,1) ; (1-Gmagvec(all_pairs(:,2))); zeros(Nimg,1)];
Nunknowns = length(Hdiag);
H = sparse(1:Nunknowns, 1:Nunknowns, Hdiag);
f = zeros(Nunknowns,1);

DiDj = sparse(1:Npairs, all_pairs(:,1), ones(Npairs, 1)) + sparse(1:Npairs, all_pairs(:,2), -ones(Npairs, 1));
Sij = sparse(1:Npairs, 1:Npairs, ones(Npairs,1));
Si = sparse(1:Nimg, 1:1:Nimg, -ones(Nimg,1));
Di = sparse(1:Nimg, 1:1:Nimg, ones(Nimg,1));
Aeq = [ [sparse(Npairs, Nimg);Si] [Sij;sparse(Nimg, Npairs)] [ DiDj ; Di] ] ;
beq = [sparse(Npairs,1); d];

%%
tic
X = quadprog(H, f, [], [], Aeq, beq);
toc

%%
figure
imshowpair(reshape(d/max(d), 240, 320), reshape(X((end-Nimg+1):end)/max(X((end-Nimg+1):end)), 240, 320), 'montage');
title(sprintf('given depth estimations                                      after quad prog with edges constraints'), 'Interpreter','Tex');
colormap(flipud(parula));shg
% fname_jpg = sprintf('depth_quad_prog_rho_%2.3f_C_%2.3f.jpg', rho, C);
% saveas(gcf, fname_jpg);

% %% 
% fname_gif = sprintf('depth_quad_prog_rho_%2.3f.gif', rho);
% 
% figure
% imshow(reshape(d/max(d), 240, 320));shg
% xlabel(sprintf('\\rho = %2.3f\ngiven depth estimations', rho),'Interpreter','Tex');
% colormap(flipud(parula));shg
% drawnow
% save_animated_gif_frame(fname_gif, true);
% save_animated_gif_frame(fname_gif, false);
% save_animated_gif_frame(fname_gif, false);
% save_animated_gif_frame(fname_gif, false);
% pause(2)
% imshow(reshape(X/max(X), 240, 320));shg
% xlabel(sprintf('\\rho = %2.3f\nafter quad prog with edges constraints', rho),'Interpreter','Tex');
% colormap(flipud(parula));shg
% drawnow
% save_animated_gif_frame(fname_gif, false);
% save_animated_gif_frame(fname_gif, false);
% save_animated_gif_frame(fname_gif, false);
% save_animated_gif_frame(fname_gif, false);
% 
% 
