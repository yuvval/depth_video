close all
clear



%% This solver is for a single frame with a QP with Objective & constraints

%%
load preprocessed_videos/approaching_toward_fsmp_15_ppvid.mat

t=6;

imrgb = ppvid.frames{t};
im = double(rgb2gray(ppvid.frames{t}));
Nimg = numel(im);

candidates = im2mcg(imrgb, 'fast');
%%
mask_img = get_mcg_mask_from_candidates(candidates, 5);
[bbox, mask_pxl_ids, mask_pxl_coords] = mask2bbox(mask_img);

uOF = ppvid.uvOFs{t}(bbox(1):bbox(2), bbox(3):bbox(4),1);
vOF = ppvid.uvOFs{t}(bbox(1):bbox(2), bbox(3):bbox(4),2);



divOF = divergence(uOF, vOF);
% imshow(abs(divOF));shg
mask_crop = (mask_img(bbox(1):bbox(2), bbox(3):bbox(4)));
zz = (divOF).*mask_crop;
% zz = abs(divOF);
figure(1);
imshow(mask_img); 
figure(2);
imshow(zz); caxis([-1 1]); shg
figure(3)
hist(divOF(:));shg
[mean(divOF(:)), std(divOF(:))]
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
