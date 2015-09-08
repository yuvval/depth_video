close all
clear


%%
load preprocessed_videos/approaching_toward_fsmp_15_ppvid.mat

t=6;
im = double(rgb2gray(ppvid.frames{t}));
%% image gradients (edges)
% [ppvid.Gx{t}, ppvid.Gy{t}] = imgradientxy(im, 'IntermediateDifference');
% 
% [Gmag, Gdir] = imgradient(ppvid.Gx{t}, ppvid.Gy{t});
Gmag =0* double(edge(im, 'canny')) +10;
imshow(Gmag/max(Gmag(:)));shg
% return
d = double(ppvid.depths_pxl{t}(:));
d = d/mean(d(:));
% d = double(ppvid.depths_superpxl{t}(:));

Gmag2 = Gmag.^2;

Gmag2vec = Gmag2(:);

all_pairs = get_inds_of_all_pixels_neighbours([size(im) 1]);

%%
% rho = 0.01;
rho = 10;
% imtmp = exp((-1/rho)*Gmag2);
% imshow(imtmp/max(imtmp(:)))

N = length(d);
% C = mean(d(:));
C = 1;
C = 1/sqrt(rho);
all_pairs_edge_weights = C*exp((-1/rho)*Gmag2vec(all_pairs(:,2)));
Wintra = sparse(all_pairs(:,1), all_pairs(:,2), all_pairs_edge_weights);
W = Wintra;
Wc = sparse(1:N, 1:N, sum(W,1));
Wr = sparse(1:N, 1:N, sum(W,2));
Wq = sparse(Wr + Wc -2*W +speye(size(W)));

%%
tic
X = quadprog(Wq, -2*d);
toc

%%
figure
imshowpair(reshape(d/max(d), 240, 320), reshape(X/max(X), 240, 320), 'montage');
title(sprintf('given depth estimations                                      \\rho = %2.3f, C = %2.3f, after quad prog with edges constraints', rho, C), 'Interpreter','Tex');
colormap(flipud(parula));shg
fname_jpg = sprintf('depth_quad_prog_rho_%2.3f_C_%2.3f.jpg', rho, C);
saveas(gcf, fname_jpg);

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
