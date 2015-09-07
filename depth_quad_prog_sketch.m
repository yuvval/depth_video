close all
clear


%%
load preprocessed_videos/approaching_toward_fsmp_15_ppvid.mat

t=7;
im = rgb2gray(ppvid.frames{t});

[Gmag, Gdir] = imgradient(ppvid.Gx{t}, ppvid.Gy{t});
d = double(ppvid.depths_pxl{t}(:));

Gmag2 = Gmag.^2;
Gmag2vec = Gmag2(:);
% Gmag2vec(Gmag2vec<1.5) = 0;
% Gmag2vec(Gmag2vec>=1.5) = 1;

all_pairs = get_inds_of_all_pixels_neighbours(size(im));
%%
rho = 1;
N = length(d);
all_pairs_edge_weights = exp((-1/rho)*Gmag2vec(all_pairs(:,2)));
W = sparse(all_pairs(:,1), all_pairs(:,2), all_pairs_edge_weights);
Wc = sparse(1:N, 1:N, sum(W,1));
Wr = sparse(1:N, 1:N, sum(W,2));
Wq = sparse(Wr + Wc -2*W +speye(size(W)));

%%
tic
X = quadprog(Wq, -2*d);
toc


%% 
close all
fname_gif = sprintf('depth_quad_prog_rho_%2.3f.gif', rho);

imshow(reshape(d/max(d), 240, 320));shg
xlabel(sprintf('\\rho = %2.3f\ngiven depth estimations', rho),'Interpreter','Tex');
colormap(flipud(parula));shg
drawnow
save_animated_gif_frame(fname_gif, true);
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);

imshow(reshape(X/max(X), 240, 320));shg
xlabel(sprintf('\\rho = %2.3f\nafter quad prog with edges constraints', rho),'Interpreter','Tex');
colormap(flipud(parula));shg
drawnow
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);


