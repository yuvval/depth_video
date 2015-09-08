close all
clear


%%
load preprocessed_videos/approaching_toward_fsmp_15_ppvid.mat

imsize = size(ppvid.frames{1}(:,:,1));
Nframes = length(ppvid.frames);
Nimg = prod(imsize);
N = Nimg*Nframes;
%% image gradients (edges)

[Gmag_frames, depth_frames] = deal(zeros([imsize Nframes]));
[uOFframes, vOFframes] =  deal(zeros([imsize (Nframes-1)]));
for t=1:Nframes
    im = double(rgb2gray(ppvid.frames{t}));
    [ppvid.Gx{t}, ppvid.Gy{t}] = imgradientxy(im, 'IntermediateDifference');
    Gmag_frames(:,:,t) = imgradient(ppvid.Gx{t}, ppvid.Gy{t});
    depth_frames(:,:,t) = ppvid.depths_pxl{t};
    if t<Nframes
        uOFframes(:,:,t) = ppvid.uvOFs{t}(:,:,1);
        vOFframes(:,:,t) = ppvid.uvOFs{t}(:,:,2);
    end
end

%%
[XX0, YY0] = meshgrid(1:imsize(2), 1:imsize(1));
OFmag = zeros([imsize (Nframes-1)]);
OF_ids_map = zeros(Nimg, (Nframes -1));
for t=1:(Nframes-1)
    XX2 = round(XX0 + uOFframes(:,:,t));
    XX2(XX2<1) = 1;
    XX2(XX2>imsize(2)) = imsize(2);
    
    YY2 = round(YY0 + vOFframes(:,:,t));
    YY2(YY2<1) = 1;
    YY2(YY2>imsize(1)) = imsize(1);
    OF_ids_map(:,t) = (t-1)*Nimg + sub2ind(imsize, YY2(:), XX2(:));
    
    OFmag(:,:,t) = sqrt(uOFframes(:,:,t).^2 + vOFframes(:,:,t).^2);
end
%%
rhoOF = 1e-3;
OFweights = exp((-1/rhoOF)*OFmag(:).^2);
Winter = sparse(1:(N-Nimg), OF_ids_map(:), OFweights, N, N);
%%
d = double(depth_frames(:));

Gmag2 = Gmag_frames.^2;
Gmag2vec = Gmag2(:);

all_pairs = get_inds_of_all_pixels_neighbours(size(depth_frames));

%%
rhoE = 0.1;
% imtmp = exp((-1/rhoE)*Gmag2);
% imshow(imtmp/max(imtmp(:)))
all_pairs_edge_weights = exp((-1/rhoE)*Gmag2vec(all_pairs(:,2)));
Wintra = sparse(all_pairs(:,1), all_pairs(:,2), all_pairs_edge_weights, N, N);
W = Wintra + Winter;
Wc = sparse(1:N, 1:N, sum(W,1));
Wr = sparse(1:N, 1:N, sum(W,2));
Wq = sparse(Wr + Wc -2*W +speye(size(W)));

%%
tic
X = quadprog(Wq, -2*d);
toc


%% 
qp_depth_vid = reshape(X/10, size(depth_frames,1), size(depth_frames,2), size(depth_frames,3));
fname_gif = sprintf('depth_video_quad_prog_rhoE_%2.1f_rhoOF_%2.3f.gif', rhoE, rhoOF);

for t=1:Nframes
    imshow(qp_depth_vid(:,:,t));shg
    xlabel(sprintf('\\rhoE = %2.3f, \\rhoOF = %2.3f\nafter quad prog with edges constraints', rhoE, rhoOF),'Interpreter','Tex');
    colormap(flipud(parula));shg
    drawnow
    save_animated_gif_frame(fname_gif, t==1);
end
save_animated_gif_frame(fname_gif, false);
save_animated_gif_frame(fname_gif, false);


