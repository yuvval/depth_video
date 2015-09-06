function [fname_depth_gif, fname_frames_gif] = visualize_preproc(vid_name, ppvid)
N = length(ppvid.depths_pxl);
fsmp_str = [ '_fsmp_' num2str(ppvid.frame_sample_interval)];

fname_depth_gif = sprintf('depth_%s%s.gif', vid_name, fsmp_str);
fname_frames_gif = sprintf('frames_%s%s.gif', vid_name, fsmp_str);
fname_frames_gray_gif = sprintf('gray_%s%s.gif', vid_name, fsmp_str);
fname_Gmag_gif = sprintf('Gmag_%s%s.gif', vid_name, fsmp_str);
fname_frames_gray_edges_gif = sprintf('gray_egdes_%s%s.gif', vid_name, fsmp_str);

max_depth = max(max(cell2mat(ppvid.depths_pxl)));
for t=1:N
    figure(1)
    imshow(ppvid.depths_pxl{t}/max_depth);
    colormap(flipud(parula));shg
%     colorbar;shg
    drawnow
    
    save_animated_gif_frame(fname_depth_gif, t==1);
    if t==N % saving last frame again
        save_animated_gif_frame(fname_depth_gif, t==1);
        save_animated_gif_frame(fname_depth_gif, t==1);
    end
    
    figure(2)
    imshow(ppvid.frames{t});
    drawnow
    save_animated_gif_frame(fname_frames_gif, t==1);
    if t==N % saving last frame again
        save_animated_gif_frame(fname_frames_gif, t==1);
        save_animated_gif_frame(fname_frames_gif, t==1);
    end

    figure(3)
    gifname = fname_frames_gray_gif;
    gray_im = rgb2gray(ppvid.frames{t});
    imshow(gray_im);
    drawnow
    save_animated_gif_frame(gifname, t==1);
    if t==N % saving last frame again
        save_animated_gif_frame(gifname, t==1);
        save_animated_gif_frame(gifname, t==1);
    end
    
    figure(4)
    gifname = fname_Gmag_gif;
    gim = imgradient(ppvid.Gx{t}, ppvid.Gy{t});
%     gim = abs(ppvid.Gx{t});
    imshow(gim);
    drawnow
    save_animated_gif_frame(gifname, t==1);
    if t==N % saving last frame again
        save_animated_gif_frame(gifname, t==1);
        save_animated_gif_frame(gifname, t==1);
    end

    figure(5)
    gifname = fname_frames_gray_edges_gif;
    gim = imgradient(ppvid.Gx{t}, ppvid.Gy{t});
    max_gim = max(gim(:));
    gim = 1-gim/max_gim;
%     gim = 1-gim;
    im = double(gray_im).*gim;
    imshow(im/256);
    drawnow
    save_animated_gif_frame(gifname, t==1);
    if t==N % saving last frame again
        save_animated_gif_frame(gifname, t==1);
        save_animated_gif_frame(gifname, t==1);
    end
    
end

