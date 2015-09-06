function save_animated_gif_frame(fname, is_first_frame, h_fig)

% adopted from:
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/275433

    if nargin<3
        h_fig = gcf;
    end
    
    %% saving as animated gif
    frame = getframe(h_fig);    
    im = frame2im(frame);    
    [imind,cm] = rgb2ind(im,256);    
    if is_first_frame;        
        imwrite(imind,cm,fname,'gif', 'Loopcount',inf);        
    else        
        imwrite(imind,cm,fname,'gif','WriteMode','append');        
    end
