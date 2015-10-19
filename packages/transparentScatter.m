function scatterPoints = transparentScatter(x,y,sizeOfCirlce,opacity)
% usage example:
% scatterPoints = transparentScatter(randn(5000,1),randn(5000,1),0.1,0.05);
% set(scatterPoints,'FaceColor',[1,0,0]);

	scale_x = max(x) - min(x);
	scale_y = max(y) - min(y);
    defaultColors = get(0,'DefaultAxesColorOrder');
    assert(size(x,2)  == 1 && size(y,2)  == 1 , 'x and y should be column vectors');
    t= 0:pi/10:2*pi;

    rep_x = repmat(x',[size(t,2),1]);
    rep_y = repmat(y',[size(t,2),1]);
    rep_t = repmat(t',[ 1, size(x,1)]);

    scatterPoints = patch((sizeOfCirlce*sin(rep_t ) * scale_x+ rep_x) ,(sizeOfCirlce*cos(rep_t )* scale_y +rep_y),defaultColors(1,:),'edgecolor','none');
    alpha(scatterPoints,opacity);

end
