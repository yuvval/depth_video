function O = ctf_population_flow(II,n_scales,th,th2,vel)

%S4: paper: Image comunication SI

n_frames = size(II,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change Image Size for Pyramid Construction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[sy1 sx1 st] = size(II);

fac = 2^(n_scales-1);

sy2 = ceil(sy1 ./ fac) .* fac; % target resolution
sx2 = ceil(sx1 ./ fac) .* fac; % target resolution

II = [ II ; repmat(II(end,:,:),[sy2-sy1 1 1]) ]; % replicate border row
II = [ II repmat(II(:,end,:),[1 sx2-sx1 1]) ]; % replicate border column


%%%%%%%%%%%%%%%%%
% Image Pyramid %
%%%%%%%%%%%%%%%%%


[II,x_pix,y_pix] = image_pyramid(II,n_frames,n_scales);


%%%%%%%%%%%%%%%%%%%%%%%%%
% Level 1 full velocity %
%%%%%%%%%%%%%%%%%%%%%%%%%

F = filt_gabor_space(II{1});
F = filt_gabor_time(F,vel);

[sy sx n_frames]=size(II{1});

Os= V1_MT(F,II{1},th,th2,vel);

O=Os;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coarse-to-fine Estimation and Merging %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for scale = 2:n_scales
    
    O = expand(O.*2);
    
    
    F = filt_gabor_space(II{scale});
    V = distribute_optic(O,n_frames);
    F = warp_sequence(F,V);
    F{1,1}(isnan(F{1,1}))=0;
    F{1,2}(isnan(F{1,2}))=0;
    F = filt_gabor_time(F,vel);

    Os= V1_MT(F,II{scale},th,th2,vel);
    clear F;
    
    O = merge_flow(O,Os);
    
end



% Remove all flow that has not been updated (confirmed) on the lowest
% scale


IND = isnan(sum(Os,3));%filled
O(cat(3,IND,IND)) = NaN;

% Remove rows and columns that were added to construct the pyramid

O(end-(sy2-sy1-1):end,:,:) = [];
O(:,end-(sx2-sx1-1):end,:) = [];









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [II,x_pix,y_pix] = image_pyramid(II,n_frames,n_scales)


[sy sx st] = size(II);

x_pix = cell(1,n_scales);
y_pix = cell(1,n_scales);

[ x_pix{n_scales} y_pix{n_scales} ] = meshgrid(1:sx,1:sy);

lpf = [1 4 6 4 1]/16;

tmp = II;
II = cell(1,n_scales);
II{n_scales} = tmp;

for scale = n_scales-1:-1:1
    for frame = 1:n_frames
        tmp(:,:,frame) = conv2b(conv2b(tmp(:,:,frame),lpf,3),lpf',3);
    end
    [Ny Nx dummy] = size(tmp);

    tmp = tmp(1:2:Ny,1:2:Nx,:);
    II{scale} = tmp;
    x_pix{scale} = x_pix{scale+1}(1:2:Ny,1:2:Nx);
    y_pix{scale} = y_pix{scale+1}(1:2:Ny,1:2:Nx);
end



%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


function O = expand(O)


sy = size(O,1);
sx = size(O,2);

[X Y] = meshgrid(1:(sx-1)/(2*sx-1):sx, ...
    1:(sy-1)/(2*sy-1):sy);


% Repeat edge pixel for border handling

O = [ O O(:,end,:) ];
O = [ O ; O(end,:,:) ];

O = bilin_interp(O,X,Y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function G_warped = warp_sequence(G,V)


[sy sx n_orient n_frames] = size(G{1});
G_warped = cell(1,2);

[X Y] = meshgrid(1:sx,1:sy);
tmp = V;
tmp(isnan(tmp)) = 0;
Vx = tmp(:,:,:,1);
Vy = tmp(:,:,:,2) ;

for frame = 1:n_frames

    Xn = X-Vx(:,:,frame);
    Yn = Y-Vy(:,:,frame);
    G_warped{1}(:,:,:,frame) = (bilin_interp((squeeze(G{1}(:,:,:,frame))),Xn,Yn));
    G_warped{2}(:,:,:,frame) = (bilin_interp((squeeze(G{2}(:,:,:,frame))),Xn,Yn));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function O = merge_flow(O1,O2)


invalid1 = isnan(sum(O1,3));
invalid2 = isnan(sum(O2,3));

O1(cat(3,invalid1,invalid1)) = 0;
O2(cat(3,invalid2,invalid2)) = 0;

invalid = invalid1 & invalid2;

textured=0;
if ~textured
    invalid = invalid1 | invalid2;
end

O = O1 + O2;
O(cat(3,invalid,invalid)) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function V = distribute_optic(O,n_frames)

[sy sx dummy]=size(O);
V = zeros(sy,sx,n_frames,2);
O(isnan(O)) = 0;

sh=[2 1 0 -1 -2];
for ii=1:n_frames
    V(:,:,ii,1) =O(:,:,1)*sh(ii);   V(:,:,ii,2) = O(:,:,2)*sh(ii);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [O]= V1_MT(F,II,th,th2,vel)

n_vel=size(F,2);


[sy sx n_orient]=size(F{1,1});
of=ones(sy,sx)*NaN;
for n_v=1:n_vel
    E(:,:,:,n_v)= sqrt(F{1,n_v}.^2+F{2,n_v}.^2);
end

E=E.^0.5;

E=E/max(max(max(max(E))));
mask=E>th;
E=E.*mask;

O = NaN(sy, sx, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%normalization: tmp matrix used later
tmp=zeros(sy,sx,n_vel);
for v=1:n_vel
    for o=1:n_orient
        tmp(:,:,v)=tmp(:,:,v)+E(:,:,o,v);
    end
end


filt = fspecial('gaussian', [5 5], 10/6);
filt=filt./sum(sum(fspecial('gaussian', [5 5], 10/6)));


for v=1:n_vel
    for o=1:n_orient
        E_v1(:,:,o,v)=E(:,:,o,v)./(tmp(:,:,v)+1e-9);%normalization
        E_v1(:,:,o,v) = conv2(E_v1(:,:,o,v), filt, 'same');%V1 spatial pooling
        mask=ones(sy,sx);
        mask(tmp(:,:,v)<th2)=NaN;%threshold (unreliable pixels)
        E_v1(:,:,o,v)=E_v1(:,:,o,v).*mask;
    end
end

E_v1=0.25*E_v1/max(max(max(max(E_v1))));%exponential gain

for v=1:n_vel
    csum=zeros(sy,sx);
    ssum=zeros(sy,sx);
    for o=1:n_orient
        theta=(o-1)*(2*pi/n_orient);
        csum=csum + cos(theta)*E_v1(:,:,o,v);%MT: orientation pooling
        ssum=ssum + sin(theta)*E_v1(:,:,o,v);
    end
    E_mtc(:,:,v)=exp(csum);%MT non-linearity
    E_mts(:,:,v)=exp(ssum);

end

vx=zeros(sy,sx);
vy=zeros(sy,sx);


    %%%%%%%%%interpolation: filling-in of borders and unreliable pixels
    MMi=max(max(max(II(:,:,3))));
    mmi=min(min(min(II(:,:,3))));
    [xx1,xx2,xx3]=size(E_mtc(:,:,1));
    if xx1<20 && (~isequal(isnan(E_mtc),zeros(size(E_mtc))) ||~isequal(isnan(E_mts),zeros(size(E_mts))))
        for v=1:n_vel
            Os(:,:,1)=E_mtc(:,:,v);
            Os(:,:,2)=E_mts(:,:,v);
           % Os=myfillin_ppp(Os,II(:,:,3),(MMi-mmi)/6);
            Os=myfillin_pp(Os);%,II(:,:,3),(MMi-mmi)/6);
            E_mtc(:,:,v)=Os(:,:,1);
            E_mts(:,:,v)=Os(:,:,2);
        end
    end

    if xx1>=20
        for v=1:n_vel
            Os(:,:,1)=E_mtc(:,:,v);
            Os(:,:,2)=E_mts(:,:,v);
            Os(1:5,:,:)=NaN; Os((end-4):end,:,:)=NaN;  Os(:,1:5,:)=NaN; Os(:,(end-4):end,:)=NaN;
            %Os=myfillin_ppp(Os,II(:,:,3),(MMi-mmi)/6);
            Os=myfillin_pp(Os);
            E_mtc(:,:,v)=Os(:,:,1);
            E_mts(:,:,v)=Os(:,:,2);
        end
    end
    



%%%%%%MT decoding
for v=1:n_vel
    vx=vx+E_mtc(:,:,v)*vel(v);
    vy=vy+E_mts(:,:,v)*vel(v);

end

O(:,:,1)=vx;
O(:,:,2)=vy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IF = filt_gabor_space(I)

DC_thr=0;
n_filters=12;

[nr,nc,n_frames] = size(I);

IF{1} = zeros(nr,nc,n_filters,n_frames);
IF{2} = zeros(nr,nc,n_filters,n_frames);


w=5;
f0=1/3.8;
sigma_s=1.2*sqrt(2*log(2)) ./ (2*pi*(f0/3));


[X,Y] = meshgrid(-w:w,-w:w);
theta=0:2*pi/n_filters:(2*pi-2*pi/n_filters);

G = exp(-(X.^2+Y.^2)/(2*sigma_s^2));

for ii=1:length(theta)
    XT=cos(theta(ii))*X+sin(theta(ii))*Y;
    GC=G.*cos(2*pi*f0*XT);
    GCB{ii}=GC-sum(sum(GC))/(2*w+1)^2;%DC
    GS=G.*sin(2*pi*f0*XT);
    GSB{ii}=GS-sum(sum(GS))/(2*w+1)^2;%DC
end


for frame = 1:n_frames

    for ii=1:n_filters/2        
        even=conv2b(I(:,:,frame),GCB{ii},3);
        odd=conv2b(I(:,:,frame),GSB{ii},3);
        
        IF{1}(:,:,ii,frame) = even;
        IF{1}(:,:,ii+n_filters/2,frame) = even;
        
        IF{2}(:,:,ii,frame) = odd;
        IF{2}(:,:,ii+n_filters/2,frame) = -odd;
    end
 
end



% invalid = (abs(real(IF))<DC_thr) | (abs(imag(IF))<DC_thr);
% IF(invalid) = NaN;

% FI_n{1}=real(IF);
% FI_n{2}=imag(IF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function IF = filt_gabor_time(F,v)


n_vel=size(v,2);
[nr nc n_orient n_frames]=size(F{1});

for n_v=1:n_vel
    
    w0=-v(n_v)/3.8;
    
    %     t=0:n_frames-1;
    %     f0=1/3.8;
    %     B=f0/2.5;
    %     sigma = sqrt(2*log(2)) ./ (2*pi*B);
    %     g = exp(-t ./ (2.*sigma.^2));
    %     Fts = g.*sin(2*pi*w0*t);
    %     Ftc = g.*cos(2*pi*w0*t);
    %     Ftc=Ftc';
    %     Fts=Fts';
    
        t=0:n_frames-1;
        f0=1/3.8;
        B=f0/2.5;
        sigma = sqrt(2*log(2)) ./ (2*pi*B);
        g = exp(-t ./ (2.*sigma.^2));
        Fts = g.*sin(2*pi*w0*t);
        Ftc = g.*cos(2*pi*w0*t);
        Ftc=Ftc';
        Fts=Fts';
    
    
    %     t=0:n_frames-1;
    %     f0=1/3.8;
    %     sigma=1.2*sqrt(2*log(2)) ./ (2*pi*(f0/3));
    %     g = exp(-t.^2 ./ (2.*sigma.^2));
    %     Fts = g.*sin(2*pi*w0*t);
    %     Ftc = g.*cos(2*pi*w0*t);
    %
    %     Ftc=Ftc';
    %     Fts=Fts';
    
    
    
    G_even_tmp=F{1};
    G_odd_tmp=F{2};
    
    G_even3d=zeros(nr,nc,n_orient);
    G_odd3d=zeros(nr,nc,n_orient);
    
    for orient=1:n_orient
        for i=1:nr
            
            G_even3d(i,:,orient) = (conv2(squeeze(G_even_tmp(i,:,orient,:))',Ftc,'valid')-conv2(squeeze(G_odd_tmp(i,:,orient,:))',Fts,'valid'))';
            G_odd3d(i,:,orient) = (conv2(squeeze(G_even_tmp(i,:,orient,:))',Fts,'valid')+conv2(squeeze(G_odd_tmp(i,:,orient,:))',Ftc,'valid'))';
            
        end
    end
    
    IF{1,n_v}= squeeze(G_even3d);
    IF{2,n_v} = squeeze(G_odd3d);
end












