cd .\..\..
cd Acquisizioni_iCub
cd(seq)
cd Images
filename1=sprintf('%8.8d.ppm', k);
filename2=sprintf('%8.8d.ppm', k+1);
filename3=sprintf('%8.8d.ppm', k+2);
filename4=sprintf('%8.8d.ppm', k+3);
filename5=sprintf('%8.8d.ppm', k+4);

im1 = double(imread(filename1));			%im1 è una matrice MxNx3 di double
im2 = double(imread(filename2));			%im2 è una matrice MxNx3 di double
im3 = double(imread(filename3));
im4 = double(imread(filename4));
im5 = double(imread(filename5));


%im1gray=mat2gray(im1);  			
%im2gray=mat2gray(im2);
%im3gray=mat2gray(im3);  			
%im4gray=mat2gray(im4);
%im5gray=mat2gray(im5);  			

%&mat2gray:converte una matrice in un’immagine in scala di grigi, quindi
%trasforma la matrice in un’immagine grigia con valori nel range [0,1] ->
%però rimane una matrice MXNX3 (mentre a noi serve ad un solo canale in quanto % ba richiede come input 2 immagini grigie), e quindi convertiamo:
%im1gray=rgb2gray(im1gray);  		
%im2gray=rgb2gray(im2gray);
%im3gray=rgb2gray(im3gray);  		
%im4gray=rgb2gray(im4gray);
%im5gray=rgb2gray(im5gray);  		

%rgb2gray: convert RGB image or colormap to grayscale 

% II(:,:,1)=im1gray;
% II(:,:,2)=im2gray;
% II(:,:,3)=im3gray;
% II(:,:,4)=im4gray;
% II(:,:,5)=im5gray;


II(:,:,1)=im1;
II(:,:,2)=im2;
II(:,:,3)=im3;
II(:,:,4)=im4;
II(:,:,5)=im5;

% %%% 3- Separable filters - Gabor in space (zero DC, energy balanced) + Exp in time%%%
%O(:,:,:) = ctf_population_flow(II,n_scales,th,n_filters,th2,vel);%2 gabor in time
cd .\..\..\..
cd algorithms\Solari_Chessa
O(:,:,:) = ctf_population_flow(II,n_scales,th,th2,vel);%2 gabor in time

uv(:,:,1)=O(:,:,1);
uv(:,:,2)=O(:,:,2);
%figure, imagesc(O(:,:,1)), title('Vx'); %Vx
%figure, imagesc(O(:,:,2)), title('Vy'); %Vy



