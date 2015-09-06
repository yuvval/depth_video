for  k = 77:85
    
    cd Acquisizioni_iCub\Act_Videos\PointingDiagon\Images    
    filename1=sprintf('%8.8d.ppm', k);
    filename2=sprintf('%8.8d.ppm', k+1);

    im1 = double(imread(filename1));			%im1 è una matrice MxNx3 di double
    im2 = double(imread(filename2));

    im1gray=mat2gray(im1);  			
    im2gray=mat2gray(im2);

    im1gray=rgb2gray(im1gray);  		
    im2gray=rgb2gray(im2gray);
    
    figure
    h=imshow(im1gray);
    title(filename1)
    %&mat2gray:converte una matrice in un’immagine in scala di grigi, quindi
    %trasforma la matrice in un’immagine grigia con valori nel range [0,1] ->
    %però rimane una matrice MXNX3 (mentre a noi serve ad un solo canale in quanto %Lucas_Kanade richiede come input 2 immagini grigie), e quindi convertiamo:
    %rgb2gray: convert RGB image or colormap to grayscale         
    cd .\..\..\..\..
end