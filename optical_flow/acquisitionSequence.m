cd videos

im1=video(:,:,:,k);
im2=video(:,:,:,k+frame_sample_interval);

im1gray=mat2gray(im1);  			
im2gray=mat2gray(im2);

im1gray=rgb2gray(im1gray);  		
im2gray=rgb2gray(im2gray);

close all
cd ./..