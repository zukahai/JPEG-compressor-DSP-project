%%
%reads the image
image = imread('Einstein.jpg','jpeg');
%% 
%define the quantization tables
qtable1=[1 1 1 1 1 2 2 4;1 1 1 1 1 2 2 4;1 1 1 1 2 2 2 4,...
;1 1 1 1 2 2 4 8;1 1 2 2 2 2 4 8;2 2 2 2 2 4 8 8;...
2 2 2 4 4 8 8 16;4 4 4 4 8 8 16 16];
qtable2=[1 2 4 8 16 32 64 128;2 4 4 8 16 32 64 128;...
4 4 8 16 32 64 128 128;8 8 16 32 64 128 128 256;...
16 16 32 64 128 128 256 256;32 32 64 128 128 256 256 256;...
64 64 128 128 256 256 256 256;128 128 128 256 256 256 256 256];
%% 
%encode and decode the image using the first quantization table then save
%it to a file
[encoded,dict,allindices]=jpegenco(image,qtable1);
compressed_image_q1=jpegdeco(encoded,dict,allindices,qtable1);
imshow(compressed_image_q1);
imwrite(compressed_image_q1,'Compressed_Image_q1.jpg','jpeg');
%% 
%encode and decode the image using the second quantization table then save
%it to a file
[encoded,dict,allindices]=jpegenco(image,qtable2);
compressed_image_q2=jpegdeco(encoded,dict,allindices,qtable2);
imshow(compressed_image_q2);
imwrite(compressed_image_q2,'Compressed_Image_q2.jpg','jpeg');

