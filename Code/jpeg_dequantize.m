function [image] = jpeg_dequantize(quant_image,qtable)
%Inputs
% 1- quant_image which is the quantised image
% 2- qtable which is the table that was used to quantise the image
%Outputs
% 1- image which is the image after dequantization
[r,c]=size(quant_image);
image=zeros(r,c);
%splits the image into 8*8 blocks and multiplies each block element
%by element with the qtable
for i=1:r/8
    for j=1:c/8
        image((i-1)*8+1:i*8,(j-1)*8+1:j*8)=round(quant_image((i-1)*8+1:i*8,(j-1)*8+1:j*8).*qtable);
    end
end
end
