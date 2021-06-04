function [quant_image] = jpeg_quantize(image,qtable)
%Inputs
% 1- image which is the image that is going to be quantised
% 2- qtable which is the table that will be used to quantise the image
%Outputs
% 1- quant_image which is the image after quantization
[r,c]=size(image);
quant_image=zeros(r,c);
%splits the image into 8*8 blocks and divides each block element
%by element with the qtable
for i=1:r/8
    for j=1:c/8
        quant_image((i-1)*8+1:i*8,(j-1)*8+1:j*8)=round(image((i-1)*8+1:i*8,(j-1)*8+1:j*8)./qtable);
    end
end
end

