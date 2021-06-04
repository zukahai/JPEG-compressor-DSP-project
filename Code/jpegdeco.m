function [decodedimage] = jpegdeco(encoded,dict,allindices,qtable)
%Inputs
% 1- encoded which is the jpeg encoded image
% 2- dict which is the huffman dictionary that was used for encoding (to
%be used in decoding)
% 3- allindices is going to be used to reverse the zigzag pattern (to be
% used in decoding)
% 4- qtable which is the table that was used to quantise the image
%Outputs
% 1- decoded image which is the after decoding
%decode the image using huffman decoding, then Runlength decoding
decoded=huffmandeco(encoded,dict);
decoded=DecodeRL(decoded);
L=length(decoded);
decodedimage=zeros(sqrt(L),sqrt(L));
%use allindices to assign the values of the linear array to their
%appropriate position in the 2d iamge array
j=1;
for i=1:L
    decodedimage(allindices(j),allindices(j+1))=decoded(i);
    j=j+2;
end
%dequantize the image and then perform IDCT
decodedimage=jpeg_dequantize(decodedimage,qtable);
decodedimage=myIDCT(decodedimage);
end

