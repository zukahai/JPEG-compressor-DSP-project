function [encoded,dict,allindices] = jpegenco(image,qtable)
%
% Inputs:
% 1- image which is the image that is going to be compressed
% 2- qtable which is the quantisation table that is going to be used
% Outputs:
% 1- encoded which is the jpeg encoded image
% 2- dict which is the huffman dictionary that was used for encoding (to
%be used in decoding)
% 3- allindices is going to be used to reverse the zigzag pattern (to be
% used in decoding)
%transform and quantise the image using myDCT 
transformed=myDCT(image);
transformed_q=jpeg_quantize(transformed,qtable);
[r,c]=size(transformed_q);
encoded=0;
allindices=0;
%Split the transformed image into 8*8 blocks then convert each 8*8 block
%into a linear array using a zigzag pattern
%encode each linear array using run-length encoding
for i=1:r/8
    for j=1:c/8
        block=transformed_q((i-1)*8+1:i*8,(j-1)*8+1:j*8);
        [temp,temp1]=zigzag_pattern(block);
        for k=1:2:length(temp1)-1
            temp1(k)=temp1(k)+(i-1)*8;
            temp1(k+1)=temp1(k+1)+(j-1)*8;
        end
        encoded=horzcat(encoded,EncodeRL(temp));
        allindices=horzcat(allindices,temp1);
    end
end
%Encode the Runlength encoded linear array using Huffman encoding
%use the probdist implemented function in order to find the distribution of
%the symbols of the linear array
encoded(1)=[];
allindices(1)=[];
[symbols,prob]=probdist(encoded);
dict=huffmandict(symbols,prob);
encoded=huffmanenco(encoded,dict);
end

