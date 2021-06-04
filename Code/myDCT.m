function [transformed] = myDCT(image)
%Inputs:
%image which is the image that is going to be transformed
%Outputs:
%transformed which is the transformed image
%Note: this is not a generic DCT, it is tailored for JPEG operation as it
%splits the input image into 8*8 blocks and perform DCT on each block
[r,c]=size(image);
transformed=zeros(r,c);
for i=1:r/8
    for j=1:c/8
        block=double(image((i-1)*8+1:i*8,(j-1)*8+1:j*8));
        newblock=zeros(8,8);
        for p=0:7
            for q=0:7
                if p==0 && q==0
                    alpha=1/8;
                elseif p==0 || q==0
                    alpha=sqrt(2)/8;
                else
                    alpha=1/4;
                end
                for x=0:7
                    for y=0:7 
                        basis=cos((2*x+1)*p*pi/16)*cos((2*y+1)*q*pi/16);
                        newblock(p+1,q+1)=newblock(p+1,q+1)+block(x+1,y+1)*alpha*basis;
                    end
                end
            end
        end
        transformed((i-1)*8+1:i*8,(j-1)*8+1:j*8)=newblock;
    end
end
end

