function [image] = myIDCT(transformed)
%Inputs:
%transformed which is the transformed image
%Outputs:
%image which is the image after reversing the DCT transform
%Note: this is not a generic IDCT, it is tailored for JPEG operation as it
%splits the input transformed image into 8*8 blocks and perform IDCT on each block
[r,c]=size(transformed);
image=uint8(zeros(r,c));
for i=1:r/8
    for j=1:c/8
        block=transformed((i-1)*8+1:i*8,(j-1)*8+1:j*8);
        newblock=zeros(8,8);
        for x=0:7
            for y=0:7
                for p=0:7
                    for q=0:7 
                        if p==0 && q==0
                            alpha=1/8;
                        elseif p==0 || q==0
                            alpha=sqrt(2)/8;
                        else
                            alpha=1/4;
                        end
                        basis=cos((2*x+1)*p*pi/16)*cos((2*y+1)*q*pi/16);
                        newblock(x+1,y+1)=newblock(x+1,y+1)+block(p+1,q+1)*alpha*basis;
                    end
                end
            end
        end
        image((i-1)*8+1:i*8,(j-1)*8+1:j*8)=uint8(newblock);
    end
end
end

