function [decoded] = DecodeRL(str)
%Inputs:
% 1- str which is the runlength encoded stream
%Outputs:
% 1- decoded which is stream after runlength decoding
i=1;
j=1;
while(i<=length(str))
    if(str(i)==0)
        for k=j:j+str(i+1)-1
            decoded(k)=0;
            j=j+1;
        end
       i=i+2;
    else
        decoded(j)=str(i);
        i=i+1;
        j=j+1;
    end
end
end
