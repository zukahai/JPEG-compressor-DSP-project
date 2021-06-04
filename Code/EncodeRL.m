function [out] = EncodeRL(str)
%Inputs:
%str which is the string that is going to be encoded
%Outputs:
%out which is the runlength encoded stream
counter=0;
j=1;
for i=1:length(str)
  if(str(i)==0 && counter~=9)
      counter=counter+1;
  elseif(str(i)~=0 &&counter~=0)
      out(j)=0;
      out(j+1)=counter;
      out(j+2)=str(i);
      j=j+3;
      counter=0;
  elseif(str(i)~=0 &&counter==0)
      out(j)=str(i);
      j=j+1;
  end
  if(counter==9)
      out(j)=0;
      out(j+1)=counter;
      counter=0;
      j=j+2;
  end
if(i==length(str) && counter~=0)
        out(j)=0;
        out(j+1)=counter;
    end
end


