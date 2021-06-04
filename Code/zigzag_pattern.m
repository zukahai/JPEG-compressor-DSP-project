function [linear_image,allindices] = zigzag_pattern(block)
%Inputs:
%1- block which is an 8*8 block
%Outputs:
%1- linear_image which is the 8*8 block converted into a linear array using
% a zigzag pattern
%2- allindices which are the actual indices of each element in the linear array in
%the block in order to be able to reverse the zigzag operation
[r,c]=size(block);
linear_image=zeros(1,r*c);
index=1;
index_1=1;
loopflag=-1;%determines whether to loop on columns first or rows first
allindices=zeros(1,2*r*c);
for l=2:16
   break_flag=0;%used to break the loops if we reach the maximum number
 %of elements in each diagonal
 %loops on every diagonal adding each element to the linear array and saves
 %its indices in the 8*8 block in the allindices array
   if (loopflag==1)
       for r=1:8
            for c=1:8
                if r+c==l && break_flag~=l-1
                    linear_image(index)=block(r,c);
                    index=index+1;
                    break_flag=break_flag+1;
                    allindices(index_1)=r;
                    allindices(index_1+1)=c;
                    index_1=index_1+2;
                end
                if(break_flag==l-1)
                    break;
                end
            end
             if(break_flag==l-1)
                break;
            end
       end
       loopflag=-loopflag;
   elseif (loopflag==-1)
       for c=1:8
            for r=1:8
                if r+c==l && break_flag~=l-1
                    linear_image(index)=block(r,c);
                    index=index+1;
                    break_flag=break_flag+1;
                    allindices(index_1)=r;
                    allindices(index_1+1)=c;
                    index_1=index_1+2;
                end
                if(break_flag==l-1)
                    break;
                end
            end
            if(break_flag==l-1)
                    break;
            end
       end
               loopflag=-loopflag;
   end
end
end

