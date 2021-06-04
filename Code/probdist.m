function [symbols,prob] = probdist(source)
%Inputs:
%source which is the sample input
%Outputs:
%symbols which is an array of the all the unique characters in source
%prob which is an array whose elements are the probability of occurence 
%of each character of symbols
[~,indices]=unique(source);
L=length(indices);
symbols=zeros(1,L);
prob=zeros(1,L);
for i=1:length(indices)
    symbols(i)=source(indices(i));
    prob(i)=mean(source==source(indices(i)));
end
end

