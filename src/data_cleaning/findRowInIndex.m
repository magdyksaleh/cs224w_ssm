function [ location ] = findRowInIndex(row,index) 
%findRowInIndex takes inputs of the segment and the reference and returns the
%location of the row in the index. (Faster version of ismember)

location = (find(any(all(bsxfun(@eq,reshape(row.',1,2,[]),index),2),3)));

end

