function [ location ] = findRowInIndexV(row,index,dim) 
%findRowInIndexV takes inputs of the segment and the reference and returns the
%location of the row in the index. (Faster version of ismember). Allows for
%Variable dimensions
   %    Inputs are the row, index, and the dimensions of the row and index 

location = (find(any(all(bsxfun(@eq,reshape(row.',1,dim,[]),index),2),dim)));

end

