function [seg_filtered,seg_filtered_Index,node_filtered,zeroFlowVessels_Index] = filterSegments2(seg,seg_Index,data)
%filterSegments2 removes vessels that have a zero pressure drop due to them
%being n2sn vessels (node 2 same node).
%   INPUTS: seg structure produced by the function sortSegments
%   OUTPUTS: seg structure identical to what is produced by the function
%   sortSegments, but with the n2sn segments filtered out

%declaring coutner for n2sn vessels
n2sn_cntr = 0;
%position counter to avoid extra zero rows
pos_cntr = 1;
%counter for zero flow vessels that aren't n2sn
zero_cntr = 1;
%preallocating index
zeroFlowVessels_Index = [];

for i = 1:numel(seg)
    %filtering out n2sn nodes
    if seg(i).ref(1,1) == seg(i).ref(1,2)
        n2sn_cntr = n2sn_cntr+1;
        continue
    end
    
    %identifying other zero flow regions for refrence when ratios are
    %calculated
    if (numel(data) >= 7 &&(seg(i).pressureDrop == 0)||(seg(i).Flux == 0))
        zeroFlowVessels_Index(zero_cntr,:) = seg(i).ref(1,:);
        zero_cntr = zero_cntr+1;
    end 
   
    
    seg_filtered(pos_cntr) =  seg(i);
    pos_cntr = pos_cntr + 1;
end

seg_filtered_Index = cat(1,seg_filtered.ref);


node_filtered = sortNodes(data,seg_filtered,seg_filtered_Index);
fprintf(strcat('\nNumber of filtered segments:\t',num2str(n2sn_cntr),'\n'))
fprintf('N2SN Filtering:\tComplete\n')
end