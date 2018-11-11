function [Seg_in, Inlet] = findInlets(node,seg,seg_Index)
%This function takes as inputs the data structure node, produced by the
%function sortNodes and identifies all inlets to the network, by finding
%all nodes that have no inflowing segments and characterising the
%outflowing segments as Inlets
%   INPUTS data strucutre node, see function sortNodes
%   OUTPUTS 2D Arry nx3 containing the refrences of the inlet segments,
%   Inlet a structure conatining the properties of each inlet

cntr  = 1;
for i = 1:numel(node)
    if isempty(node(i).connectionIn)
        for j = 1:numel(node(i).connectionOut)/3
            Seg_in(cntr,:) = node(i).connectionOut(j,:);
            loc = findRowInIndexV(node(i).connectionOut(j,:),seg_Index,3);
            Inlet(cntr) = seg(loc);
            cntr = cntr + 1;
        end
    end
end
    disp('Inlets found')
end 
