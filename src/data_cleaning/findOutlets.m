function [Seg_out, Outlet] = findOutlets(node,seg,seg_Index)
%This function takes as inputs the data structure node, produced by the
%function sortNodes and identifies all inlets to the network, by finding
%all nodes that have no inflowing segments and characterising the
%outflowing segments as Inlets
%   INPUTS data strucutre node, see function sortNodes
%   OUTPUTS 2D Arry nx3 containing the refrences of the inlet segments,
%   Inlet a structure conatining the properties of each inlet

cntr  = 1;
for i = 1:numel(node)
    if isempty(node(i).connectionOut)
        for j = 1:numel(node(i).connectionIn)/3
            Seg_out(cntr,:) = node(i).connectionIn(j,:);
            loc = findRowInIndexV(node(i).connectionIn(j,:),seg_Index,3);
            Outlet(cntr) = seg(loc);
            cntr = cntr + 1;
        end
    end
end

for i=1:numel(Outlet)
    loc = findRowInIndexV(Outlet(i).ref,seg_Index,3);
    Outlet(i).locSeg = loc;
end 

    disp('Outlets found')
end 
