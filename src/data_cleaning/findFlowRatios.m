function [scale, scale_Index] = findFlowRatios(node,seg_Index,seg)
%findFlowRatios calculated the ratio of the flow of each daughter vessel to
%the flow of its mother vessel.
%   INPUTS:  node - structure, which is the result of the sortNodes
%   function. seg_Index - matrix, holds the references for all the vessels
%   in the network. seg - structure, which is the results of the
%   sortSegments function. 
%   OUTPUTS: scale - array containing the flow ratios. scale_Index -
%   matrix, containing the references to the segments in the order that
%   they appear

%variable preallocation
sumFlux = 0;
cntr = 1;



for i  = 1:numel(node)
    if isempty(node(i).connectionIn)
        %defining the scale for the inlet vessels 
        for j = 1:numel(node(i).connectionOut)/3
            scale(cntr).Ref = node(i).connectionOut(j,:);
            scale(cntr).Val = 1;
            cntr = cntr +1;
        end
    else
        %calculating the sum of flows in the node
        for j = 1:numel(node(i).connectionIn)/3
            curSeg = node(i).connectionIn(j,:);
            locSeg = findRowInIndexV(curSeg,seg_Index,3);
            sumFlux = seg(locSeg).Flux + sumFlux;
        end
        %calculating the flow ratio (scale) 
        for j = 1:numel(node(i).connectionOut)/3
            curSeg = node(i).connectionOut(j,:);
            locSeg = findRowInIndexV(curSeg,seg_Index,3);
            scale(cntr).Ref = node(i).connectionOut(j,:);
            if sumFlux == 0
                scale(cntr).Val = 0;
            else
                scale(cntr).Val = seg(locSeg).Flux/sumFlux;
            end
            cntr = cntr +1;
        end
    end
    sumFlux = 0;
end
scale_Index = cat(1,scale.Ref);

disp('Ratios found')
end