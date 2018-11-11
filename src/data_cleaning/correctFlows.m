function [node_filtered,seg_filtered,seg_filtered_Index] = correctFlows(node_filtered,seg_filtered,seg_filtered_Index,data)
%correctFlows corrects the  flow directions in the networks to optimise
%flow conservation
%   INPUTS:  node_filtered - structure, which is the result of the sortNodes
%   function. seg_filtered_Index - matrix, holds the references for all the vessels
%   in the network. seg_filtered- structure, which is the results of the
%   sortSegments function. data - structure, raw data structure from the
%   readAmira function
%   OUTPUTS: Updated versions of the node_filtered, seg_fileterd, and
%   seg_filetered_Index inputs. 


%merging extra parameters into the structure to reduce need for excessive
%searching. 
for i = 1:numel(node_filtered)
    for j = 1:numel(node_filtered(i).connection)/3
        curSeg = node_filtered(i).connection(j,:);
        locSeg = findRowInIndexV(curSeg,seg_filtered_Index,3);
        node_filtered(i).Flux(j,1) = seg_filtered(locSeg).Flux;
        node_filtered(i).locSeg(j,1) = locSeg;
        if isfield(seg_filtered,'pressureDrop')
            node_filtered(i).pressureDrop(j,1) = seg_filtered(locSeg).pressureDrop;
        end
    end
end


%Variable Preallocation
sumFluxIn = 0;
sumFluxOut = 0;
ErrorCheck = zeros(numel(node_filtered),1);


%identifing the difference between flow in and flow out at each node
for i = 1:numel(node_filtered)
    if isempty(node_filtered(i).connectionIn)||isempty(node_filtered(i).connectionOut)
        continue
    end
    for j = 1:numel(node_filtered(i).connectionIn)/3
        loc = findRowInIndexV(node_filtered(i).connectionIn(j,:),seg_filtered_Index,3);
        sumFluxIn = sumFluxIn + seg_filtered(loc).Flux;
    end
    for j = 1:numel(node_filtered(i).connectionOut)/3
        loc = findRowInIndexV(node_filtered(i).connectionOut(j,:),seg_filtered_Index,3);
        sumFluxOut = sumFluxOut + seg_filtered(loc).Flux;
    end
    ErrorCheck(i) = (sumFluxOut - sumFluxIn);
    sumFluxIn = 0;
    sumFluxOut = 0;
end

%setting minimum threshold for nodes where the error is clearly not due to
%working precision differences
ErrorFlag = (abs(ErrorCheck)>1e-4);
ErrorSize = sum(ErrorFlag);
%Identifing nodes with large errors
Error_Node = node_filtered(ErrorFlag);

for i = 1:ErrorSize
    ref = Error_Node(i).ref;
    Error_Node(i).Error = ErrorCheck(ref+1);
    %Flagging edges with zero pressure drop but
    %non-zero flow in the network, as this tends to be an erronous entery
    %in the data and causes a lack of flow conservation
    for j = 1:numel(Error_Node(i).connection)/3
        seg_filtered(Error_Node(i).locSeg(j)).flipRow = (Error_Node(i).pressureDrop(j) == 0);
    end
end

    %Switiching the direction for these edges
if isfield(seg_filtered, 'flipRow')
    for i = 1:numel(seg_filtered)
        if seg_filtered(i).flipRow == 1
            tempVar1 = seg_filtered(i).ref(1,1);
            tempVar2 = seg_filtered(i).StartNode;
            tempVar3 = seg_filtered(i).startNodePressure;
            tempVar4 = seg_filtered(i).SNPos;
            
            seg_filtered(i).ref(1,1) = seg_filtered(i).ref(1,2);
            seg_filtered(i).StartNode = seg_filtered(i).EndNode;
            seg_filtered(i).startNodePressure = seg_filtered(i).endNodePressure;
            seg_filtered(i).SNPos = seg_filtered(i).ENPos;
            
            seg_filtered(i).ref(1,2) = tempVar1;
            seg_filtered(i).EndNode = tempVar2;
            seg_filtered(i).endNodePressure = tempVar3;
            seg_filtered(i).ENPos = tempVar4;
        end
    end
end


% Reanalysing the network for new duplicate segments
%this is due to the changes in naming that occured when the flow direction
%was flipped
seg_filtered_Index = cat(1,seg_filtered.ref);

%Ignoring the previous third values in the edge reference which had refered
%to the presence of a duplicate, but needs to updated
seg_filtered_Index_firsttwo = seg_filtered_Index(:,1:2);
%Here, ic is an array that holds the position of the segments including the
%repititions for repeated segments
[unique_Seg, ~, ic] = unique(seg_filtered_Index_firsttwo,'rows');
%identify number of duplicate segments
num_duplicate_segments = length(seg_filtered_Index_firsttwo) - length(unique_Seg);

%construct an array for duplicate segments locations, preallocating it
findduplicate_Index = zeros(num_duplicate_segments,1);
%rep_index will hold the number of repition for each segment, preallocating
%it
rep_Index = findduplicate_Index;
%defining loop counter that is incremented under certain conditions
counter = 1;
for i = 1:numel(ic)
    %identify the number of times the segment is repeated in the array
    numb_of_repitions = numel(find(ic == i));
    %for each repeated segment
    if numb_of_repitions > 1
        for j = 1:(numb_of_repitions - 1); %1 is subtracted to ignore the first instance of the segment
            %fill rep_index at the corresponding position with the number
            %of times the segment is repeated in the structure
            rep_Index(counter) = (numb_of_repitions);
            %The findduplicate_Index is an array holding the location of
            %the first occurence of segments in seg_Index, that are
            %duplicated throughout seg_Index
            findduplicate_Index(counter,1) = i;
            counter = counter+1;
        end
    end
end

%duplicate_Segments is now indexed to hold all the duplicated segments
duplicate_Segments = unique_Seg(findduplicate_Index,:);
%duplicate_index is a one dimensional array with size equal to the number
%of segments. For every unique segment the corresponding position in the
%duplicate index will hold a zero. If it occurs multiple times then the
%first segment will have a 1 one in its corresponding positon, the second
%occurence will have a 2 etc..
duplicate_Index = zeros(numel(seg_filtered),1);
for i = 1:length(duplicate_Segments)
    %location identifies the location of duplicate segments in seg_Index
    location =  findRowInIndexV(duplicate_Segments(i,:),seg_filtered_Index_firsttwo,2);
    %the lines below are used if you want location to ignore the first
    %occurence of the segment
    %location(1) = 0;
    %location = location(location~=0);
    %         seg(i).duplicateLocation = location;
    rep_Counter = 1;
    for j = 1:numel(location)
        duplicate_Index(location(j)) = rep_Counter;
        rep_Counter = rep_Counter +1;
    end
end


for i = 1:numel(seg_filtered)
    seg_filtered(i).duplicate = duplicate_Index(i);
    seg_filtered(i).ref(1,3) = duplicate_Index(i);
end


%Repopulate seg_Index to include the duplicate indicator
seg_filtered_Index = cat(1,seg_filtered.ref);

node_filtered= sortNodes(data,seg_filtered,seg_filtered_Index);

for i = 1:numel(node_filtered)
    for j = 1:numel(node_filtered(i).connection)/3
        curSeg = node_filtered(i).connection(j,:);
        locSeg = findRowInIndexV(curSeg,seg_filtered_Index,3);
        node_filtered(i).Flux(j,1) = seg_filtered(locSeg).Flux;
        node_filtered(i).locSeg(j,1) = locSeg;
        if isfield(seg_filtered,'pressureDrop')
            node_filtered(i).pressureDrop(j,1) = seg_filtered(locSeg).pressureDrop;
        end
    end
end
end