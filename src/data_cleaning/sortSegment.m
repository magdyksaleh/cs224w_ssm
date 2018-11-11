function [seg, seg_Index] = sortSegment(data)

%sortSegment sorts the structure data into segments and organises them by
%locationn, flux, radius, start and end node, and duplicate effects are
%considered here. 
%   Input is data structure received from the READ_SPATIAL_GRAPH code of an amira file
%   Output is mat structure containing the information sorted by segments

%declaring the position tracker for the number of points
pos_tracker = 0;
%Preallocating the structure
seg(data(2).NumEl).ref(1) = 0;
for i = 1:data(2).NumEl
    seg(i).ref = data(2).Val(i,:);
    seg(i).NumPoints = data(3).Val(i);
    seg(i).StartNode = data(2).Val(i,1);
    %tracking the parameter position across the amira file
    seg(i).SNPos = pos_tracker + 1;
    seg(i).EndNode = data(2).Val(i,2);
    pos_tracker = pos_tracker + seg(i).NumPoints;
    seg(i).ENPos = pos_tracker;
    if numel(data) >= 7
        seg(i).startNodePressure = data(6).Val(seg(i).SNPos);
        seg(i).endNodePressure = data(6).Val(seg(i).ENPos);
        seg(i).pressureDrop = seg(i).endNodePressure - seg(i).startNodePressure;
    end
    
    if numel(data) > 6
        seg(i).Flux = data(7).Val(pos_tracker);
    else
        seg(i).Flux = data(6).Val(pos_tracker);
    end
    %if the pressure drop is negative, the start and end nodes are switched
    %around to represent the direction of flow
    if numel(data) >= 7 && seg(i).pressureDrop < 0
        seg(i).ref = [data(2).Val(i,2), data(2).Val(i,1)];
        temp = seg(i).StartNode;
        seg(i).StartNode = seg(i).EndNode;
        seg(i).EndNode = temp;
        seg(i).pressureDrop = abs(seg(i).pressureDrop);
        seg(i).Flux = abs(seg(i).Flux);
    elseif numel(data) >= 7 && seg(i).pressureDrop == 0
        %seg(i).Flux = 0;
    elseif seg(i).Flux < 0
        seg(i).ref = [data(2).Val(i,2), data(2).Val(i,1)];
        temp = seg(i).StartNode;
        seg(i).StartNode = seg(i).EndNode;
        seg(i).EndNode = temp;
        seg(i).Flux = abs(seg(i).Flux);
    end
end

%  BELOW IS CODE USED TO CALCULATE TIME DELAY
fluxConversion = 1e+6/60; %converting from nL/min to um^3/s

%preallocation
seg(numel(seg)).Radius = 0;
seg(numel(seg)).Velocity = 0;
seg(numel(seg)).Delay = 0;

%defining the centre of the network
net_center = [mean(data(1).Val(:,1)) mean(data(1).Val(:,2)) mean(data(1).Val(:,3))];


for i = 1:numel(seg)
    %finding the distance between the start and end node of the segment
    x1 = [data(4).Val(seg(i).SNPos,:);
        data(4).Val(seg(i).ENPos,:)];
    %distance in micrometers
    seg(i).LengthSimp = pdist(x1);
    length_Array = zeros(seg(i).NumPoints-2+1,1);
    for j = 1:(seg(i).NumPoints)-2+1
        %finding the distance between consecutive points of the segment
        x1 = [data(4).Val(seg(i).SNPos+j-1,:);
            data(4).Val(seg(i).SNPos+j,:)];
        dist = pdist(x1);
        length_Array(j,1) = dist;
    end
    seg(i).Length = sum(length_Array);
    
    %Finding the average radius of the segment in micrometers
    seg(i).Radius = mean([data(5).Val(seg(i).SNPos),data(5).Val(seg(i).ENPos)])/2;
    %Calculating the velocity based on the known flux and cross section
    %Converting the flux to (um^3/s)
    %Velocity in (um/s)
    seg(i).Velocity =(seg(i).Flux*fluxConversion)/((seg(i).Radius)^2*pi);
    %Delay in s
    if seg(i).Velocity == 0
        seg(i).Delay = 0;
    else
        seg(i).Delay= seg(i).Length/seg(i).Velocity;
    end
    seg(i).Volume = pi*(seg(i).Radius)^2*seg(i).Length;
    %Calculating the distance from the center
    startNodeDist = pdist([data(4).Val(seg(i).SNPos,:); net_center]);
    endNodeDist = pdist([data(4).Val(seg(i).ENPos,:); net_center]);
    seg(i).distCenter = mean([startNodeDist,endNodeDist]);
    
    
end

%  END OF CODE USED TO CALCULATE TIME DELAY


%IDENTIFYING DUPLICATE SEGMENTS
%create an index for the segments in the order in which they appear in the
%structure of seg
seg_Index = zeros(numel(seg),2);
for i = 1:numel(seg)
    seg_Index(i,:) = seg(i).ref(1,:);
end

%Here, ic is an array that holds the position of the segments including the
%repititions for repeated segments
[unique_Seg, ~, ic] = unique(seg_Index,'rows');
%identify number of duplicate segments
num_duplicate_segments = length(seg_Index) - length(unique_Seg);
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
duplicate_Index = zeros(numel(seg),1);
for i = 1:length(duplicate_Segments)
    %location identifies the location of duplicate segments in seg_Index
    location =  findRowInIndex(duplicate_Segments(i,:),seg_Index);
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


for i = 1:numel(seg)
    seg(i).duplicate = duplicate_Index(i);
    seg(i).ref = [seg(i).ref duplicate_Index(i)];
end


%Repopulate seg_Index to include the duplicate indicator
seg_Index = cat(1,seg.ref);



disp('Segment Sorting: Complete')



end

