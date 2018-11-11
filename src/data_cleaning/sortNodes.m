function nodes = sortNodes(data,seg,seg_Index)
%sortNodes takes in the n by 2 seg and seg_Index data. includes duplicate
%analysis
%find the length of the segment index that holds all the vessel information
len = length(seg_Index);
%looking mainly at the index based on start and end node
seg_Index_ref = seg_Index(:,1:2);

%define the center of the network as the median of the coordinates

net_center = [median(data(1).Val(:,1)) median(data(1).Val(:,2)) median(data(1).Val(:,3))];

%Preallocating structure to hold interim node information
b(data(1).NumEl).connection(1) = 1;
%looping through number of nodes
for i = 1:data(1).NumEl
    %adjusting for zero indexing
    nodes(i).ref = i-1;
    %storing the coordinates of the node
    nodes(i).coordinates = data(1).Val(i,:);
    
    %DEPTH FUNCTION
    %Measures the euclidian distance of a point from the center of the network
    tempVar = [nodes(i).coordinates;net_center];
    nodes(i).distCenter = pdist(tempVar);
    
    %finding number of connection that involve the aforementioned node
    b(i).connection = find(seg_Index_ref == (i-1));
    %declaring counter
    bcounter = 1;
    for j = 1:length(b(i).connection)
        %this condition checks if the occurence of the node is in the first
        %or second columns of the index, if in first, then it is initially
        %an outbound flowing connection and vice versa.
        if b(i).connection(j) <= len
            %since maltab indexes fields in a matrix as a long 1d array
            %stacked on top of each other, the length is used to jump
            %between columns
            nodes(i).connection(j,1) = seg_Index(b(i).connection(j));%first column entry
            nodes(i).connection(j,2) = seg_Index(b(i).connection(j)+len);%second column entry
            nodes(i).connection(j,3) = seg_Index(b(i).connection(j)+2*len);%third column entry
            nodes(i).connectionOut(j,1) = seg_Index(b(i).connection(j));%first column entry
            nodes(i).connectionOut(j,2) = seg_Index(b(i).connection(j)+len);%second column entry
            nodes(i).connectionOut(j,3) = seg_Index(b(i).connection(j)+2*len);%third column entry
        elseif and(b(i).connection(j) >= len,b(i).connection(j) <= 2*len)
            nodes(i).connection(j,1) = seg_Index(b(i).connection(j)-len);%first column entry
            nodes(i).connection(j,2) = seg_Index(b(i).connection(j));%second column entry
            nodes(i).connection(j,3) = seg_Index(b(i).connection(j)+len);%third column entry
            nodes(i).connectionIn(bcounter,1) = seg_Index(b(i).connection(j)-len);%first column entry
            nodes(i).connectionIn(bcounter,2) = seg_Index(b(i).connection(j));%second column entry
            nodes(i).connectionIn(bcounter,3) = seg_Index(b(i).connection(j)+len);%third column entry
            bcounter = bcounter+1;
        end
    end
end
disp('Connection & Node Sorting: Complete1')
pos = 1;
%finding the position of nodal information based on amira's segment-centric
%file write
for i = 1:numel(seg)
    nodeloc1(i).ref = seg(i).StartNode;
    nodeloc1(i).pos = pos;
    nodeloc2(i).ref = seg(i).EndNode;
    pos = pos + seg(i).NumPoints - 1;
    nodeloc2(i).pos = pos;
    pos = pos + 1;
end
%position locator needs to be reviewed, and code needs commenting
A = vertcat(nodeloc1.ref,nodeloc2.ref);
B = vertcat(nodeloc1.pos,nodeloc2.pos);
nodeloc = horzcat(A,B);
b(data(1).NumEl).pos = 1;
for i = 1:data(1).NumEl;
    b(i).pos = find(A == (i-1));
    for j = 1:length(b(i).pos);
        nodes(i).position(j,1) = B(b(i).pos(j));
    end
end


end

