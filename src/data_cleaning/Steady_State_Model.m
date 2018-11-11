%This program calculates the cumulative scaling effect due to the flux
%bifurcation in a vessel network

%WRITTEN BY: Magdy Saleh
%DATE LAST EDIT: Feb 25, 2017

%STATUS: working

clearvars
clc

%LOADING IN THE DATA
filename = 'Sample_File.txt'

[data, definition, vdefinition] = readAmira(filename);

%SORTING THE DATA


%check how many parameters given in data
numParam = numel(data);

%Sort segment by location, flux, radius, and time delay
[seg_unfiltered, seg_Index_unfiltered]= sortSegment(data);

%Filtering out node 2 same node connections
[seg_filtered,seg_filtered_Index, node_filtered,zeroFlowVessels_Index] = filterSegments2(seg_unfiltered,seg_Index_unfiltered,data);

%correcting flows for improved flow conservation 
[node_filtered,seg_filtered,seg_filtered_Index] = correctFlows(node_filtered,seg_filtered,seg_filtered_Index,data);
disp('Segments filtered, zero flow regions identifies')

%Identifing inlets to the network
[Inlet_Index, Inlet] = findInlets(node_filtered,seg_filtered,seg_filtered_Index);
%Identifing outlets to the network
[Outlet_Index, Outlet] = findOutlets(node_filtered,seg_filtered,seg_filtered_Index);

%Find the flow rate factors
[scale, scale_Index] = findFlowRatios(node_filtered,seg_filtered_Index,seg_filtered);

%MP_Constant is the membrance permeability constant across the network
MP_Constant = 1e-5;
seg_filtered = getPermeability(MP_Constant, seg_filtered);

%Summing the total Flux into the network as the sum of the individual inlet
%fluxes
total_Inflowing_Flux = sum(cat(1,Inlet.Flux));

%defining the ratio of fluxes for each inlet as a percentage
for i = 1:numel(Inlet)
    Inlet(i).Scale = 100*(Inlet(i).Flux/total_Inflowing_Flux);
end

disp('Inlet Adjustment Complete');


%Merging structures with extra properties
for i = 1:numel(node_filtered)
    for j = 1:numel(node_filtered(i).connectionOut)/3
        curSeg = node_filtered(i).connection(j,:);
        locSeg = findRowInIndexV(curSeg,seg_filtered_Index,3);
        locInlet = findRowInIndexV(curSeg,Inlet_Index,3);
        if (~isempty(locInlet))
            seg_filtered(locSeg).InletScale(1,1) = Inlet(locInlet).Scale(1,1);
        end
    end
end

disp('Preprocessing Complete');
%END OF PREPROCESSING

tic
num_Unknowns = numel(scale_Index)/3 - numel(Inlet_Index)/3;
A = speye(num_Unknowns);
b = sparse(num_Unknowns,1);
counter = 1;
line_Counter = 1;
%USING AN OUTFLOWING IN TERMS OF INFLOWING CONVENTION

%identifying vessels whose cumulative distribution is not originally known
ref_Unknown = zeros(num_Unknowns,3);
for i = 1:numel(node_filtered)
    if (~isempty(node_filtered(i).connectionOut))&&(~isempty(node_filtered(i).connectionIn))
        for j = 1:numel(node_filtered(i).connectionOut)/3
            ref_Unknown(line_Counter,:) = node_filtered(i).connectionOut(j,:);
            line_Counter = line_Counter +1;
        end
    end
end

fprintf(strcat('Number of Unknown Variables:\t',num2str(length(ref_Unknown)),'\n'))


%Populating the matrix
for i = 1:numel(node_filtered)
    inletArray = zeros(2,1);
    inletScaleArray = inletArray;
    %checking if the node is at the periphery of the network, and skipping
    %it
    if isempty(node_filtered(i).connectionIn)||isempty(node_filtered(i).connectionOut)
        continue
    end
    %checking and identifying which if any of the inflowing segments are
    %known inlets to the network
    for j = 1:numel(node_filtered(i).connectionIn)/3
        inletArray(j,1) = 0;
        inletScaleArray(j,1) = 0;
        inletchecker = findRowInIndexV(node_filtered(i).connectionIn(j,:),Inlet_Index,3);
        if (~isempty(inletchecker))
            inletArray(j,1) = 1;
            inletScaleArray(j,1) = Inlet(inletchecker).Scale;
        end
    end
    %looping through the outflowing segments to start populating the
    %matricies. variable j is purposfully overwritten
    for j = 1:numel(node_filtered(i).connectionOut)/3
        curSegOut = node_filtered(i).connectionOut(j,:);
        loc_curSegOut = findRowInIndexV(curSegOut,ref_Unknown,3);
        loc_Scale = findRowInIndexV(curSegOut,scale_Index,3);
        
        %looping through inlets populating the A matrix
        for k = 1:numel(node_filtered(i).connectionIn)/3
            curSegIn = node_filtered(i).connectionIn(k,:);
            loc_curSegIn = findRowInIndexV(curSegIn,ref_Unknown,3);
            A(loc_curSegOut,loc_curSegIn) = -scale(loc_Scale).Val*not(inletArray(k))*(1-seg_filtered(loc_curSegOut).MP);
        end
        %looping through inlets populating the b vector
        for k = 1:numel(node_filtered(i).connectionIn)/3
            inletScale = scale(loc_Scale).Val*inletScaleArray(k)*inletArray(k)*(1-seg_filtered(loc_curSegOut).MP);
            b(loc_curSegOut,1) = (b(loc_curSegOut,1) + inletScale);
        end
    end
end
disp('Matrix Population Complete')
SOLUTION = A\b;
disp('Matrix Calculation Complete')
toc

% write to data & to seg_filtered the solution
data(numParam+1).Name = 'Cumulative_Dilution';
data(numParam+1).Marker=strcat('@',num2str(numParam+1));
data(numParam+1).Dim=1;
zero_Tracker  = 1;
inlet_Tracker = 1;
pos_Tracker = 0;
missing_seg_tracker = 0;


%adding the numbers in the correct format
for i= 1:length(seg_Index_unfiltered)
    curSeg = seg_Index_unfiltered(i,:);
    numPoints = seg_unfiltered(i).NumPoints;
    test1 = findRowInIndexV(curSeg,ref_Unknown,3);
    if isempty(test1)
        test2 = findRowInIndexV(curSeg,seg_filtered_Index,3);
        if isempty(test2)
            for j = pos_Tracker+1:pos_Tracker+numPoints
                data(numParam+1).Val(j,1) = 0;
            end
            zero_Tracker  = zero_Tracker  +1;
            pos_Tracker = pos_Tracker+numPoints;
            
        else
            for j = pos_Tracker+1:pos_Tracker+numPoints
                data(numParam+1).Val(j,1) = seg_filtered(test2).InletScale;
            end
            inlet_Tracker = inlet_Tracker + 1;
            pos_Tracker = pos_Tracker+numPoints;
        end
    else
        for j = pos_Tracker+1:pos_Tracker+numPoints
            data(numParam+1).Val(j,1) = SOLUTION(test1,:);
        end
        pos_Tracker = pos_Tracker+numPoints;
    end
end
data(numParam+1).NumEl = numel(data(numParam+1).Val);
vdefinition(numParam+1).vname = 'Cumulative_Dilution';
vdefinition(numParam+1).vtype = 'POINT';
vdefinition(numParam+1).vdatatype = 'float';
vdefinition(numParam+1).vdatadim = 1;
vdefinition(numParam+1).vmarker = strcat('@',num2str(numParam+1));
disp('Data write-in complete');
disp('Done!')
