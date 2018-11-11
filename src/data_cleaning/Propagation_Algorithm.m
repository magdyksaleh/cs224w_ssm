%This program calculates propagates through the network calculating a time
% dependant drug distribution

%WRITTEN BY: Magdy Saleh
%DATE LAST EDIT: March 10, 2017

%STATUS: working


clearvars
clc

%LOADING IN THE DATA
filename = '/Users/magdy/Desktop/Stanford/Fall18/224w/project/data/og_files/SW122_spatialGraph_RIN.txt'

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
disp('Segments filtered, zero flow regions identified')

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
    Inlet(i).Scale = Inlet(i).Flux/total_Inflowing_Flux;
end

disp('Inlet Adjustment Complete');

%merging the segment and node information into one structure to eliminate

for i = 1:numel(node_filtered)
    for j = 1:numel(node_filtered(i).connectionOut)/3
        curSeg = node_filtered(i).connection(j,:);
        locSeg = findRowInIndexV(curSeg,seg_filtered_Index,3);
        locScale = findRowInIndexV(curSeg,scale_Index,3);
        locInlet = findRowInIndexV(curSeg,Inlet_Index,3);
        node_filtered(i).Scale(j,1) = scale(locScale).Val;
        node_filtered(i).Flux(j,1) = seg_filtered(locSeg).Flux;
        node_filtered(i).Delay(j,1) = seg_filtered(locSeg).Delay;
        node_filtered(i).Volume(j,1) = seg_filtered(locSeg).Volume;
        node_filtered(i).locSeg(j,1) = locSeg;
        if (~isempty(locInlet)) 
            seg_filtered(locSeg).InletScale(1,1) = Inlet(locInlet).Scale(1,1);
        end
        if isfield(seg_filtered,'pressureDrop')
            node_filtered(i).pressureDrop(j,1) = seg_filtered(locSeg).pressureDrop;
        end
    end
end
disp('Preprocessing Complete');
