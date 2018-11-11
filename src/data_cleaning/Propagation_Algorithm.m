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
disp('Preprocessing Complete');
