%This program calculates propagates through the network calculating a time
% dependant drug distribution

%WRITTEN BY: Magdy Saleh
%DATE LAST EDIT: March 10, 2017

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

%END OF PREPROCESSING
%DEFINING TIME SCALE
t = 0:0.1:60;
C_Func_select = 1;
C_in = cfunction(t,C_Func_select);

%Define threshold value for propagation
ignore_Threshold = 1e-7;

Concentration_matrix = zeros(numel(seg1),length(t));
cntr2 = 1;
cntr3 = 1;
tic

%iterating through all inlets
for k = 1:numel(Inlet_Index)/3
    %current Inlet
    curSegIn = Inlet_Index(k,:);
    %current Segment
    curSeg = curSegIn;
    BreakCheck = 2;
    cntr = 1;
    
    %Defining the arrays for cumulative analysis of permeability, scale and
    %delay
    CumPermeabilityInitial = [1 1];
    CumFluxScaleInitial = [1 1];
    CumDelayInitial = [0 0];
    locationTracker = [1 0 0];
    
    %This algorithm works by iteratin thorough two steps according to a
    %propagating front, where the current generation of edges considered 
    %and their properties are defined. Then the next generation of the 
    %edges and their properties is defined. After every iteration the next
    %current generation is replaced with the next generation and it's next
    %generation defined etc..
    current_Nodes = [];
    next_Segments =[];
    Concentration_matrix_temp = zeros(numel(seg1),length(t));
    while BreakCheck ~= 0
        P_Counter = 1;
        if curSeg(1,:) == curSegIn
            loc = findRowInIndexV(curSeg,seg_filtered_Index,3);
            %if inlet, just add the input function to the inlet row
            CumFluxScaleInitial = seg_filtered(loc).InletScale;
            CumDelay = 0;
            CumPermeabilityInitial = (1-seg_filtered(loc).MP);
            Concentration_matrix_temp(loc,:) = CumFluxScaleInitial.*C_in*CumPermeabilityInitial;
        end
        %storing all the final nodes for all the edges in curSeg, one is
        %added to adjust for MATLAB's non zero indexing
        for i = 1:numel(curSeg)/3
            current_Nodes(i,1) = curSeg(i,2) + 1;
        end
        
        for i = 1:numel(curSeg)/3
            ref = curSeg(i,:);
            for j = 1:numel(node_filtered(ref(1,1)+1).connectionOut)/3
                %this loop is used to locate the segment property from the
                %node file, by looping through the node's outflowing edges
                %and identifying the one that corresponds to the currently
                %considered edge
                if (node_filtered(ref(1,1)+1).connectionOut(j,:) == ref(1,:))
                    loc1 = node_filtered(ref(1,1)+1).locSeg(j,:);
                end
            end
            if (~isempty(node1(current_Nodes(i,1)).connectionOut))
                for j = 1:numel(node_filtered(current_Nodes(i,1)).connectionOut)/3
                    next_Segments(cntr,:) = node_filtered(current_Nodes(i,1)).connectionOut(j,:);
                    %loc3 is the location of the vessel that is connected
                    %to the current node in the vessel (segment) index
                    loc3 = node_filtered(current_Nodes(i,1)).locSeg(j,:);
                    
                    %locationTracker stores the location of the
                    %inflowing edge considered to the node, so that we
                    %can track cumlative effects experienced by each
                    %edge
                    locationTracker(cntr) = i;
                    
                    CumFluxScaleNew(cntr) = CumFluxScaleInitial(locationTracker(cntr))*node_filtered(current_Nodes(i,1)).Scale(j,:);
                    CumPermeabilityNew(cntr) = CumPermeabilityInitial(locationTracker(cntr))*(1-seg_filtered(loc3).MP);
                    CumDelayNew(cntr) = CumDelayInitial(locationTracker(cntr))+seg_filtered(loc3).Delay;
                    
                    
                    tempRow = CumFluxScaleNew(cntr).*cfunction((t-CumDelayNew(cntr)),C_Func_select)*CumPermeabilityNew(cntr);
                    %This iterative addition means that we superimpose all
                    %the paths leading to the vessel 
                    Concentration_matrix_temp(loc3,:) = Concentration_matrix_temp(loc3,:) + tempRow;
                    
                    %check if the cumulative scale of the path is smaller
                    %than the path threshold, and in that case that path is
                    %ignored
                    if CumFluxScaleInitial(i) < ignore_Threshold
                        next_Segments(cntr,:) = [];
                        CumFluxScaleNew(cntr) = [];
                        CumDelayNew(cntr) = [];
                        cntr = cntr - 1;
                        cntr2 = cntr2 - 1;
                    end
                    cntr = cntr + 1;
                    cntr2 = cntr2 + 1;
                end
            end
        end
        
        cntr = 1;
        
        %Once the propagating front is completely empty, which means its
        %done with propagation, then the breakCheck is changed so that on
        %the next iteration it exits the loop
        if isempty(next_Segments)
            BreakCheck = 0;
        end
       
        %track propagating front size
        PF_Size(cntr3,1) = length(next_Segments);
        
        %Display code progress
        if mod(cntr3,100) == 0
            disp(PF_Size(cntr3,1));
            progress = (3*k/numel(Inlet_Index))*100;
            disp(strcat('Progress:',num2str(progress),'%%'));
        end
        
        
        %Updating the values for the next iteration 
        CumPermeabilityInitial = CumPermeabilityNew;
        CumFluxScaleInitial = CumFluxScaleNew;
        CumDelayInitial = CumDelayNew;
        curSeg = next_Segments;
        next_Segments =[];
        CumFluxScaleNew = [];
        CumDelayNew = [];
        CumDelay = [];
        locationTracker = 1;
        
        P_Counter = P_Counter + numel(curSeg)/3;
        cntr3 = cntr3 +1;
    end
%after each inlet iteration, the temp concentration matrix is superimposed
%unto the total concentration matrix 
    Concentration_matrix = Concentration_matrix + Concentration_matrix_temp;
end
toc
Prop_result = sum(Concentration_matrix,2);
save prop_result.mat Prop_result
toc

%Taking snapshots of the data at desired intervals for visualisation

timePoint = 1:2:100;

for k = 1:numel(timePoint)
    % write to data & to seg_filtered the solution
    data(numParam+1).Name = strcat('Cumulative_Dilution_',num2str(t(timePoint(k))),'s');
    data(numParam+1).Marker = strcat('@',num2str(numParam+1));
    data(numParam+1).Dim=1;
    zero_Tracker  = 1;
    inlet_Tracker = 1;
    pos_Tracker = 0;
    missing_seg_tracker = 0;
    disp(k)
    tic
    %adding the numbers in the correct format
    for i= 1:length(seg_Index1)
        curSeg = seg_Index1(i,:);
        numPoints = seg1(i).NumPoints;
        locSeg = findRowInIndexV(curSeg,seg_filtered_Index,3);
        if isempty(locSeg)
            for j = pos_Tracker+1:pos_Tracker+numPoints
                data(numParam+1).Val(j,1) = 0;
            end
            pos_Tracker = pos_Tracker+numPoints;
        else
            for j = pos_Tracker+1:pos_Tracker+numPoints
                data(numParam+1).Val(j,1) = Concentration_matrix(locSeg,timePoint(k));
            end
            pos_Tracker = pos_Tracker+numPoints;
        end
    end
    data(numParam+1).NumEl = numel(data(numParam+1).Val);
    vdefinition(numParam+1).vname = strcat('Cumulative_Dilution_',num2str(t(timePoint(k))),'s');
    vdefinition(numParam+1).vtype = 'POINT';
    vdefinition(numParam+1).vdatatype = 'float';
    vdefinition(numParam+1).vdatadim = 1;
    vdefinition(numParam+1).vmarker = strcat('@',num2str(numParam+1));
    numParam = numParam+1;
    toc
end
disp('Data write-in complete');