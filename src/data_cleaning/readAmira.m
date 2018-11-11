function [ data, definition, vdefinition ] = readAmira(filename)
%readAmira takes input the string name of the amira txt file in the same
%path of the funtion file and reads it to provide structure with the parsed
%data from the txt file
%   Data: the different data encoded in the file 
%   definition: the definitons of the data types 
%   vdefinition: the different variables that are in the file and what
%   their type is
%   Input: Filename as a string, File must be a valid amira spatialGraph
%   file converted to a .txt file 

%This program reads Amira Graph files


fileID = fopen(filename,'r');
curLine = fgetl(fileID);
mgc = '# AmiraMesh 3D ASCII 2.0';
if strcmp(mgc,curLine) == 0
    error_Message = 'Error, READ_SPAITALGRAPH_AM: Not a valid Amira file';
    error(error_Message)
    fclose(fileID);
    return
end
frewind(fileID)
%Reading the header
curLine = fgetl(fileID);
i = 1;
j = 1;
vdefinition(j).vname = [];
definition(i).type = [];
while ischar(curLine)
    
    %Defintion Line
    if ~isempty(strfind(curLine,'define'));
        spl = strsplit(curLine,' ');
        nspl = numel(spl);
        if nspl ~= 3
            error('Error, READ_SPATIALGRAPH_AM: Unexpected formatting on defintion line');
            return
        else
            definition(i).type = spl(2);
            definition(i).len = spl(3);
            i = i+1;
        end
    end
    
    %Parameter Line
    if ~isempty(strfind(curLine,'Parameters'))
        paramLine = fgetl(fileID);
        pspl = strsplit(paramLine,' ');
        npspl = numel(pspl);
        mgc = '"HxSpatialGraph"';
        if strcmp(pspl(npspl),mgc) == 0
            error('Error, READ_SPATIALGRAPH_AM: Not a SpatialGraph file!')
            return
        end
    end
    
    %Variable Definition line
    if ~isempty(strfind(curLine,'@'))
        vspl = strsplit(curLine,' ');
        nvspl = numel(vspl);
        if nvspl == 6
            vdefinition(j).vtype = vspl{1};
            tmp = strsplit(vspl{3},{'[',']'});
            ntmp = numel(tmp);
            if ntmp == 3
                vdefinition(j).vdatatype = tmp{1};
                vdefinition(j).vdatadim = str2double(tmp{2});
            else
                vdefinition(j).vdatatype = tmp{1};
                vdefinition(j).vdatadim = 1;
            end
            vdefinition(j).vname = vspl{4};
            vdefinition(j).vmarker = vspl{6};
        end
        j= j+1;
    end
    
    curLine = strtrim(fgetl(fileID));
    %End of header line
    if ~isempty(strfind(curLine,'#'))
        break
        curLine = 0;
    end
    
end

for l = 1:j-1
    if isempty(vdefinition(l).vname)
        error('Error, READ_SPATIALGRAPH_AM: No variables defined in header')
        return
    end
end

fseek(fileID, 2, 0);


%Reading in the data section
for i = 1:numel(vdefinition)
    data(i).Name = vdefinition(i).vname; 
    data(i).Marker = vdefinition(i).vmarker;
    data(i).Dim = vdefinition(i).vdatadim;
    switch i
        case {1, 4}
            A = fscanf(fileID,'%f', [3,Inf]);
            data(i).NumEl = numel(A)/3;
        case 2
            A = fscanf(fileID,'%f', [2,Inf]);
            data(i).NumEl = numel(A)/2;
        otherwise
            A = fscanf(fileID,'%f', [1,Inf]);
            data(i).NumEl = numel(A);
    end
    data(i).Val = transpose(A);
    position = ftell(fileID);
    fseek(fileID, 3, 0);
    disp('Done with stage: ')
    disp(i)
end
fclose(fileID);
disp('File Read: Complete')


end

