function [ fileName1 ] = write2Amira( data,vdefinition,definition, fileName );
%write2Amira Writes data to AMira file format
%   INPUTS: DATA structure containing fields for amira
%   definition structure containing number of nodes, edges and edge points
%   vdefinition structure containing variable definitons 

% WRITE_SPATIALGRAPH_AM.m 
% This program writes Amira Graph files 
clc

% Either Import the data straight from the workspace, or save the variables
% by using the following line  'save Amira_Data.mat data definition vdefinition -v7.3;'
% from the command line and and then uncomment the next line to load it. 


%load('Amira_Data.mat')

fid = fopen( fileName, 'w+' );
% Writing the Header
Str_Opener = '# AmiraMesh 3D ASCII 2.0';
fprintf(fid, '%s', Str_Opener);

% Writing the defintion section

def1 = strcat('define',{' '}, definition(1).type,{' '}, definition(1).len);
fprintf(fid, '\n\n\n%s', char(def1));

def1 = strcat('define',{' '}, definition(2).type,{' '}, definition(2).len);
fprintf(fid, '\n%s', char(def1));

def1 = strcat('define',{' '}, definition(3).type,{' '}, definition(3).len);
fprintf(fid, '\n%s', char(def1));

% Writing the parameter section 

fprintf(fid, '\n\n%s', 'Parameters {');
fprintf(fid, '\n\t%s', 'ContentType "HxSpatialGraph"');
fprintf(fid, '\n%s\n', '}');

% Writing the variable definition section 
for i = 1:numel(data)
     if data(i).Dim == 1
            vdef = strcat(vdefinition(i).vtype,{' {'},{' '},vdefinition(i).vdatatype,{' '},vdefinition(i).vname,{' } '},vdefinition(i).vmarker);
            fprintf(fid, '\n%s', char(vdef));
     else
            vdef = strcat(vdefinition(i).vtype,{' {'},{' '},vdefinition(i).vdatatype,{'['},num2str(vdefinition(i).vdatadim),{'] '},vdefinition(i).vname,{' } '},vdefinition(i).vmarker);
            fprintf(fid, '\n%s', char(vdef));
    end
end
disp('Variable Definition Section: Complete')
% Writing end of header line

fprintf(fid, '\n\n%s', '# Data section follows');

% Writing data section

for i = 1:numel(data)
    tic
    fprintf(fid, '\n%s', char(data(i).Marker));
    switch i
        case {1,4}
            for j = 1:(data(i).NumEl)
                d1 = strcat(num2str(data(i).Val(j,1)),{' '},num2str(data(i).Val(j,2)),{' '},num2str(data(i).Val(j,3)));
                fprintf(fid, '\n%s', char(d1));
            end
            fprintf(fid, '\n%s','');
        case 2
            for j = 1:str2double(char(definition(i).len))
                d1 = strcat(num2str(data(i).Val(j,1)),{' '},num2str(data(i).Val(j,2)));
                fprintf(fid, '\n%s', char(d1));
            end
            fprintf(fid, '\n%s','');
        otherwise
            for j = 1:(data(i).NumEl)
                d1 = num2str(data(i).Val(j));
                fprintf(fid, '\n%s', char(d1));
            end
            fprintf(fid, '\n%s','');                 
    end
disp(strcat('Data Field',num2str(i),': complete'))
toc
end
disp('File Writing Complete')

fileName1 = fileName



end

