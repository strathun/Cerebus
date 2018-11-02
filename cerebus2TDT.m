function [orderedData] = cerebus2TDT(cerebusData,dataType)
%cerebus2TDT This function will make a new rearranged array that accounts for
%the differences in channel naming conventions between the cerebus and TDT
% [orderedData] = cerebus2TDT(cerebusData,dataType);
%   Inputs : 
%       cerebusData : data should be organized with rows = channel numbers
%       dataType : 0 for array; 1 for cell array
%   Outputs : 
%       orderedData : cerebusData reorderd so that it matches TDT
%       convention.
%
% Channels
% Cerebus : 1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16
% TDT     : 1   9   2   10  3   11  4   12  5   13   6  14   7  15   8  16
%%%%

tdtConverter = [1 9 2 10 3 11 4 12 5 13 6 14 7 15 8 16];

[row, column] = size(cerebusData);

if dataType == 0
    orderedData = NaN(row, column);
    
    for i = 1:length(tdtConverter)
        orderedData(tdtConverter(i),:) = cerebusData(i,:);
    end
end

if dataType == 1
    orderedData = cell( row, 1 );
    
    for i = 1:length( tdtConverter )
        orderedData{i} = cerebusData(tdtConverter(i),:);
    end
    
end
end

