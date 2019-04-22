function [ dataFiltered ] = comAvgRef( data )
%comAvgRef Takes input data (channels = rows) and takes mean of each data
%point (column) and subtracts this from each row. 
%   Inputs:
%       data : data arranged with rows as the individual channels
%
%   Outputs
%       dataFiltered : common average of each column subtracted from every
%       channel.

[ rowNum, columnNum ] = size(data);

cOffset = mean( data, 2);
data = data - cOffset;
cAvg = mean( data );
dataFiltered = NaN( rowNum, columnNum );

for i = 1:rowNum
    dataFiltered( i,: ) = data( i,: ) - cAvg;
end

end

