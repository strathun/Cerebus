function [ dataMeanSub ] = meanSubtraction( data )
%Takes the mean of each row and subtracts it from itself to remove offsets.
%   Inputs:
%        data : data file to undergo subtraction. Individual channels
%               should be in rows.
%   Outputs: 
%        dataMeanSub : filtered data

[ numRows , ~ ] = size( data );

data1 = mean( data,2 );

for i = 1:numRows
    dataMeanSub(i,:) = data( i,: ) - data1( i );
end

end

