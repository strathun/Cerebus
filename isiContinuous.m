function [isiArray,xAxis] = isiContinuous(spikeTimes,binsize,totalTime)
%   [isiArray,xAxis] = isiContinuous(spikeTimes,binsize,totalTime)
%   Computes average isi for specified intervals of continuous data.
%   Inputs :
%       spikeTimes : contains spike times (in ms)
%       binsize : ms / bin
%       totalTime : time in seconds

cap = totalTime * 1e3 ;
binNum = cap/binsize ;
timeInc = binsize/2 ; 
xAxis = timeInc:binsize:cap;
i = 1;

while ( i * binsize ) <= cap
    ranger{i} = spikeTimes( ( spikeTimes >= ( (i-1) * binsize ) ) & ( spikeTimes <= ( (i) * binsize ) ) );
    i = i + 1;
end

for n = 1:binNum
    temp = diff(ranger{n});
    isiArray(n) = mean(temp);
end

isiArray(isnan(isiArray))=0;

end

