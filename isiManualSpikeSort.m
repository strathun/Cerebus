function [spikesNew] = isiManualSpikeSort(continuousData,spikes, Fs)
%UNTITLED3 Summary of this function goes here
%   Inputs :
%           
%           spikes : should be index values

%% move to after rejection from other function!!!

spikesNew = spikes;
spikeTime = ( spikes / Fs ) * (1e3);
isiArray = diff(spikeTime);

uncertainEvents = find( isiArray < 1 );
nextSpike = uncertainEvents + 1;
spikeIndex = spikes(uncertainEvents);
nextspikeIndex = spikes(nextSpike);

events = length(uncertainEvents);
x = (( 1:1:(Fs*(2e-3)) ) / Fs ) * 1e3;

for i = 1:events
    subplot(2,1,1)
    be = spikeIndex(i) - ( ( 1e-3 ) * Fs );
    en = spikeIndex(i) + ( ( 1e-3 ) * Fs );
    plot( continuousData(be:en))
    
    subplot(2,1,2)
    be = nextspikeIndex(i) - ( ( 1e-3 ) * Fs );
    en = nextspikeIndex(i) + ( ( 1e-3 ) * Fs );
    plot( continuousData(be:en))
    
    prompt = 'Keep top (1) or bottom (2) or neither (0)?' ;
    g = input(prompt);
    tempIndex = uncertainEvents(i);
    
    if g == 0
        spikesNew(tempIndex) = NaN;
        spikesNew(tempIndex+1) = NaN;
    elseif g == 1
        spikesNew(tempIndex+1) = NaN;
    elseif g == 2
        spikesNew(tempIndex) = NaN;
    end
    
    close
end

end

