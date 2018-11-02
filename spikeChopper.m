function [choppedData, fat] = spikeChopper(rawdata,eventarray,Fs,eventLocation,duration)
%UNTITLED Summary of this function goes here
%   Inputs :
%       rawdata: row of raw voltage to be chopped up.
%       eventarray : array containing spike events (in time). Arrange as column.
%       Fs : sampling frequency
%       eventLocation : 'Beginning' or 'Threshold'. Beginning indicates that the
%           timing in the eventarray is for the beginning of the event,
%           Threshold indicates it is at either the peak or trough.
%       duration : duration of the event of interest. (in ms).
%   Outputs :
%       choppedData : original data with sections removed
%       fat : chopped out bits

[ events, ~ ] = size( eventarray ); 
endofData = length(rawdata);
choppedData = rawdata;
fat = [];
eventDur = ( duration / (1e3) ) * Fs ; 

for i = 1:events
    index = floor( eventarray(i) * Fs );
    temp = [];
    if index-eventDur >= 1 && index+eventDur <= endofData
        if eventLocation == 'Beginning'
            temp = rawdata( index :(index+eventDur) );
            choppedData( index :(index+eventDur) ) = NaN;
        elseif eventLocation == 'Threshold'
            temp = rawdata( index - ( eventDur / 2 ) :( index + ( eventDur / 2 ) ) );
            choppedData( index - ( eventDur / 2 ) :( index + ( eventDur / 2 ) ) ) = NaN;
        end
        fat = [fat temp];
    end
end

choppedData = choppedData(~isnan(choppedData));

end



