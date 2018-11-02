%% Spike code to find the mid point of spikes
% Created by Aaron Choong and Jesslyn Sutanti
% Modified 10/8/2017
% Supervised by Kostas, Yan & Masoud
% Last modified: 05/18/2018 by Hunter Strathman (University of Utah)
%   Added ability to set threshold input as voltage value or rms multiplier
%   (voltORrms).

function [spikes_index, thresholdVal] = spike_detection(spike_data,threshold, voltORrms, crossORmax)
%% This spike code finds the mid point of the spikes instead of the peaks
% Faster computation than spike_times (~5x faster for big data)

%-----------------------------%
% INPUT VARIABLES DEFINITIONS %
%-----------------------------%

% 1) spike_data holds the raw neural data with voltage values (the index will be its time value)
% 2) threshold is the voltage level where the spikes are thresholded
%       against. Can be positive or negative.
% 3) voltORrms: 0 means threshold will be interpreted as a 
        % voltage, 1 means as a factor multipled by the RMS. (RMS will be
        % calculated by this function.

%------------------------------%
% OUTPUT VARIABLES DEFINITIONS %
%------------------------------%

% 1) spikes is an array of all the time stamps of the middle of the spikes

%% Initialisation of variables
if voltORrms == 1
    Vrms = rms( spike_data );
    threshold_temp = Vrms * threshold;
    threshold = threshold_temp ;
    thresholdVal = threshold;
end

if threshold > 0 

spike_index = find((spike_data)>threshold);    % Temporary variable to hold trace of the spikes (all values above the spike threshold is recorded)
spike_index_array = zeros(1,length(spike_index));    % To hold all the midpoint of spikes
i=1;                                    % Index of spike_index matrix


% Start of spike sorting 
while ( i < length(spike_index))
    n=0; % To hold how much values it should skip
    temp_var = 0;
    while((spike_index(i) == spike_index(i+n)-n) && i+n<length(spike_index)) % This loop checks if the next element in the array corresponds to the same spike                                             
        if(temp_var<spike_data(spike_index(i+n)))
            temp_var = spike_data(spike_index(i+n));    
            spike_index_array(1,i) = spike_index(i+n);  % Holds the peak of the spike
        end
        n=n+1; % Increments if the next raw neural data corresponds to the same spike   
    end

    i=i+n;                              % Increment i by the skip value
end                                     % End of while-loop

% Allocating output matrix 
spikes_index = spike_index_array(spike_index_array~=0);   % Allocate the output as a non-zero array 

else
spike_index = find((spike_data)<threshold);    % Temporary variable to hold trace of the spikes (all values above the spike threshold is recorded)
spike_index_array = zeros(1,length(spike_index));    % To hold all the midpoint of spikes
i=1;                                    % Index of spike_index matrix


% Start of spike sorting 
while ( i < length(spike_index))
    n=0; % To hold how much values it should skip _corrected_hs
    temp_var = 0;
    while((spike_index(i) == spike_index(i+n)-n) && i+n<length(spike_index)) % This loop checks if the next element in the array corresponds to the same spike                                  
        if(temp_var > spike_data(spike_index(i+n)))
            temp_var = spike_data(spike_index(i+n)); 
            if n > 0 
                if crossORmax == 0
                    n = n+1;
                    continue
                elseif crossORmax == 1
                    spike_index_array(1,i) = spike_index(i+n);  % Holds the peak of the spike
                end
            else
                spike_index_array(1,i) = spike_index(i+n);  % Holds the peak of the spike
            end
        end 
        n=n+1; % Increments if the next raw neural data corresponds to the same spike           
    end 

    i=i+n;                              % Increment i by the skip value
end                                     % End of while-loop

% Allocating output matrix 
spikes_index = spike_index_array(spike_index_array~=0);   % Allocate the output as a non-zero array 
end

end 