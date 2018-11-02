function [neuron1,neuron2] = neuronSorter(data, threshold, dataType)
%neuronSorter takes an array of isolated spiking events and separates traces 
% of different sizes (possible different neurons). UNFINISHED
%   Inputs : 
%       data : should be organized with rows = single firing event. Specify double vs cell
%       threshold : value used to separate spikes
%       dataType : 1 = double; 2 = cell

[~, events] = size(data);
ii = 1;
iii = 1;

if dataType == 2
    data = cell2array(data);
end

for i = 1:events
    if min(data(:,i)) > threshold
        neuron1(ii,:) = data(:,i);
        ii = ii + 1;
    else
        neuron2(iii,:) = data(:,i);
        iii = iii + 1;
    end
end

end

