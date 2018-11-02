function [ combedData ] = voltageGrabber( data, highLow )
%voltageGrabber Allows the user to selectively pull data from an array
%using the activity of a reference channel (must be placed in final row of
%array) as the pointer. 
%
%   Inputs:
%       data : this is your data array with each row assigned to a channel.
%           The final row must be the reference channel (high/low channel).
%           SHOULD USE RAWWW DATA!
%       highLow : set to 0 to pull data when switch or stimulus channel is
%           low and 1 for when it is active (high).
%
%   Outputs: 
%       combedData : will return the original data (arranged as before)
%           with only the data points sampled when the stimulus channel was
%           on/off.
%       

data = double(data);
a = size(data);
x = max(data(a(1),:))*.6;

eventArray = zeros(1,a(2));

%combs through data looking for sample numbers of on/off foot pedal events.
for i = 1:a(2)
    if data(a(1),i) > x
        eventArray(i) = 1;
    end
end

if highLow == 0 
    eventArray = abs(eventArray - 1);
end

data(data==0) = 1e40;

data1 = zeros( a(1), a(2) );

for ii = 1:a(1)-1
    data1(ii,:) = data(ii,:).*eventArray;
end

rows = a(1)-1;
columns = sum(eventArray(:)==1);

data1(data1==0) = [];
data1(data1==1e40) = 0;

combedData = reshape(data1,rows,columns);

end