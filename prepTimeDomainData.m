% This script will filter our raw time domain data

Fs = 30e3;              % Sampling Frequency
passBandF = 500; % 180; %750;        % Frequency in Hz of high passband
bandstop = 26;
order = 5;              % Order for Butterworth Filter
filterType = 1;         % 0 for lowpass; 1 for high pass butter; 2 for bandstop; 3 to bypass filter 
plotyn = 1 ;            % 0 does not plot ; 1 does

V = double( rawdata );
% V = double( Vbreathing );

% Cerebus uses .25 uV per bit. Data comes in as bits. 
V = ( V * ( 1e-6 ) )/4;

%cutting out bad electrode (13) and foot pedal;
V1 = V( 1:12,: );
V2 = V( 14:16,: );
V = [ V1 ; V2 ];

%all electrodes
%V = V(1:16,:);

%% Filtering
%Removing offset
dataMeanSub = meanSubtraction( V );

if filterType == 1
    %Butterworth filter _ high pass
    [ B, A ] = butter( order, passBandF / ( Fs/2 ), 'high');
    dataMeanSub = dataMeanSub.';    %filter takes column as channels, not rows
    dataHighPass = filter( B, A, dataMeanSub);
    dataHighPass = dataHighPass.';
elseif filterType == 0
    %Butterworth filter _ low pass
    [ B, A ] = butter( order, passBandF / ( Fs/2 ));
    dataMeanSub = dataMeanSub.';    %filter takes column as channels, not rows
    dataHighPass = filter( B, A, dataMeanSub);
    dataHighPass = dataHighPass.';
 elseif filterType == 2
    %Butterworth filter _ bandstop
    d = designfilt('bandstopiir','FilterOrder',6, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);
    dataMeanSub = dataMeanSub.';
    dataHighPass = filtfilt(d,dataMeanSub);    %bad naming convention...this isn't really high pass.
    dataHighPass = dataHighPass.';
else
    dataHighPass = dataMeanSub ;
end

%Common Average Referencing
Vfiltered = comAvgRef( dataHighPass ); 

%% Plot data

if plotyn == 1
    figure
    plot( time, V )

    figure
    plot( time, Vfiltered )
end
