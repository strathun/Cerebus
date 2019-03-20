%% Pulls most relevant data from Cerebus files.
% Use this to pull commonly used data chunks from .ns5 and .nev files
% Run in the directory that contains the .nev and .ns5 file

fileNameGEN = 'TDT15_Iso_1_stim10s_noVitals_2';       % filename (no extension) of datafiles you want to examine

% Appends for other file types
fileNameNS5 = [fileNameGEN '.ns5'];    %
fileNameNEV = [fileNameGEN '.nev'];    %

%%%%%%%%%_Paramaters_%%%%%%%%%%%%%
Fs = 30e3;              % Sampling Frequency, 30e3 for both Cerebus and Ripple
threshold = -4.0 ;      % vrms multiplier for spike detection
fullData = 1;           % analyze full dataset? yes 1, no 0.
startTime = 20;         % if fullData = 0, specificy analysis start and stop time
stopTime = 24;
channels = 16;          % Number of signal recording channels
rejectMod = 2;          % reject waveforms > (rejectMod x mean waveform) at any location
filterType = 1;         % 1 = low pass -> high pass filter
passBandF = 250;        % Frequency in Hz of high passband
passBandFL = 7500;
order = 3;              % Order for filter
chSelect = 13;           % Channel to analyze
ARP = .001;             % Absolute refractory period. Any detected spikes with isi < ARP will be rejected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% raw data

[ NSxFileArray, NSxbasicHeader, NSxchannelHeader, NSxTimeStamps ] = NSxGetMemMapFile( fileNameNS5 );
time = double( NSxTimeStamps ) / 30000;
rawdata = NSxFileArray.Data.NSxData;

V = double( rawdata );

% Cerebus uses .25 uV per bit. Data comes in as bits. Here we're leaving
% units as uV
V = ( V  )/4;

Vfiltered = comAvgRef( V(1:16,:) ); %17th channel is foot switch, Cerebus
% Vfiltered = V(1:16,:);
% Vfiltered = V;

%% Waveform

location = pwd;
[status, waveformdata ] = WFPlot(fileNameGEN,location,[],-4);
[spikeEvents, stimulus, units] = nev2MatSpikesOnly( fileNameNEV );

%%  Prep data

if filterType == 0
    if fullData == 0
        Vfiltered = (Vfiltered(:,(startTime*Fs):(stopTime*Fs)));
    end
    Vfilteredtest = Vfiltered;
    Vfiltered = bandpass(Vfilteredtest(chSelect,:),[750 4000],30e3);
    % Vfilered = Vfiltered(:, ( Fs*.05 ) : ( end - (Fs*.05)) );    %Cuts out filter artefacts at beginning and end of filtered data.

elseif filterType == 1
    %Butterworth filter _ high pass
    [ B, A ] = butter( order, passBandF / ( Fs/2 ), 'high');
    if fullData == 0
        dataMeanSub = Vfiltered(:,(startTime*Fs):(stopTime*Fs)).';    %filter takes column as channels, not rows
    else
        dataMeanSub = Vfiltered.';
    end
    dataHighPass = filter( B, A, dataMeanSub);

    %Butterworth filter _ low pass
    [ B, A ] = butter( order, passBandFL / ( Fs/2 ));
    dataHighPass = filter( B, A, dataHighPass);
    dataHighPass = dataHighPass.';
    Vfiltered = dataHighPass(chSelect,:);
end

%%

x = 1:1:length(Vfiltered);
time = (1:1:length(Vfiltered))/Fs;
timefilt = time;
eventsNew = [];

for i = 1:1
    %grabbing data
    [spikesIndex, threshVal] = spike_detection(Vfiltered(i,:),threshold,1,0);
    [waveforms{i}, timeWave, spikesIndex] = waveformGrabber(Vfiltered(i,:), spikesIndex, 1.6, Fs);
    
%     spikesTemp = cell2mat(spikesIndex);
    spikesTemp = (spikesIndex);
    spikesTime = spikesTemp/Fs;
    spikes{i} = spikesTime;
    
    %insert chopper here
    spikesUse = spikes{1,i}.';
    [choppedData, fat] = spikeChopper(Vfiltered(i,:),spikesUse,Fs,'Threshold',1.6);
    
    Vrms(i) = rms(choppedData);
    
    figure
    plot(timefilt,Vfiltered(i,:))
    hold on
    
    %adds threshold line to plot
    threshLine = ones(1,length(Vfiltered)) * threshVal;
%     plot(timefilt, threshLine,'--')
    title([ 'Background Noise Vrms : ' num2str(Vrms(i)) ' uV'])
    
    %plots reset events onto graph
    for iii = 1:length(eventsNew)
        xVar = eventsNew(iii)/Fs;
        line([xVar (xVar +(xVar/1e5))],[-80 80],'Color','black')
    end
    %xlim([0 (timefilt(end))])
    xlabel('Time (s)');
    
    [events, ~] = size(waveforms{i}) ; 
    figure
    waveform = waveforms{i} ;
    [waveform, spikeEventsNew] = templateMatcher(waveform,rejectMod, spikesIndex,  ARP, Fs); %removes "bad" spikes
    [eventsMod, ~] = size(waveform) ; 
    [threshCount(i+1), ~] = size(waveform); %+1 to add zero at beginning for below...
    spikeEventsRaster{i} = spikeEventsNew;
    
    %SNR calculation (Per RC Kelly (2007) J Neurosci 27:261)
    waveformMean = mean(waveform);
    waveformNoise = waveform - waveformMean;
    waveformNoiseSD(i) = std2(waveformNoise);
    SNRKelly(i) = ( max(waveformMean) - min(waveformMean) ) / ( 2 * waveformNoiseSD(i) );

    %SNR calculation (Per K Ludwig (2009) J Neurophys 0022-3077)
    PkPkNoiseFloor = 6 * ( std( choppedData ) );
    SNRLudwig(i) = ( max(waveformMean) - min(waveformMean) ) / PkPkNoiseFloor;
    
    for ii = 1:eventsMod
        plot(timeWave*1e3, waveform(ii,:),'Color',[.5 .5 .5], 'LineWidth', 1.2)
        hold on
    end
    threshLine = threshLine(1:(length(timeWave*1e3)));
    plot(timeWave*1e3,threshLine,'--')
    
    title([ 'SNRKelly: ' num2str(SNRKelly(i)) ', SNRLudwig: ' num2str(SNRLudwig(i)) ', SpikeCount: ' num2str(threshCount(i+1))])
    
    meanWave = mean(waveform) ;
    plot(timeWave*1e3, meanWave, 'LineWidth', 3)
    xlabel('Time (ms)');
%     ylim([-80 60])
end

% Raster Plot
figure
LineFormat.LineWidth = 0.8;
[xpoints, ypoints] = plotSpikeRaster(spikeEventsRaster,'PlotType','vertline','LineFormat',LineFormat);
hold on
    for iii = 1:length(eventsNew)
        xVar = eventsNew(iii)/Fs;
        line([xVar (xVar +(xVar/1e5))],[-80 80],'Color','blue')
    end
xlabel('Time (s)');

%% puts Raster data on filtered data graph

xpoints = (xpoints.')/Fs;
ypoints = ypoints.';
threshCount = threshCount * 3;

xpointsSorted = xpoints(threshCount(1) + 1: (threshCount(1)+ threshCount(2)));
ypointsSorted = ypoints(threshCount(1) + 1: (threshCount(1)+ threshCount(2)));
figure(2)
ylim([-150 150])
plot(xpointsSorted,((ypointsSorted)*50)-(150),'k')
% xlim([0 4])

%% Plot Cerebus detected waveforms
x = (1:1:length(waveformdata(1).waveforms(:,1))) / Fs ;
figure
[~, commercialCount] = size(waveformdata(1).waveforms);
for i = 1:commercialCount
plot( x, waveformdata(1).waveforms(:,i),'Color',[.5 .5 .5], 'LineWidth', 1.2)
hold on
end

threshLine = ones(1,length(x)) * -28;
threshLine = threshLine(1:(length(x)));
plot(x,threshLine,'--')

meanWave = mean(waveformdata(1).waveforms,2);
plot(x, meanWave, 'LineWidth', 3)
title('Detected Waveforms (Commercial System)')
ylim([-80 60])
%% Plot of raw data for gut check

figure(300)
title('Raw Channels')
for i = 1:1
    figure
    plot(time(1,:),V(i,:))
    xlabel('Time (s)');
end

%% Plots PSD
avgs = 64;

figure
for i = 1:channels
    [pxx1,f] = psdWalker(Vfiltered,avgs,Fs);
    loglog(f,(pxx1))
    hold on
end
ylabel('Noise Voltage ( nV / \surd Hz )','Interpreter','tex')
xlabel('Frequency (Hz)')

figure
for i = 1:channels
    [pxx1,f] = psdWalker(V(chSelect,(20*Fs):(24*Fs)),avgs,Fs);
    loglog(f,(pxx1))
    hold on
end
ylabel('Noise Voltage ( nV / \surd Hz )','Interpreter','tex')
xlabel('Frequency (Hz)')

%% ISI data

spikeEventsDif = diff((spikeEventsNew/Fs)*(1e3));
% [h, edges] = histcounts(spikeEventsDif,64);
edges = 0:1:200;
figure
histogram(spikeEventsDif,edges)
xlim([0 100])
xlabel('Inter-Spike Interval (ms)')
ylabel('Number of Counts')


[isiArray, isiX] = isiContinuous((spikeEventsNew/Fs)*(1e3),100,4);
figure
plot(isiX/(1e3),isiArray)
xlabel('Time (s)')
ylabel('Inter-Spike Interval (ms)')