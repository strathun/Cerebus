% Plots peristimulus histogram 
%NOTE: load spiking data and rawdata and time before running

% Channels
% Cerebus : 1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16
% TDT     : 1   9   2   10  3   11  4   12  5   13   6  14   7  15   8  16


%% sepearates spike times into individual channels
%searches for event for given channel
ch1 = spikes(:,1) == 1;
ch2 = spikes(:,1) == 2;
ch3 = spikes(:,1) == 3;
ch4 = spikes(:,1) == 4;
ch5 = spikes(:,1) == 5;
ch6 = spikes(:,1) == 6;
ch7 = spikes(:,1) == 7;
ch8 = spikes(:,1) == 8;
ch9 = spikes(:,1) == 9;
ch10 = spikes(:,1) == 10;
ch11 = spikes(:,1) == 11;
ch12 = spikes(:,1) == 12;
ch13 = spikes(:,1) == 13;
ch14 = spikes(:,1) == 14;
ch15 = spikes(:,1) == 15;
ch16 = spikes(:,1) == 16;

%%
% spiking = {};
%pulls spike times of given channel to new array spikes#
spiking{1} = spikes(ch1,2);
spiking{2} = spikes(ch2,2);
spiking{3} = spikes(ch3,2);
spiking{4} = spikes(ch4,2);
spiking{5} = spikes(ch5,2);
spiking{6} = spikes(ch6,2);
spiking{7} = spikes(ch7,2);
spiking{8} = spikes(ch8,2);
spiking{9} = spikes(ch9,2);
spiking{10} = spikes(ch10,2);
spiking{11} = spikes(ch11,2);
spiking{12} = spikes(ch12,2);
spiking{13} = spikes(ch13,2);
spiking{14} = spikes(ch14,2);
spiking{15} = spikes(ch15,2);
spiking{16} = spikes(ch16,2);

spiking = spiking.';

%%
figure
plotSpikeRaster(spiking,'PlotType','vertline')
% xlim([60 180])    %used to view high activity period of waking files

%% load stim interval on top of raster. Must load rawdata and time first
hold on
yyaxis right
area(time,rawdata(17,:),'FaceAlpha',.2)
ylim([.2e4 1e4]);
set(gca,'YTickLabel',[]);
%%
channel = 14; %channel used to build histogram
time = 180; %seconds
window = 1; %desired window length in seconds
bins = time/window;    % total number of bins
h = histcounts(spiking{channel},bins);
% 
% figure
% histogram(spiking{channel},bins);

figure(2)
plot(h,'LineWidth',1.2)
hold on
