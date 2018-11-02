%% Waveform Plot

channel = 2;    %channel of interest
Fs = 30e3;      %sampling Freq

[length, traces] = size(waveformdata(channel).waveforms);
x = ( ( 1:1:length ) / Fs ) * 1e3 ; 

figure
for i = 1:length
    plot( x, waveformdata(channel).waveforms(:,i) )
    hold on
end

meanWave = mean(waveformdata(channel).waveforms,2);
plot( x, meanWave , 'LineWidth', 3)

ylabel('(uV)')
xlabel('Time (ms)')
