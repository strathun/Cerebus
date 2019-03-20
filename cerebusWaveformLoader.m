%cerebusWaveformLoader

fileNameGen = 'datafile1138';   %set to desired file. Must be in directory
Fs = 30e3;
[timeWaveform, ~] = size(waveformdata(1).waveforms);
timeWaveform = ( 1:1:timeWaveform ) / Fs ;

[status, waveformdata ] = WFPlot(fileNameGEN,location,[],4);
[channels, ~ ] = size(waveformdata);

for i = 1:channels
    [~, waveforms ] = size(waveformdata(i).waveforms);
    figure
    for ii = 1:waveforms
        plot(timeWaveform, waveformdata(i).waveforms(:,ii))
        hold on
    end
    meanWaveform = mean(waveformdata(i).waveforms,2)
    plot(timeWaveform, meanWaveform, 'LineWidth', 3)
end
