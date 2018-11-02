% spectrogram settings that look the best so far...
% currently using rawdataWaking03. 

limits = [-142 -138] ;
channel = 4;

figure

%%
% subplot(3,1,2);
subplot(2,1,2)
spectrogram( Vfiltered( channel,: ), 200, 20, 3e3, 30e3, 'yaxis' ); %choose channel here
%colorbar off


% 
% limits = [-120 -80]; %breathing
% caxis(limits);
% xlim([1 3]);
% ylim([.15 3]);
ylim([0 .2])
xlim([2 2.2])

colormap(jet);

%%
subplot(2,1,1)
plot(time, rawdata(17,:))
xlim([120 132])
% subplot(3,1,1);
% plot(time/60,Vfiltered(channel,: ) )
% xlim([1 3]);

