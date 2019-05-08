function [ ] = rasterColorCoded( spikesByItem, pictureLabels, psthParams, ISI, channel_i, unit_i, colors )
%RASTER makes a raster plot in the current figure
%    Note: if no spikes on any trial of an image, that image will not appear in legend

preAlign = psthParams.psthPre;
postAlign = psthParams.psthPost;
imDur = psthParams.psthImDur;

xlim([-preAlign,imDur+postAlign]);
hold on;
axis ij

yLevel = 0; % accumulated trial index, sets height in raster
legendHandles = [];
for item_i = 1:length(spikesByItem)  
  yLevelStart = yLevel;
  for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
    yLevel = yLevel + 1;
    trialSpikes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i);
    for spike_i = 1:length(trialSpikes.times)
      plot([trialSpikes.times(spike_i) trialSpikes.times(spike_i)],[yLevel-0.4 yLevel+0.4],'color', colors(mod(item_i,size(colors,1)) + 1,:));% + 1),'-'));
    end
  end
  h = plot([0 0],[yLevelStart yLevel],'color', colors(mod(item_i,size(colors,1)) + 1, :));
  legendHandles = vertcat(legendHandles,h);
  plot([imDur imDur],[yLevelStart yLevel],'color', colors(mod(item_i,size(colors,1)) + 1, :));
end
xlimits = xlim();
ylim([0,yLevel]);
if xlimits(2) < imDur + ISI
  plot([imDur+ISI imDur+ISI],[0 yLevel],'b-');
end
xlabel('Time after stimulus onset (ms)');
ylabel('single trials');
set(gca,'YTick',[]);
legend(legendHandles, pictureLabels, 'Location','northeastoutside');
hold off;
end



