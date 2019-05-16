function [ ] = raster( spikesByItem, pictureLabels, psthParams, ISI, channel_i, unit_i, colors )
%RASTER makes a raster plot in the current figure
%    Note: if no spikes on any trial of an image, that image will not appear in legend

preAlign = psthParams.psthPre;
postAlign = psthParams.psthPost;
imDur = psthParams.psthImDur;

xlim([-preAlign,imDur+postAlign]); 
hold on;
axis ij
yLevel = -0.5; % accumulated trial index, sets height in raster
legendHandles = [];

for item_i = 1:length(spikesByItem)  
  yLevelStart = yLevel;
  for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
    yLevel = yLevel + 1;
    trialSpikes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i);
    for spike_i = 1:length(trialSpikes.times)
      plot([trialSpikes.times(spike_i) trialSpikes.times(spike_i)],[yLevel-0.4 yLevel+0.4],'color', colors(mod(item_i,size(colors,1)) + 1, :));
    end
  end
  %Plot relevant stimulus onset and offset marker
  h = plot([0 0],[yLevelStart+0.5 yLevel+0.5],'color', colors(mod(item_i,size(colors,1)) + 1, :));
  legendHandles = vertcat(legendHandles,h);
  plot([imDur imDur],[yLevelStart+0.5 yLevel+0.5],'color', colors(mod(item_i,size(colors,1)) + 1, :));
  if item_i ~= length(spikesByItem)
    plot([xlim()],[yLevel+0.5 yLevel+0.5],'color','black', 'LineWidth', 2, 'LineStyle', '--')
  end
end
xlimits = xlim();
ylim([0,yLevel+0.5]);
if xlimits(2) < imDur + ISI
  plot([imDur+ISI imDur+ISI],[0 yLevel],'b-');
end
xlabel('Time after stimulus onset (ms)');
ylabel('single trials');
set(gca,'YTick',[]);
set(gca,'YTickLabels',[]);
legend(legendHandles, pictureLabels, 'Location','northeastoutside');
hold off;
end



