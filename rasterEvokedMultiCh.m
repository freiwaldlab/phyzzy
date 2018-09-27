function [ ] = rasterEvokedMultiCh( spikesByItem, lfpByItem, pictureLabels, preAlign, postAlign, imDur, ISI, lfpPaddedBy, channels, channelNames, colors )
%RASTER makes a raster-evoked potential overlay in the current figure
%    Note: if no spikes on any trial of an image, that image will not appear in legend
%    Note: if first on top == 1, preferred stimulus will appear at top of plot 
if length(channels) > 4
  disp('sorry, I cannot display more than four channels at once, because I will run out of linestyles');
  return
end
normalize = 1;
xlim([-preAlign,imDur+postAlign]); 
hold on;
legendHandles = [];
yMax = 0;
if normalize
  spikeHeight = 0.025;
else
  spikeHeight = 0.075*max(max(squeeze(lfpByItem{1}(1,channels(1),:,lfpPaddedBy+1:end-lfpPaddedBy))));
end
unitColors = {'r-','g-','b-','v-'};
channelLinestyles = {'-','--','.-','.'};
legendHandles = [];
yTickPositions = zeros(length(spikesByItem),1);
yOffset = -1;
for item_i = length(spikesByItem):-1:1
  if normalize
    yOffset = yOffset + 1;
  else
    yOffset = yMax - min(min(min(squeeze(lfpByItem{item_i}(1,:,:,lfpPaddedBy+1:end-lfpPaddedBy)))));
  end
  yTickPositions(item_i) = yOffset;
  for channel_i = 1:length(channels)
    channel_index = channels(channel_i); %allows us to plot only a subset of channels, if we want
    if normalize
      lfpNormMult = max(max(squeeze(lfpByItem{item_i}(1,channel_index,:,lfpPaddedBy+1:end-lfpPaddedBy)))) - min(min(squeeze(lfpByItem{item_i}(1,channel_index,:,lfpPaddedBy+1:end-lfpPaddedBy))));
      lfpNormAdd = mean(mean(squeeze(lfpByItem{item_i}(1,channel_index,:,lfpPaddedBy+1:end-lfpPaddedBy))));
    end
    for trial_i = 1:length(spikesByItem{item_i}{channel_index}{1})
      trialLfp = squeeze(lfpByItem{item_i}(1,channel_index,trial_i,lfpPaddedBy+1:end-lfpPaddedBy))';
      if normalize
        trialLfp = (trialLfp - lfpNormAdd)/lfpNormMult;
      end
      trialLfp = trialLfp + yOffset;
      yMax = max(yMax, max(trialLfp));
      h = plot(-preAlign:imDur+postAlign,trialLfp,'color', colors(mod(item_i,size(colors,1)) + 1,:),channelLinestyles{channel_i});
      if item_i == 1 && trial_i == 1
        legendHandles = vertcat(legendHandles, h);
      end
      for unit_i = 1:length(spikesByItem{1}{channel_i}) - 1
        trialSpikes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i);
        for spike_i = 1:length(trialSpikes.times)
          if ~(trialSpikes.times(spike_i) < (-preAlign+1) || trialSpikes.times(spike_i) > imDur + postAlign)
            plot([trialSpikes.times(spike_i) trialSpikes.times(spike_i)],[trialLfp(round(trialSpikes.times(spike_i)+preAlign))-spikeHeight trialLfp(round(trialSpikes.times(spike_i)+preAlign))+spikeHeight],unitColors{unit_i});
          end
        end
      end
    end
  end
end
plot([0 0],ylim(),'b-');
plot([imDur imDur],ylim(),'b-');
xlimits = xlim();
ylimits = ylim();
if xlimits(2) < imDur + ISI
  plot([imDur+ISI imDur+ISI],[0 ylimits(2)],'b-');
end
xlabel('Time after stimulus onset (ms)');
ylabel('evoked potential (au)');
disp(flip(yTickPositions));
set(gca,'YTick',flip(yTickPositions),'YTicklabel',flip(pictureLabels),...
    'box','off','TickDir','out','FontSize',14,'TickLength',[.012 .012],'LineWidth',1.3);
position = get(gca,'position'); % these two lines, and the lines after 'hold off', are a gross hack to get two legends
ax = gca;
legend(legendHandles, channelNames);
hold off;
xlimits = xlim(); ylimits = ylim();
axes('position',position);
hold on;
plot(1,unitColors{1});
plot(1,unitColors{2});
plot(1,unitColors{3});
set(gca,'visible','off');
xlim(xlimits); ylim(ylimits);
legend({'hash','unit 1','unit2'},'location','southeastoutside');
%plot(1,unitColors{4});
axes(ax);
end



