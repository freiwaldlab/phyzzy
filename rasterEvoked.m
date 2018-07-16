function [ ] = rasterEvoked( spikesByItem, lfpByItem, pictureLabels, preAlign, postAlign, imDur, ISI, lfpPaddedBy, channel_i, colors, firstOnTop )
%RASTER makes a raster-evoked potential overlay in the current figure
%    Note: if no spikes on any trial of an image, that image will not appear in legend
%    Note: if first on top == 1, first stimulus will appear at top of plot 
xlim([-preAlign,imDur+postAlign]); 
hold on;
legendHandles = [];
yMax = 0;
spikeHeight = 0.075*max(max(squeeze(lfpByItem{1}(1,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy))));
unitColors = {'r-','g-','b-','v-'}; %todo: fix hardcode!!!!!!!
if ~firstOnTop
  for item_i = 1:length(spikesByItem)  
    yOffset = yMax - min(min(squeeze(lfpByItem{item_i}(1,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy))));
    for trial_i = 1:length(spikesByItem{item_i}{channel_i}{1})
      trialLfp = squeeze(lfpByItem{item_i}(1,channel_i,trial_i,lfpPaddedBy+1:end-lfpPaddedBy))' + yOffset;
      yMax = max(yMax, max(trialLfp));
      h = plot(-preAlign:imDur+postAlign,trialLfp, strcat(colors(mod(item_i,length(colors)) + 1),'-'));
      if trial_i == 1
        legendHandles = vertcat(legendHandles,h);
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
else
  for item_i = length(spikesByItem):-1:1  
    yOffset = yMax - min(min(squeeze(lfpByItem{item_i}(1,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy))));
    for trial_i = 1:length(spikesByItem{item_i}{channel_i}{1})
      trialLfp = squeeze(lfpByItem{item_i}(1,channel_i,trial_i,lfpPaddedBy+1:end-lfpPaddedBy))' + yOffset;
      yMax = max(yMax, max(trialLfp));
      h = plot(-preAlign:imDur+postAlign,trialLfp, strcat(colors(mod(item_i,length(colors)) + 1),'-'));
      if trial_i == 1
        legendHandles = vertcat(legendHandles,h);
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
set(gca,'YTick',[]);
position = get(gca,'position'); % these two lines, and the lines after 'hold off', are a gross hack to get two legends
ax = gca;
if firstOnTop
  legend(flipud(legendHandles), pictureLabels, 'Location','northeastoutside');
else
  legend(legendHandles, pictureLabels, 'Location','northeastoutside');
end
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



