function [ ] = rasterColorCoded(spikesByItem, pictureLabels, psthParams, ISI, channel_i, unit_i, attendedObjData)
%RASTER makes a raster plot in the current figure of the specified unit and
%channel data. Uses information from the attendedObjData structure to color
%code the spikes based on what object the subject was looking at. 
%    Note: if no spikes on any trial of an image, that image will not appear in legend

preAlign = psthParams.psthPre;
postAlign = psthParams.psthPost;
imDur = psthParams.psthImDur;

xlim([-preAlign,imDur+postAlign]);
hold on;
axis ij

%Reshape attendedObjData to be correctly indexed.
attendObjInd = zeros(length(pictureLabels),1);
for obj_ind = 1:length(pictureLabels)
  attendObjInd(obj_ind) = find(strcmp(attendedObjData.eventIDs, pictureLabels{obj_ind}));
end

attendedObjData.eventIDs = attendedObjData.eventIDs(attendObjInd);
attendedObjData.attendedObjVect = attendedObjData.attendedObjVect(attendObjInd);
attendedObjData.frameStartInd = attendedObjData.frameStartInd(attendObjInd);
attendedObjData.frameEndInd = attendedObjData.frameEndInd(attendObjInd);
%attendedObjData.colorCode{end+1} = [0, 0, 0]; %Black for non-stim spikes.
uniqueColorTrials = zeros(length(attendedObjData.colorCode),1);
uniqueColorInds = zeros(length(attendedObjData.colorCode),3);

%Create the color index for all the spikes
for item_i = 1:length(spikesByItem)
  for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
    %Initialize a vector for colors
    spikeColorTmp = zeros(length(spikesByItem{item_i}{channel_i}{unit_i}(trial_i).times),3);
    spikeTimes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i).times;
    spikeTimes = spikeTimes - psthParams.latency; %Shift by latency.

    %Assign remanining spikes colors
    for spike_i = 1:length(spikeTimes)
      if ~(spikeTimes(spike_i) < 0 || spikeTimes(spike_i) > imDur)     %Pre and Post stimuli spikes are black.
        objInd = find(spikeTimes(spike_i)-attendedObjData.frameStartInd{item_i}>0,1,'last');
        objAtt = attendedObjData.attendedObjVect{item_i}(trial_i, objInd);
        colorInd = strcmp(attendedObjData.objList, objAtt);
        spikeColorTmp(spike_i,:) = attendedObjData.colorCode{colorInd};
        if sum(logical(uniqueColorTrials)) ~= sum(logical(uniqueColorTrials+colorInd)) %If adding the index makes a new logical vecotr, this is a new color.
          uniqueColorInds(colorInd,:) = [item_i trial_i spike_i]; %Used later to get handles for legend.
          uniqueColorTrials = uniqueColorTrials + colorInd;
        end
      end
    end
    spikeColor{item_i}{trial_i}.color = spikeColorTmp;
  end
end

%Plot spikes
yLevel = 0; % accumulated trial index, sets height in raster
legendHandles = [];
spikeHandles = [];
spikeLabels = [];
for item_i = 1:length(spikesByItem)
  yLevelStart = yLevel;
  for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
    yLevel = yLevel + 1;
    trialSpikes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i);
    trialColors = spikeColor{item_i}{trial_i}.color;
    for spike_i = 1:length(trialSpikes.times)
      if any(sum(uniqueColorInds == [item_i trial_i spike_i], 2) == 3)
        handleInd = sum(uniqueColorInds == [item_i trial_i spike_i], 2) == 3;
        spikeLabels = [spikeLabels; attendedObjData.objList(handleInd)];
        h = plot([trialSpikes.times(spike_i) trialSpikes.times(spike_i)],[yLevel-0.4 yLevel+0.4],'color', trialColors(spike_i,:),'LineWidth',2);
        spikeHandles = vertcat(spikeHandles, h);
      else
        plot([trialSpikes.times(spike_i) trialSpikes.times(spike_i)],[yLevel-0.4 yLevel+0.4],'color', trialColors(spike_i,:),'LineWidth',2);
      end
    end
  end
  %Plot Stim start, stim end, and Bottom dividing bar
  h = plot([0 0],[yLevelStart yLevel+0.5],'color', 'black', 'LineWidth', 2); % Stimulus Onset
  plot([imDur imDur],[yLevelStart yLevel+0.5],'color', 'black', 'LineWidth', 2); % Stimulus Offset
  plot([xlim()],[yLevel+0.5 yLevel+0.5],'color','black', 'LineWidth', 2, 'LineStyle', '--'); %Stimuli dividing line
  
  legendHandles = vertcat(legendHandles,h);
end
%Add Spike handles to legend;
legendHandles = vertcat(legendHandles, spikeHandles);
legendLabels = vertcat(pictureLabels, spikeLabels);
xlimits = xlim();
ylim([0,yLevel]);
if xlimits(2) < imDur + ISI
  plot([imDur+ISI imDur+ISI],[0 yLevel],'b-');
end
xlabel('Time after stimulus onset (ms)');
ylabel('single trials');
set(gca,'YTick',[]);
legend(legendHandles, legendLabels, 'Location','northeastoutside');
hold off;
end



