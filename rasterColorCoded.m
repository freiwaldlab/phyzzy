function [ ] = rasterColorCoded(figHandle, spikesByItem, pictureLabels, psthParams, ISI, channel_i, unit_i, attendedObjData)
%RASTER makes a raster plot in the current figure of the specified unit and
%channel data. Uses information from the attendedObjData structure to color
%code the spikes based on what object the subject was looking at. 
%    Note: if no spikes on any trial of an image, that image will not appear in legend

%Color Spikes vs Color Shaded Region
colorType = 2; %1 = Spikes, 2 = Shaded Regions

%Set the figure handle as the current handle
set(0, 'CurrentFigure', figHandle)
%Attach resizing function
%figHandle.ResizeFcn = @reAlignColorLegend;

%Unpack Variables
preAlign = psthParams.psthPre;
postAlign = psthParams.psthPost;
imDur = psthParams.psthImDur;
rasterAxes = axes();

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
attendedObjData.tracePlotData = attendedObjData.tracePlotData(attendObjInd);

%attendedObjData.colorCode{end+1} = [0, 0, 0]; %Black for non-stim spikes.
uniqueColorTrials = zeros(length(attendedObjData.colorCode),1);
uniqueColorInds = zeros(length(attendedObjData.colorCode),3);

if colorType == 1
  %Create the color index for all the spikes
  for item_i = 1:length(spikesByItem)
    for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
      %Initialize a vector for colors
      spikeTimes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i).times;
      spikeColorTmp = zeros(length(spikeTimes),3);
      spikeTimes = spikeTimes - psthParams.latency; %Shift by latency.
      
      %Assign remanining spikes colors
      for spike_i = 1:length(spikeTimes)
        if ~(spikeTimes(spike_i) < 0 || spikeTimes(spike_i) > imDur)     %Pre and Post stimuli spikes are black.
          frameInd = find((spikeTimes(spike_i)-attendedObjData.frameStartInd{item_i}) > 0,1,'last');
          objAtt = attendedObjData.attendedObjVect{item_i}(trial_i, frameInd);
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
end

%Plot the Shaded area, if appropriate
if colorType == 2
  attObjData = cat(1, attendedObjData.tracePlotData{1:end});
  im = image(attObjData);
  im.XData(2) = psthParams.psthImDur;
  im.AlphaData = 0.5;
  im.YData = [im.YData(1)-0.5 im.YData(2)-0.5];
  figHandle.UserData.shadedAreaHandle = im;
end

%Plot spikes
yLevel = -0.5; % accumulated trial index, sets height in raster
legendHandles = [];
% spikeHandles = [];
trialLabels = [];

for item_i = 1:length(spikesByItem)
  yLevelStart = yLevel;
  for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
    yLevel = yLevel + 1;
    trialSpikes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i);
    trialLabels = [trialLabels 1:length(spikesByItem{item_i}{channel_i}{unit_i})];
    for spike_i = 1:length(trialSpikes.times)
      if colorType == 1
        trialColors = spikeColor{item_i}{trial_i}.color;
        plot([trialSpikes.times(spike_i) trialSpikes.times(spike_i)],[yLevel-0.4 yLevel+0.4],'color', trialColors(spike_i,:));
      else
        plot([trialSpikes.times(spike_i) trialSpikes.times(spike_i)],[yLevel-0.4 yLevel+0.4],'color', 'black');
      end
    end
  end
  %Plot Stim start, stim end, and Bottom dividing bar
  h = plot([0 0],[yLevelStart+0.5 yLevel+0.5],'color', 'black', 'LineWidth', 2); % Stimulus Onset
  plot([imDur imDur],[yLevelStart+0.5 yLevel+0.5],'color', 'black', 'LineWidth', 2); % Stimulus Offset
  if item_i ~= length(spikesByItem)
    plot([xlim()],[yLevel+0.5 yLevel+0.5],'color','black', 'LineWidth', 2, 'LineStyle', '--'); %Stimuli dividing line
  end
  legendHandles = vertcat(legendHandles,h);
end

title(figHandle.Name)
xlimits = xlim();
ylim([0,yLevel+0.5]);
if xlimits(2) < imDur + ISI
  plot([imDur+ISI imDur+ISI],[0 yLevel],'b-');
end
trialLabelsTmp = cell(1, length(trialLabels));
for tick_i = 1:length(trialLabels)
  trialLabelsTmp{tick_i} = num2str(trialLabels(tick_i));
end

set(gca,'YTickLabels',trialLabelsTmp);
set(gca,'TickLength',[0 0]);
yticks((1:length(trialLabelsTmp))-0.5)
xlabel('Time after stimulus onset (ms)');
ylabel('single trials');
legend(legendHandles, pictureLabels, 'Location','northeastoutside');
hold off;

%Create a legend which gives the color code for the shaded region or the
%individual spikes. Invisible axes with bar plots for each color.
invisAx = axes('Color','none','XColor','none');

handleToThisBarSeries = gobjects(length(attendedObjData.objList),1);
for b = 1 : length(attendedObjData.objList)
	handleToThisBarSeries(b) = bar(nan, nan, 'BarWidth', 0.5);
	set(handleToThisBarSeries(b), 'FaceColor', attendedObjData.colorCode{b});
	hold on;
end
legend(attendedObjData.objList, 'Location', 'southeastoutside');
invisAx.Visible = 'off';
linkprop([rasterAxes invisAx],'Position');

end

