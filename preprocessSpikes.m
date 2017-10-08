function [ spikesByChannel, taskTriggers, channelUnitNames ] = preprocessSpikes( spikeFilename, params )
%UNTITLED5 Summary of this function goes here
%   params is struct with fields
%   - spikeChannels: same length as LFP channels, in the same order, if analyzing both
%   - cPtCal: conversion from spike sample indices to timestep of decimated LFP
%   returns:
%   - spikesByChannel: nChannels x 1 array of structs with fields:
%     times: nSpikes x 1, spike times in ms since start of recording
%     units: nSpikes x 1, unit number of each spike (0 for unsorted)
%     waveforms: nspikes x nSpikeSamples array of spike waveforms
%   - taskTriggers: nPackets x 1 array of structs; serial-digital IO port log
%     
Output.VERBOSE('loading blackrock event file');
assert(logical(exist(spikeFilename,'file')),'The spike-event file you requested does not exist.');
NEV = openNEV(spikeFilename,'read','nosave','noparse'); %note: add param 'report' for verbose
Output.VERBOSE('parsing serial IO packets');
taskTriggers = NEV.Data.SerialDigitalIO;

%%%%% remove spike data from non-spike channels (e.g. reference electrodes), unsort low quality units, and remove noise units
spikesByChannel = repmat(struct('times',[],'units',[],'waveforms',[]),length(params.spikeChannels),1);
unitNames = {'unsorted', 'unit 1','unit 2','unit 3','unit 4','unit 5'};
channelUnitNames = cell(length(params.spikeChannels),1);
for channel_i = 1:length(params.spikeChannels)
  %change units from sample index to ms; type from int32 to double
  tmp.times = params.cPtCal*double(NEV.Data.Spikes.Timestamps(NEV.Data.Spikes.Electrode == params.spikeChannels(channel_i)));
  tmp.units = NEV.Data.Spikes.Unit(NEV.Data.Spikes.Electrode == params.spikeChannels(channel_i));
  tmp.waveforms = NEV.Data.Spikes.Waveform(NEV.Data.Spikes.Electrode == params.spikeChannels(channel_i),:);
  unitNamesTmp = unitNames(1:length(unique(tmp.units)));
  for discard_i = 1:length(params.unitsToDiscard{channel_i})
    tmp.times = tmp.times(tmp.units ~= params.unitsToDiscard{channel_i}(discard_i));
    tmp.waveforms = tmp.waveforms(tmp.units ~= params.unitsToDiscard{channel_i}(discard_i));
    tmp.units = tmp.units(tmp.units ~= params.unitsToDiscard{channel_i}(discard_i));
  end
  assert(~ismember(0,params.unitsToUnsort{channel_i}),'0 cannot appear in params.unitsToUnsort: cannot unsort unsorted');
  for unsort_i = 1:length(params.unitsToUnsort{channel_i})
    tmp.units(tmp.units == params.unitsToUnsort{channel_i}(unsort_i)) = 0;
  end
  spikesByChannel(channel_i) = tmp;
  if ~isempty(tmp.units)
    channelUnitNames{channel_i} = [unitNamesTmp(setdiff(0:(length(unitNamesTmp)-1),union(params.unitsToUnsort{channel_i},params.unitsToDiscard{channel_i}))+1),{'MUA'}];
  end
  Output.VERBOSE(channelUnitNames{channel_i});
end

if params.shiftSpikeWaveforms
  rawSpikes = spikesByChannel;
  for channel_i = 1:length(params.spikeChannels)
    if length(channelUnitNames{channel_i}) == 2
      continue
    end
    for unit_i = 1:(length(channelUnitNames{channel_i})-2)
      unitWaveforms = spikesByChannel(channel_i).waveforms(spikesByChannel(channel_i).units == unit_i,:);
      meanWaveform = mean(unitWaveforms,1);
      meanWaveform = meanWaveform - mean(meanWaveform,2);
      meanWaveform = repmat(meanWaveform,size(unitWaveforms,1),1);
      unitWaveforms = unitWaveforms - mean(unitWaveforms,2);
      shiftedWaveforms = zeros(size(unitWaveforms));
      shiftQuality = zeros(size(unitWaveforms,1),11);
      shifts = -5:5;
      for shift_i = 1:length(shifts)
        shift = shifts(shift_i);
        shiftQuality(:,shift_i) = sum(meanWaveform(:,6+shift:end-5+shift).*unitWaveforms(:,6-shift:end-5-shift),2);
      end
      [~,bestShifts] = max(shiftQuality,[],2);
      figure();
      disp('opening figure');
      hist(bestShifts - 6);
      for spike_i = 1:size(unitWaveforms,1)
        shift = shifts(bestShifts(spike_i));
        shiftedWaveforms(spike_i,6+shift:end-5+shift) = unitWaveforms(spike_i,6-shift:end-5-shift);
        shiftedWaveforms(spike_i,1:6+shift) = shiftedWaveforms(spike_i,6+shift);
        shiftedWaveforms(spike_i,end-5+shift:end) = shiftedWaveforms(spike_i,end-5+shift);
      end
      spikesByChannel(channel_i).waveforms(spikesByChannel(channel_i).units == unit_i,:) = shiftedWaveforms;
    end
  end
end

colors = ['k','r','c','g','b'];
if params.plotSpikeWaveforms
  endTime = 0;
  for channel_i = 1:length(params.spikeChannels)
    try %defense against unit with no spikes
      endTime = max(endTime, spikesByChannel(channel_i).times(end));
    catch
      continue
    end
  end
  halfTime = endTime/2;
  for channel_i = 1:length(params.spikeChannels)
    tmp = spikesByChannel(channel_i);
    fh = figure();
    numPlotColumns = length(unique(tmp.units)) + 1; %extra is for MUA plot
    %initialize top row 
    for subplot_i = 1:numPlotColumns-1
      subplot(3,numPlotColumns,subplot_i);
      title(sprintf('%s Unit %d',params.channelNames{channel_i},subplot_i-1));
      if subplot_i == 1
        ylabel('voltage (uV)');
      end
      hold on
    end
    subplot(3,numPlotColumns,numPlotColumns);
    title(sprintf('%s MUA',params.channelNames{channel_i}));
    hold on
    %initialize middle row 
    for subplot_i = 1:numPlotColumns-1
      subplot(3,numPlotColumns,numPlotColumns+subplot_i);
      title(sprintf('%s Unit %d early',params.channelNames{channel_i},subplot_i-1));
      if subplot_i == 1
        ylabel('voltage (uV)');
      end
      hold on
    end
    subplot(3,numPlotColumns,2*numPlotColumns);
    title(sprintf('%s MUA early',params.channelNames{channel_i}));
    hold on
    %initialize bottom row 
    for subplot_i = 1:numPlotColumns-1
      subplot(3,numPlotColumns,2*numPlotColumns+subplot_i);
      title(sprintf('%s Unit %d late',params.channelNames{channel_i},subplot_i-1));
      if subplot_i == 1
        ylabel('voltage (uV)');
      end
      xlabel('time (ms)');
      hold on
    end
    subplot(3,numPlotColumns,3*numPlotColumns);
    title(sprintf('%s MUA late',params.channelNames{channel_i}));
    xlabel('time (ms)');
    hold on
    tAxis = params.cPtCal*(1:size(tmp.waveforms,2));
    toSkipByUnit = zeros(numPlotColumns-1,1);
    spikesToPlot = 100;
    for unit_i = 1:length(toSkipByUnit)
      toSkipByUnit(unit_i) = floor(sum(tmp.units == unit_i-1)/spikesToPlot);
    end
    for unit_i = 1:numPlotColumns-1
      unitWaveforms = tmp.waveforms(tmp.units == unit_i-1,:);
      unitTimes = tmp.times(tmp.units == unit_i-1);
      unitWaveformsToPlot = unitWaveforms(1:toSkipByUnit(unit_i):size(unitWaveforms,1),:);
      unitTimesToPlot = unitTimes(1:toSkipByUnit(unit_i):size(unitWaveforms,1));
      midPoint = find(unitTimesToPlot > halfTime, 1 );
      subplot(3,numPlotColumns,unit_i);
      for spike_i = 1:length(unitWaveformsToPlot)
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors(unit_i));
      end
      subplot(3,numPlotColumns,numPlotColumns); %MUA plot
      for spike_i = 1:length(unitWaveformsToPlot)
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors(unit_i));
      end
      subplot(3,numPlotColumns,numPlotColumns+unit_i) %second row of the subplot
      for spike_i = 1:midPoint - 1
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors(unit_i));
      end    
      subplot(3,numPlotColumns,2*numPlotColumns); %MUA plot
      for spike_i = 1:midPoint - 1
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors(unit_i));
      end
      subplot(3,numPlotColumns,2*numPlotColumns+unit_i) %third row of the subplot
      for spike_i = midPoint:length(unitWaveformsToPlot)
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors(unit_i));
      end
      subplot(3,numPlotColumns,3*numPlotColumns); %MUA plot
      for spike_i = midPoint:length(unitWaveformsToPlot)
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors(unit_i));
      end
    end
    drawnow;
    if params.plotSpikeWaveforms == 1
      close(fh);
    end
  end
end
if params.spikeWaveformPca
  if ~exist('halfTime','var')
    endTime = 0;
    for channel_i = 1:length(params.spikeChannels)
      endTime = max(endTime, spikesByChannel(channel_i).times(end));
    end
    halfTime = endTime/2;
  end
  for channel_i = 1:length(params.spikeChannels)
    tmp = spikesByChannel(channel_i);
    if size(tmp.waveforms,1) < 3
      continue
    end
    [~,score] = pca(tmp.waveforms,'NumComponents',3);
    numUnits = length(unique(tmp.units));
    fh = figure();
    subplot(3,3,1);
    title(sprintf('%s',params.channelNames{channel_i}));
    xlabel('1st PC coefficient');
    ylabel('2nd PC coefficient');
    hold on
    subplot(3,3,4);
    title(sprintf('%s early',params.channelNames{channel_i}));
    xlabel('1st PC coefficient');
    ylabel('2nd PC coefficient');
    hold on
    subplot(3,3,7);
    title(sprintf('%s late',params.channelNames{channel_i}));
    xlabel('1st PC coefficient');
    ylabel('2nd PC coefficient');
    hold on
    % 1 vs 3
    subplot(3,3,2);
    title(sprintf('%s',params.channelNames{channel_i}));
    xlabel('1st PC coefficient');
    ylabel('3rd PC coefficient');
    hold on
    subplot(3,3,5);
    title(sprintf('%s early',params.channelNames{channel_i}));
    xlabel('1st PC coefficient');
    ylabel('3rd PC coefficient');
    hold on
    subplot(3,3,8);
    title(sprintf('%s late',params.channelNames{channel_i}));
    xlabel('1st PC coefficient');
    ylabel('3rd PC coefficient');
    hold on
    % 2 vs 3
    subplot(3,3,3);
    title(sprintf('%s',params.channelNames{channel_i}));
    xlabel('2nd PC coefficient');
    ylabel('3rd PC coefficient');
    hold on
    subplot(3,3,6);
    title(sprintf('%s early',params.channelNames{channel_i}));
    xlabel('2nd PC coefficient');
    ylabel('3rd PC coefficient');
    hold on
    subplot(3,3,9);
    title(sprintf('%s late',params.channelNames{channel_i}));
    xlabel('2nd PC coefficient');
    ylabel('3rd PC coefficient');
    hold on
    %%% now, draw scatters
    scatterHandles = gobjects(9,1);
    % 1 vs 2
    h = subplot(3,3,1);
    scatterHandles(1) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1,1),score(tmp.units == unit_i-1,2),36,colors(unit_i));
    end
    h = subplot(3,3,4);
    scatterHandles(4) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times < halfTime,1),score(tmp.units == unit_i-1 & tmp.times < halfTime,2),36,colors(unit_i));
    end
    h = subplot(3,3,7);
    scatterHandles(7) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times >= halfTime,1),score(tmp.units == unit_i-1 & tmp.times >= halfTime,2),36,colors(unit_i));
    end
    % 1 vs 3
    h = subplot(3,3,2);
    scatterHandles(2) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1,1),score(tmp.units == unit_i-1,3),36,colors(unit_i));
    end
    h = subplot(3,3,5);
    scatterHandles(5) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times < halfTime,1),score(tmp.units == unit_i-1 & tmp.times < halfTime,3),36,colors(unit_i));
    end
    h = subplot(3,3,8);
    scatterHandles(8) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times >= halfTime,1),score(tmp.units == unit_i-1 & tmp.times >= halfTime,3),36,colors(unit_i));
    end
    % 2 vs 3
    h = subplot(3,3,3);
    scatterHandles(3) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1,2),score(tmp.units == unit_i-1,3),36,colors(unit_i));
    end
    h = subplot(3,3,6);
    scatterHandles(6) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times < halfTime,2),score(tmp.units == unit_i-1 & tmp.times < halfTime,3),36,colors(unit_i));
    end
    h = subplot(3,3,9);
    scatterHandles(9) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times >= halfTime,2),score(tmp.units == unit_i-1 & tmp.times >= halfTime,3),36,colors(unit_i));
    end
    linkaxes(scatterHandles);
    %
    drawnow;
    if params.spikeWaveformPca == 1
      close(fh);
    end
  end
end
if params.shiftSpikeWaveforms
  spikesByChannel = rawSpikes;
end
clear NEV
end

