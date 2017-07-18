function [ spikesByChannel, taskTriggers ] = preprocessSpikes( spikeFilename, params )
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
NEV = openNEV(spikeFilename,'read','nosave','noparse'); %note: add param 'report' for verbose
Output.VERBOSE('parsing serial IO packets');
taskTriggers = NEV.Data.SerialDigitalIO;

%%%%% remove spike data from non-spike channels (e.g. reference electrodes), unsort low quality units, and remove noise units
spikesByChannel = repmat(struct('times',[],'units',[],'waveforms',[]),length(params.spikeChannels),1);
for channel_i = 1:length(params.spikeChannels)
  %change units from sample index to ms; type from int32 to double
  tmp.times = params.cPtCal*double(NEV.Data.Spikes.Timestamps(NEV.Data.Spikes.Electrode == params.spikeChannels(channel_i)));
  tmp.units = NEV.Data.Spikes.Unit(NEV.Data.Spikes.Electrode == params.spikeChannels(channel_i));
  tmp.waveforms = NEV.Data.Spikes.Waveform(NEV.Data.Spikes.Electrode == params.spikeChannels(channel_i),:);
  for discard_i = 1:length(params.unitsToDiscard{channel_i})
    tmp.times = tmp.times(tmp.units ~= params.unitsToDiscard{channel_i}(discard_i));
    tmp.waveforms = tmp.waveforms(tmp.units ~= params.unitsToDiscard{channel_i}(discard_i));
    tmp.units = tmp.units(tmp.units ~= params.unitsToDiscard{channel_i}(discard_i));
  end
  for unsort_i = 1:length(params.unitsToUnsort{channel_i})
    tmp.units(tmp.units == params.unitsToUnsort{channel_i}(unsort_i)) = 0;
  end
  spikesByChannel(channel_i) = tmp;
end


colors = ['k','r','c','g','b'];
if params.plotSpikeWaveforms
  endTime = 0;
  for channel_i = 1:length(params.spikeChannels)
    endTime = max(endTime, spikesByChannel(channel_i).times(end));
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
    skippedByUnit = zeros(numPlotColumns-1,1); %numPlotColumns-1 = numUnits
    toSkipByUnit = zeros(numPlotColumns-1,1);
    spikesToPlot = 100;
    for unit_i = 1:length(toSkipByUnit)
      toSkipByUnit(unit_i) = floor(sum(tmp.units == unit_i)/spikesToPlot);
    end
    for spike_i = 1:length(tmp.waveforms)
      skippedByUnit(tmp.units(spike_i)+1) = skippedByUnit(tmp.units(spike_i)+1) + 1;
      if ~(skippedByUnit(tmp.units(spike_i)+1) > toSkipByUnit(tmp.units(spike_i)+1))
        continue
      else
        skippedByUnit(tmp.units(spike_i)+1) = 0;
      end
      subplot(3,numPlotColumns,double(tmp.units(spike_i))+1);
      plot(tAxis,tmp.waveforms(spike_i,:),'color',colors(tmp.units(spike_i)+1));
      subplot(3,numPlotColumns,numPlotColumns); %MUA plot
      plot(tAxis,tmp.waveforms(spike_i,:),'color',colors(tmp.units(spike_i)+1));
      if tmp.times(spike_i) < halfTime
        subplot(3,numPlotColumns,numPlotColumns+double(tmp.units(spike_i))+1) %second row of the subplot
        plot(tAxis,tmp.waveforms(spike_i,:),'color',colors(tmp.units(spike_i)+1));
        subplot(3,numPlotColumns,2*numPlotColumns); %MUA plot
        plot(tAxis,tmp.waveforms(spike_i,:),'color',colors(tmp.units(spike_i)+1));
      else
        subplot(3,numPlotColumns,2*numPlotColumns+double(tmp.units(spike_i))+1) %third row of the subplot
        plot(tAxis,tmp.waveforms(spike_i,:),'color',colors(tmp.units(spike_i)+1));
        subplot(3,numPlotColumns,3*numPlotColumns); %MUA plot
        plot(tAxis,tmp.waveforms(spike_i,:),'color',colors(tmp.units(spike_i)+1));
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
  legendStrings = {'unsorted','unit 1','unit 2','unit 3','unit 4','unit 5'};
  for channel_i = 1:length(params.spikeChannels)
    tmp = spikesByChannel(channel_i);
    disp('starting pca');
    [coeff,score] = pca(tmp.waveforms,'NumComponents',2);
    disp('finished pca');
    numUnits = length(unique(tmp.units));
    fh = figure();
    subplot(1,3,1);
    title(sprintf('%s',params.channelNames{channel_i}));
    xlabel('1st PCA coefficient');
    ylabel('second PCA coefficient');
    hold on
    subplot(1,3,2);
    title(sprintf('%s early',params.channelNames{channel_i}));
    xlabel('1st PCA coefficient');
    ylabel('second PCA coefficient');
    hold on
    subplot(1,3,3);
    title(sprintf('%s late',params.channelNames{channel_i}));
    xlabel('1st PCA coefficient');
    ylabel('second PCA coefficient');
    hold on
    for unit_i = 1:numUnits
%       disp(size(coeff));
%       disp(unit_i);
%       disp(size(colors));
%       disp(size(tmp.units));
      scatter(score(tmp.units == unit_i & tmp.times < halfTime,1),score(tmp.units == unit_i & tmp.times < halfTime,2),36,colors(unit_i));
    end
    subplot(1,3,1);
    for unit_i = 1:numUnits
      scatter(coeff(tmp.units == unit_i & tmp.times >= halfTime,1),coeff(tmp.units == unit_i & tmp.times >= halfTime,2),36,colors(unit_i));
    end
    legend(legendStrings(1:numUnits));
    if spikeWaveformPca == 1
      close(fh);
    end
  end
end
clear NEV
end

