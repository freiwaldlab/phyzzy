function [ spikesByChannel, taskTriggers, channelUnitNames ] = preprocessSpikes(spikeFilename, params )
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
NEV = openNEV(spikeFilename,'nosave','nomat'); %note: add param 'report' for verbose

%Check to see if new openNEV (where Spikes struct has been flipped), and if
%so, flip it back.
if isfield(NEV.Data.Spikes,'WaveformUnit')
  NEV.Data.Spikes.TimeStamp = NEV.Data.Spikes.TimeStamp';
  NEV.Data.Spikes.Electrode = NEV.Data.Spikes.Electrode';
  NEV.Data.Spikes.Unit = NEV.Data.Spikes.Unit';
  NEV.Data.Spikes.Waveform = NEV.Data.Spikes.Waveform';
  NEV.Data.SerialDigitalIO.TimeStamp = NEV.Data.SerialDigitalIO.TimeStamp';
  NEV.Data.SerialDigitalIO.TimeStampSec = NEV.Data.SerialDigitalIO.TimeStampSec';
end

Output.VERBOSE('parsing serial IO packets');
taskTriggers = NEV.Data.SerialDigitalIO;
if ~params.needSpikes
  spikesByChannel = {};
  channelUnitNames = {};
  return
end

%%%%% remove spike data from non-spike channels (e.g. reference electrodes), unsort low quality units, and remove noise units
spikesByChannel = repmat(struct('times',[],'units',[],'waveforms',[]),length(params.spikeChannels),1);
unitNames = {'unsorted', 'unit 1','unit 2','unit 3','unit 4','unit 5', 'unit 6','unit 7','unit 8','unit 9','unit 10'};
channelUnitNames = cell(length(params.spikeChannels),1);

[spikeFilePath, spikeFile, ~] = fileparts(spikeFilename);

%Checks for resorted spikes, overwrites NEV structure with new unit
%assignments and time stamps.
if isfield(params,'offlineSorted') && params.offlineSorted == 1
    spikeFilenameSorted = [spikeFilePath '/' spikeFile '.xls'];
    assert(logical(exist(spikeFilenameSorted,'file')),'The Offline sorted spike file you requested does not exist.');
    spikeMat = xlsread(spikeFilenameSorted);
    %Overwrite NEV fields
    NEV.Data.Spikes.Electrode = spikeMat(:,1);
    NEV.Data.Spikes.Unit = spikeMat(:,2);
    NEV.Data.Spikes.TimeStamp = spikeMat(:,3)*30e3; %Sampling Freq should likely be a variable pulled from elsewhere.
    NEV.Data.Spikes.Waveform = spikeMat(:,4:end);
end

if isfield(params, 'waveClus') && params.waveClus
  try
    %Temporarily add directories needed for wave_clus
    addpath(genpath('dependencies/wave_clus'))
  
    %use the typical naming convention to find the contious trace (ns5)
    lfpFilename = [spikeFilePath '/' spikeFile '.ns5'];
    
    %parse the ns5, or see if they are already parsed.
    parsedData = parse_data_NSx(lfpFilename, 2); %(filename,max_memo_GB)
    
    %Find spikes in the files which matter
    for ii = 1:length(parsedData)
      tmpFilename = parsedData{ii};
      tmpCut = split(tmpFilename, ["_","."]);
      if str2double(tmpCut{end-1}(3:end)) < 129 %Connector Banks on Blackrock are channels 1 - 128.
       output_paths = Get_spikes(parsedData{ii});
      else
        delete(parsedData{ii})
      end
    end
    
    %Cluster them, based on either a specified param file or the default.
    if isfield(params, 'paramHandle')
      clusterResults = Do_clustering(output_paths, 'par', params.paramHandle);
    else
      clusterResults = Do_clustering(output_paths);
    end
    
    tmp = load(output_paths{1},'par');
    threshold = mean(tmp.par.threshold);
    if strcmp(tmp.par.detection,'neg')
      threshold = -1*threshold;
    end
    clear tmp
    
    %Cycle through cluster results (done per electrode) and load them into
    %a temporary NEV structure.
    [tmpSpikes.TimeStamp, tmpSpikes.Electrode, tmpSpikes.Unit, tmpSpikes.Waveform] = deal([]);
    for ii=1:length(clusterResults)
      WC = load(clusterResults{ii}); %Electrode should actually parse the name of the file. 
      tmpSpikes.Unit = vertcat(tmpSpikes.Unit, WC.cluster_class(:,1));
      tmpSpikes.TimeStamp = vertcat(tmpSpikes.TimeStamp, WC.cluster_class(:,2));
      tmpSpikes.Waveform = vertcat(tmpSpikes.Waveform, WC.spikes);
      tmpSpikes.Electrode = vertcat(tmpSpikes.Electrode, ones(length(WC.cluster_class), 1)*ii);
    end
    
    %TimeStamp are already in ms, so unscale them so later code works.
    tmpSpikes.TimeStamp = tmpSpikes.TimeStamp/params.cPtCal;
    
    %Conservative measure - if the mean waveform is too close to the
    %threshold, merge into MUA. First, find the mean waveform of each
    %cluster class.
    if params.waveClusReclass
      clusters = unique(tmpSpikes.Unit);
      mean_wave = nan(length(clusters)-1,size(tmpSpikes.Waveform,2));
      for cluster_i = 2:length(clusters) %start @ 2 to ignore 0th cluster.
        cluster_id = clusters(cluster_i);
        mean_wave(cluster_id,:) = mean(tmpSpikes.Waveform(tmpSpikes.Unit == cluster_id, :));
      end
      
      %Plot Average waveforms, and threshold for detection
      figure('Name','waveClusResult - AverageWaveform','Visible','On','NumberTitle','off');
      title('Average Waveforms, thresholds for detection and reclustering')
      hold on
      for wave_i = 1:size(mean_wave,1)
        plot(mean_wave(wave_i,:),'LineWidth',3)
      end
      plot([0 length(mean_wave)], [threshold threshold],'Linewidth',3,'color','r')
      
      %Reassignment the clusters within a certain fraction of the threshold
      %back to MUA (cluster 0).
      waveform_trough = min(mean_wave, [], 2);
      MUA_threshold = threshold * params.waveClusMUAThreshold;
      plot([0 length(mean_wave)], [MUA_threshold MUA_threshold],'Linewidth',3,'color','b')
      clusters_to_MUA = find(waveform_trough > MUA_threshold); %Cluster numbers that need to be 0 now.
      for ii = 1:length(clusters_to_MUA)
        tmpSpikes.Unit(tmpSpikes.Unit == clusters_to_MUA(ii)) = 0;
      end
      WC.par.unsortedClusters = clusters_to_MUA;
      
      %Now re-assign clusters as to not skip numbers.
      [sd,r] = sort(unique(tmpSpikes.Unit),'ascend');
      new_clusters = r - 1;
      for ii = 1:length(unique(tmpSpikes.Unit))
        tmpSpikes.Unit(tmpSpikes.Unit == sd(ii)) = new_clusters(ii);
      end
    end
    
    % the high dimensional feature space used to cluster can be visualized,
    % assuming the switch is selected and there are at least 2 clusters to
    % see. 
    if params.waveClusProjPlot && (length(unique(tmpSpikes.Unit)) > 2)
      PlotAllFeatures(WC);
    end
    
    %Overwrite the NEV data.
    NEV.Data.Spikes = tmpSpikes; 
    
    %Save figures
    if isfield(params, 'saveFig') && params.saveFig
      figHandles = findobj('Type', 'figure');
      for ii = 1:length(figHandles)
        if strncmp(figHandles(ii).Name, 'waveClus', 8)
          savefig(figHandles(ii), [params.outDir spikeFile '_' figHandles(ii).Name], 'compact') %Will save files
        end
      end
    end

    %Append waveClus params to the AnalysisParams file in the outDir.
    waveClusParams = WC.par;
    save([params.outDir 'AnalysisParams.mat'], 'waveClusParams', '-append');
  catch
    warning('waveClus failure - proceeding unsorted')
    error('waveClus failure')
  end
    %Clean up - Remove added paths, delete folder with files if requested.
    rmpath(genpath('dependencies/wave_clus'))
    switch params.waveClusClear
      case 1
        delete([spikeFilePath '/' spikeFile '_parsed' '/times*'])
        delete([spikeFilePath '/' spikeFile '_parsed' '/*dg*'])
        delete([spikeFilePath '/' spikeFile '_parsed' '/*.txt'])
        delete([spikeFilePath '/' spikeFile '_parsed' '/*.png'])
        delete([spikeFilePath '/' spikeFile '_parsed' '/*spikes.mat'])
      case 2
        rmdir([spikeFilePath '/' spikeFile '_parsed'], 's');
    end
end

for channel_i = 1:length(params.spikeChannels)
  %change units from sample index to ms; type from int32 to double
  tmp.times = params.cPtCal*double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode == params.spikeChannels(channel_i)));
  tmp.units = NEV.Data.Spikes.Unit(NEV.Data.Spikes.Electrode == params.spikeChannels(channel_i));
  tmp.waveforms = NEV.Data.Spikes.Waveform(NEV.Data.Spikes.Electrode == params.spikeChannels(channel_i),:);
  if min(tmp.units) > 0 && isempty(params.unitsToUnsort{channel_i})
    unitNamesTmp = unitNames(2:length(unique(tmp.units))+1);
  else
    unitNamesTmp = unitNames(1:length(unique(tmp.units)));
  end
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
      title('Waveform shift (in samples)')
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

if isfield(params, 'spikeWaveformsColors')
  if isnumeric(params.spikeWaveformsColors) %assumes a RGB array, where each row is a color
    for ii = 1:size(params.spikeWaveformsColors,1)
      colors{ii} = params.spikeWaveformsColors(ii,:);
    end
  else
      colors = params.spikeWaveformsColors;
  end
else
  colors = {'k','r','c','g','b'};
end

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
    fh = figure('Name','Spike Waveform Development','NumberTitle','off');
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
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors{mod(unit_i-1,length(colors))+1});
      end
      subplot(3,numPlotColumns,numPlotColumns); %MUA plot
      for spike_i = 1:length(unitWaveformsToPlot)
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors{mod(unit_i-1,length(colors))+1});
      end
      subplot(3,numPlotColumns,numPlotColumns+unit_i) %second row of the subplot
      for spike_i = 1:midPoint - 1
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors{mod(unit_i-1,length(colors))+1});
      end    
      subplot(3,numPlotColumns,2*numPlotColumns); %MUA plot
      for spike_i = 1:midPoint - 1
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors{mod(unit_i-1,length(colors))+1});
      end
      subplot(3,numPlotColumns,2*numPlotColumns+unit_i) %third row of the subplot
      for spike_i = midPoint:length(unitWaveformsToPlot)
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors{mod(unit_i-1,length(colors))+1});
      end
      subplot(3,numPlotColumns,3*numPlotColumns); %MUA plot
      for spike_i = midPoint:length(unitWaveformsToPlot)
        plot(tAxis,unitWaveformsToPlot(spike_i,:),'color',colors{mod(unit_i-1,length(colors))+1});
      end
    end
    drawnow;
    if isfield(params, 'saveFig') && params.saveFig
      figHandles = findobj('Type', 'figure');
      savefig(figHandles(1), [params.outDir spikeFile '_SpikeWaveforms'], 'compact')
    end
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
      scatter(score(tmp.units == unit_i-1,1),score(tmp.units == unit_i-1,2),36,colors{mod(unit_i-1,length(colors))+1});
    end
    h = subplot(3,3,4);
    scatterHandles(4) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times < halfTime,1),score(tmp.units == unit_i-1 & tmp.times < halfTime,2),36,colors{mod(unit_i-1,length(colors))+1});
    end
    h = subplot(3,3,7);
    scatterHandles(7) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times >= halfTime,1),score(tmp.units == unit_i-1 & tmp.times >= halfTime,2),36,colors{mod(unit_i-1,length(colors))+1});
    end
    % 1 vs 3
    h = subplot(3,3,2);
    scatterHandles(2) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1,1),score(tmp.units == unit_i-1,3),36,colors{mod(unit_i-1,length(colors))+1});
    end
    h = subplot(3,3,5);
    scatterHandles(5) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times < halfTime,1),score(tmp.units == unit_i-1 & tmp.times < halfTime,3),36,colors{mod(unit_i-1,length(colors))+1});
    end
    h = subplot(3,3,8);
    scatterHandles(8) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times >= halfTime,1),score(tmp.units == unit_i-1 & tmp.times >= halfTime,3),36,colors{mod(unit_i-1,length(colors))+1});
    end
    % 2 vs 3
    h = subplot(3,3,3);
    scatterHandles(3) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1,2),score(tmp.units == unit_i-1,3),36,colors{mod(unit_i-1,length(colors))+1});
    end
    h = subplot(3,3,6);
    scatterHandles(6) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times < halfTime,2),score(tmp.units == unit_i-1 & tmp.times < halfTime,3),36,colors{mod(unit_i-1,length(colors))+1});
    end
    h = subplot(3,3,9);
    scatterHandles(9) = h;
    for unit_i = 1:numUnits
      scatter(score(tmp.units == unit_i-1 & tmp.times >= halfTime,2),score(tmp.units == unit_i-1 & tmp.times >= halfTime,3),36,colors{mod(unit_i-1,length(colors))+1});
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