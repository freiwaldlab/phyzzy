function [ spikesByImage, spikesByCategory, lfpByImage, lfpByCategory, categoryList, picFiles ] = multiChSessSerene(  )
%multiChSessSerene is a top level function to analyze units, MUA, lfps, coupling, and RFs. 
%   - handles single-channel and multi-channel sessions
%   - relies on raw visiko (.log) and blackrock (.ns5) files
%   - aligns visiko events to blackrock clock (preprocessLogFile.m)
%   - excludes trials with broken fixation, fix spot flash, or (optionally)
%     juice delivery (via exludeStimuli)
%   - decimates and (optionally) filters raw lfp data (default, to 1 kHz)
%   - makes psth's and evoked potentials for (possibly overlapping)categories, 
%     specified at picParamsFilename
%   - RF analyses: spike rate, evoked power, mean spike latency
%   - performs time-frequency and spike-field coherence analyses for each
%     channel, using the multitaper method implemented in Chronux
%   - performs across channel coupling analyses: spike-spike, field-field,
%     and spike-field, using Chronux
%   - todo: perform band-resolved Granger causality analysis across bands
%   - todo: make tuning curves for parameterized stimuli (e.g. view angle)
%   Notes:
%   - currently, in graphs, MUA is top-index unit; e.g. for unsorted u0 and
%     sorted u1 and u2, MUA (u0+u1+u2 spikes) would be u3
%   - currently, assumes that 'face' and 'nonface' appear as categories in
%     picParamsFilename
%   - chronux params currently set when tf analyses start; ctl-f chr_params
%   Depends:
%   - chronux toolbox
%   - partial contents of 'dependencies' folder (details coming)
%   - R2016a (or later) if joint psth-evoked potential plots desired


addpath('dependencies');
addpath(genpath('chronux_2_12/spectral_analysis'));
addpath(genpath('chronux_2_12/locfit'));
addpath(genpath('blackrock_read'));
addpath(genpath('xml2struct'));

%%%%%%%  USER PARAMETERS %%%%%%%%
needNS2 = 0;
runNum = '004';
dateSubject = '161011ALAN';
volume = '/Volumes/OptiHDD/Freiwald/ALAN_DATA/Blackrock'; %'/Volumes/Users-1/User/Desktop%' %
logVolume = '/Volumes/OptiHDD/Freiwald/ALAN_DATA/Visiko';
picParamsFilename = '/Volumes/OptiHDD/Freiwald/AnalysisSerene/StimParamsFullFOB3.mat';
categoryListSlim = {'humanFace','monkeyFace','place','fruit','humanBody','monkeyBody','techno'}; %minimal cat list for clean plots
save_fig = 1;
verbose = 0;
spikeChannels = [1,33,34];
ns5channels = [1,33,34]; 
channel_names = {'ML','AL','AM'};
ns5channelMicroVoltMap = [8191/32764 8191/32764 8191/32764];  %TODO: not yet implemented: will hold conversion from ns5 raw values
accelChannels = [];
accelChannelNames = {};
common_ref = [0, 35, 35]; %TODO: not implemented; will allow software re-refrence across headstages
cPtCal = 1/30; % conversion from spike sample indices to milliseconds
decimateFactorPass1 = 6; %note: product of the two decimate factors should be 30, if 1 khz samples desired
decimateFactorPass2 = 5;
samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)

% high-pass filter design bool. and params; see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
N = 8;
Fpass = 1;
Astop = 100;
Apass = 0.5;
Fs = 1000;
hp1Hz = designfilt('highpassiir', 'FilterOrder',N,'PassbandFrequency',Fpass, ...
  'StopbandAttenuation',Astop,'PassbandRipple',Apass,'SampleRate',Fs);  % if 0, no filter applied and can ignore 
ephysFilter = ''; % if filtering desired, ephysFilter is a digitalFilter 

% for lfps, constrain first and (optional) last [n m] samples to 0 mean
DCSUB_SAM = [20 20]; % for no DC sub, [0 ...], for both ends, [sam sam], for start only [sam 0] 
% 
psthPre = 0; % if e.g. +200, then start psth 200ms before trial onset; todo: allow for values other than zero in chronux plots
psthImDur = 400;  % todo: get this (and pre/post?) from log file
psthPost = 200;
smoothingWidth = 10; %psth smoothing width, in ms
%
frCalcOn = 80;
frCalcOff = psthImDur+frCalcOn;
% Boolean variables to specify which computations to perform; TODO: read
% from config file, eventually with conditional on log file info
makeImPSTH = 1;
makeCatPSTH = 1;
imageTF = 0;
catTF = 1;
calcLatencyRF = 0;
crossTF = 1;
calcEvokedPowerRF = 0;
calcCoherenceRFcpt = 0;
calcCoherenceRFcc = 0;
calcCoherenceRFptpt = 0;
calcGrangerRF = 0;

%%% set paths and directories %%%
ns2Filename = sprintf('%s/%s/%s%s.ns2',volume,dateSubject,dateSubject,runNum);
ns5Filename = sprintf('%s/%s/%s%s.ns5',volume,dateSubject,dateSubject,runNum);
nevFilename = sprintf('%s/%s/%s%s.nev',volume,dateSubject,dateSubject,runNum);
logFilename = sprintf('%s/%s/%s0%s.log',logVolume,dateSubject,dateSubject,runNum);
%
figures_dir = sprintf('/Volumes/OptiHDD/Freiwald/ALAN_DATA/Analyzed/%s/',dateSubject);
if save_fig && ~exist(figures_dir,'dir')
  mkdir(figures_dir);
end
%
load('cocode2.mat');
mystyle= map;


%%%%%%%%%%%%%%%%%

%%%% first, load visiko data, transform it to blackrock clock reference, and load spike data    %
[ taskData, NEV ] = preprocessLogFile(nevFilename,logFilename);
%todo: parse ns2 data, if needed for accelerometers, eye tracking, etc.
if needNS2
  ns2obj = openNSx(ns2Filename,'report','read');
  ns2header = ns2obj.MetaTags;
  ns2data = ns2obj.Data;
end
% load 30 kSample/sec data
ns5obj = openNSx(ns5Filename,'report','read');
ns5header = ns5obj.MetaTags;
assert(ns5header.SamplingFreq == 30000, 'error: expected ns5 sampling freq 30 ks/sec');
ns5dataRaw = ns5obj.Data; %note: returns each channel as a row

% sort ns5 data so channel indexing matches indexing in ns5channels array
ns5channelMap = zeros(length(ns5channels),1);
for i = 1:length(ns5channels)
  assert(any(ns5header.ChannelID == ns5channels(i)), strcat('error: requested analysis for unrecorded ns5 channel: ',num2str(ns5channels(i))));
  ns5channelMap(i) = find(ns5header.ChannelID == ns5channels(i));
end
ns5channelMap = ns5channelMap(ns5channelMap > 0);
ns5data = ns5dataRaw(ns5channelMap,:);
ns5dataDec = zeros(size(ns5data,1),ceil(size(ns5data,2)/30));

% convert scaled units to microvolts, and decimate (note: decimation broken
% into two steps, per matlab doc recommendation
disp('decimating, scaling, and filtering LFP');
for i = 1:size(ns5data,1)
  ns5dataDec(i,:) = ns5channelMicroVoltMap(i)*decimate(decimate(ns5channelMicroVoltMap(i)*ns5data(i,:),decimateFactorPass1),decimateFactorPass2);
  if isa(ephysFilter,'digitalFilter')
    ns5dataDec(i,:) = filter(FILTER, ns5dataDec(i,:));
  end
end
ns5data = ns5dataDec;
if verbose
  disp('size ns5 data'); disp(size(ns5data));
  disp('size ns5 data after decimation'); disp(size(ns5data));
end
disp('done decimating, scaling, and filtering LFP');
clear ns5DataRaw
if verbose
  disp('ns5 channel map'); disp(ns5channelMap);
  disp('channel id info'); disp(ns5header.ChannelID);
  disp('size ns5data'); disp(size(ns5data));
end


%%%%% exclude stimuli for fixation out, flash on, frame dropped, (accel high, juice on)
taskDataAll = taskData;
% params are ( taskData, fixPre, fixPost, flashPre, flashPost, varargin )
taskData = excludeStimuli( taskData, 200, 200, 200, 200, 'juicePre',500,'juicePost', 0);

%%%%% remove spike data from non-spike channels (e.g. reference electrodes)
spikesByChannel = repmat(struct('times',[],'units',[],'waveforms',[]),length(spikeChannels),1);
for i = 1:length(spikeChannels)
  tmp.times = cPtCal*double(NEV.Data.Spikes.Timestamps(NEV.Data.Spikes.Electrode == spikeChannels(i))); %change units from sample index to ms; type from int32 to double
  tmp.units = NEV.Data.Spikes.Unit(NEV.Data.Spikes.Electrode == spikeChannels(i));
  tmp.waveforms = NEV.Data.Spikes.Waveform(NEV.Data.Spikes.Electrode == spikeChannels(i));
  spikesByChannel(i) = tmp;
end
clear NEV

%%%%% load the stimulus category information   %
tmp = load(picParamsFilename);
picCategories = tmp.paramArray;
categoryList = tmp.categoryLabels;
picFiles = unique(taskData.pictureFilenames);
onsetsByImage = cell(length(picFiles),1);
offsetsByImage = cell(length(picFiles),1);
for i = 1:length(picFiles)
  onsetsByImage{i} = taskData.pictureStartTimes(strcmp(taskData.pictureFilenames,picFiles{i}));
  offsetsByImage{i} = taskData.pictureEndTimes(strcmp(taskData.pictureFilenames,picFiles{i}));
end
if taskData.RFmap
  jumpsByImage = cell(length(picFiles),1);
  for i = 1:length(picFiles)
    jumpsByImage{i} = taskData.pictureJumps(strcmp(taskData.pictureFilenames,picFiles{i}),:);
  end
end

%%%%%  align spikes by trial, and sort by image  %

spikesByImage = cell(length(picFiles),1);
psthEmptyByImage = cell(length(picFiles),1);
for image_i = 1:length(picFiles)
  im_spikesByChannel = cell(length(spikeChannels),1); %channels x 1
  im_emptyByChannel = cell(length(spikeChannels),1);
  onsets = onsetsByImage{image_i}; 
  offsets = offsetsByImage{image_i};
  for channel_i = 1:length(spikeChannels)
    channelSpikes = spikesByChannel(channel_i);
    units = sort(unique(channelSpikes.units));
    imageChannelSpikes = cell(length(units)+1,1);
    imageChannelEmpty = cell(length(units)+1,1);
    for unit_i = 1:length(units)
      tstamps = channelSpikes.times(channelSpikes.units == units(unit_i));
      imageUnitSpikes = repmat(struct('times',[]),length(onsets),1); 
      empty = 1;
      for trial_i = 1:length(onsets)
        imageUnitSpikes(trial_i) = struct('times', tstamps(tstamps>=onsets(trial_i) - psthPre & tstamps <= offsets(trial_i) + psthPost) - onsets(trial_i));
        if ~isempty(imageUnitSpikes(trial_i).times)
          empty = 0;
        end
      end
      imageChannelSpikes{unit_i} = imageUnitSpikes;
      imageChannelEmpty{unit_i} = empty;
    end
    % channel MUA
    tstamps = channelSpikes.times;
    imageChannelMUA = repmat(struct('times',[]),length(onsets),1);
    empty = 1;
    for trial_i = 1:length(onsets)
      imageChannelMUA(trial_i) = struct('times', tstamps(tstamps>=onsets(trial_i) - psthPre & tstamps <= offsets(trial_i) + psthPost) - onsets(trial_i));
      if ~isempty(imageChannelMUA(trial_i).times)
        empty = 0;
      end
    end
    imageChannelSpikes{end} = imageChannelMUA;
    imageChannelEmpty{end} = empty;
    im_spikesByChannel{channel_i} = imageChannelSpikes;
    im_emptyByChannel{channel_i} = imageChannelEmpty;
  end
  spikesByImage{image_i} = im_spikesByChannel;
  psthEmptyByImage{image_i} = im_emptyByChannel;
end

% find which line of picCategories corresponds to each picture
picFilePicCategoryIndex = zeros(size(picFiles));
for image_i = 1:length(picFiles)
  for picCat_i = 1:length(picCategories)
    if ~isempty(regexp(picFiles{image_i},picCategories{picCat_i}{1},'ONCE'))
      picFilePicCategoryIndex(image_i) = picCat_i;
    end
  end
end

assert(all(picFilePicCategoryIndex ~= 0),'error: images missing from picture category directory');

onsetsByCategory = cell(length(categoryList));
offsetsByCategory = cell(length(categoryList));
for cat_i = 1:length(categoryList)
  catOnsets = [];
  catOffsets = [];
  for image_i = 1:length(picFiles)
    if any(strcmp(picCategories{picFilePicCategoryIndex(image_i)},categoryList{cat_i}))
      catOnsets = vertcat(catOnsets,onsetsByImage{image_i});
      catOffsets = vertcat(catOffsets, offsetsByImage{image_i});
    end
  end
  onsetsByCategory{cat_i} = catOnsets;
  offsetsByCategory{cat_i} = catOffsets;
end
if verbose
  disp('numel cat onsets');
  disp(numel(catOnsets));
end

% now, spikes by cateogry
% note: has same structure as for image, so can repmat one of those
spikesByCategory = cell(length(categoryList),1); 
psthEmptyByCategory = cell(length(categoryList),1);
for cat_i = 1:length(categoryList)
  cat_spikesByChannel = cell(length(spikeChannels),1); %channels x 1
  cat_emptyByChannel = cell(length(spikeChannels),1);
  onsets = onsetsByCategory{cat_i}; 
  offsets = offsetsByCategory{cat_i};
  for channel_i = 1:length(spikeChannels)
    channelSpikes = spikesByChannel(channel_i);
    units = sort(unique(channelSpikes.units));
    catChannelSpikes = cell(length(units)+1,1);
    catChannelEmpty = cell(length(units)+1,1);
    for unit_i = 1:length(units)
      tstamps = channelSpikes.times(channelSpikes.units == units(unit_i));
      catUnitSpikes = repmat(struct('times',[]),length(onsets),1); 
      empty = 1;
      for trial_i = 1:length(onsets)
        catUnitSpikes(trial_i) = struct('times', tstamps(tstamps>=onsets(trial_i) - psthPre & tstamps <= offsets(trial_i) + psthPost) - onsets(trial_i));
        if ~isempty(catUnitSpikes(trial_i).times)
          empty = 0;
        end
      end
      catChannelEmpty{unit_i} = empty;
      catChannelSpikes{unit_i} = catUnitSpikes;
    end
    % channel MUA
    tstamps = channelSpikes.times;
    catChannelMUA = repmat(struct('times',[]),length(onsets),1);
    empty = 1;
    for trial_i = 1:length(onsets)
      catChannelMUA(trial_i) = struct('times', tstamps(tstamps>=onsets(trial_i) - psthPre & tstamps <= offsets(trial_i) + psthPost) - onsets(trial_i));
      if ~isempty(catChannelMUA(trial_i))
        empty = 0;
      end
    end
    catChannelSpikes{end} = catChannelMUA;
    catChannelEmpty{end} = empty;
    cat_spikesByChannel{channel_i} = catChannelSpikes;
    cat_emptyByChannel{channel_i} = catChannelEmpty;
  end
  spikesByCategory{cat_i} = cat_spikesByChannel;
  psthEmptyByCategory{cat_i} = cat_emptyByChannel;
end
disp('size spikes by category line 205');
disp(size(spikesByCategory));

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

%   now, read the LFP data, align by trial, and sort by image and category  %

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
lfpByImage = cell(length(picFiles),1);
for image_i = 1:length(picFiles)  
  onsets = int32(onsetsByImage{image_i});
  lfpArray = zeros(1,length(ns5channels),length(onsets),samPerMS*(psthPre+psthImDur+psthPost)+1); %(1,channel,trial,sample)
  for trial_i = 1:length(onsets)
    lfpArray(1,:,trial_i,:) = ns5data(:,samPerMS*(onsets(trial_i)-psthPre):samPerMS*(onsets(trial_i)+psthImDur+psthPost));
    if DCSUB_SAM(1)
      for ch_i = 1:size(lfpArray,2)
        lfpArray(1,ch_i,trial_i,:) = lfpArray(1,ch_i,trial_i,:) - mean(lfpArray(1,ch_i,trial_i,1:DCSUB_SAM(1)),4); %set the first 20ms mean equal across trials and stimuli
      end
    end
    if DCSUB_SAM(2)
      for ch_i = 1:size(lfpArray,2)
        endVal = mean(lfpArray(1,ch_i,trial_i,end-DCSUB_SAM(2):end),4);
        lfpArray(1,ch_i,trial_i,:) = squeeze(lfpArray(1,ch_i,trial_i,:)) - linspace(0,endVal,size(lfpArray,4))'; %subtract the linear term to set final values to zero
      end
    end
  end
  lfpByImage{image_i} = lfpArray;
end


% now, by category
lfpByCategory = cell(length(categoryList),1);
for cat_i = 1:length(categoryList)
  onsets = int32(onsetsByCategory{cat_i});
  lfpArray = zeros(1,length(ns5channels),length(onsets),samPerMS*(psthPre+psthImDur+psthPost)+1);
  for trial_i = 1:length(onsets)
    lfpArray(1,:,trial_i,:) = ns5data(:,samPerMS*(onsets(trial_i)-psthPre):samPerMS*(onsets(trial_i)+psthImDur+psthPost));
    if DCSUB_SAM(1)
      for ch_i = 1:size(lfpArray,2)
        lfpArray(1,ch_i,trial_i,:) = lfpArray(1,ch_i,trial_i,:)-mean(lfpArray(1,ch_i,trial_i,1:20),4);
      end
    end
    if DCSUB_SAM(2)
      for ch_i = 1:size(lfpArray,2)
        endVal = mean(lfpArray(1,ch_i,trial_i,end-DCSUB_SAM(2):end),4);
        lfpArray(1,ch_i,trial_i,:) = squeeze(lfpArray(1,ch_i,trial_i,:)) - linspace(0,endVal,size(lfpArray,4))'; %subtract the linear term to set final values to zero
      end
    end
  end
  lfpByCategory{cat_i} = lfpArray;
end
for cat_i = 1:length(categoryList)
  disp(categoryList{cat_i});
  disp(size(lfpByCategory{cat_i}));
end

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

%   now, build spike and lfp psth matrices, by image, category, and orientation  %

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

% get image filenames without path, extension, or underscores
picFilesPretty = picFiles;
for i = 1:length(picFiles)
  tmp = regexp(picFiles{i}, '\','split');
  tmp2 = regexp(tmp{end},'.bm','split');
  picFilesPretty{i} = regexprep(tmp2{1},'_','');
end

%spike psth color plot
if makeImPSTH
  for channel_i = 1:length(spikesByImage{1})
    for unit_i = 1:length(spikesByImage{1}{channel_i}) %TODO: BUG! if no spikes on first stimulus, can skip unit entirely (12/12/16: still true?)
      quickImagePSTH = zeros(length(picFiles),psthPre+psthImDur+psthPost);
      for image_i = 1:length(picFiles)
        if ~psthEmptyByImage{image_i}{channel_i}{unit_i}            
          quickImagePSTH(image_i,:) = 1000*psth(spikesByImage{image_i}{channel_i}{unit_i},smoothingWidth,'n',[1 psthPre+psthImDur+psthPost],0,1:psthPre+psthImDur+psthPost);
        end
      end
      xrange= [psthPre psthImDur+psthPost];
      yrange= [1 length(picFiles)];
      figure('Name',strcat(channel_names{channel_i}, ', unit #',num2str(unit_i)),'NumberTitle','off');
      title(strcat(channel_names{channel_i}, ', unit #',num2str(unit_i))); 
      caxis([min(min(quickImagePSTH)),max(max(quickImagePSTH))]);
      imagesc(xrange, yrange, quickImagePSTH); h = colorbar; ylabel(h,'Hz','FontSize',14);
      colormap(mystyle);
      ylimits= ylim();
      yRange = ylimits(2) - ylimits(1);  %todo: remove egregious use of this variable twice, but with different caps
      hold on
      draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
      draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
      set(gca,'YTick',linspace(ylimits(1)+yRange/(2*length(picFiles)),ylimits(2)-yRange/(2*length(picFiles)),length(picFiles)),'YTicklabel',picFilesPretty,...
        'box','off','TickDir','out','FontSize',14,'TickLength',[.012 .012],'LineWidth',1.3);
      xlabel('Time from stimulus onset (ms)', 'FontSize',14);
      hold off
      if save_fig     
        savefig(strcat(figures_dir, sprintf('PSTH_%s_Unit%d.fig',channel_names{channel_i},unit_i)));
        export_fig([figures_dir sprintf('PSTH_%s_Unit%d.png',channel_names{channel_i},unit_i)],'-m1.2','-transparent','-opengl');
      end
    end
  end
end
if makeCatPSTH
  for channel_i = 1:length(spikesByImage{1})
    for unit_i = 1:length(spikesByImage{1}{channel_i})  %TODO: BUG! if no spikes on first category, can skip unit entirely (12/12/16: still a thing?)
      quickCatPSTH = zeros(length(categoryList),psthPre+psthImDur+psthPost);
      for cat_i = 1:length(categoryList)
        if ~psthEmptyByCategory{cat_i}{channel_i}{unit_i}         
          quickCatPSTH(cat_i,:) = 1000*psth(spikesByCategory{cat_i}{channel_i}{unit_i},smoothingWidth,'n',[1 psthPre+psthImDur+psthPost],0,1:psthPre+psthImDur+psthPost);
        end
      end
      xrange= [psthPre psthImDur+psthPost];
      yrange= [1 length(categoryList)];
      figure('Name',strcat(channel_names{channel_i}, ', unit #',num2str(unit_i)),'NumberTitle','off');
      title(strcat(channel_names{channel_i}, ', unit #',num2str(unit_i))); 
      if verbose
        disp([min(min(quickCatPSTH)),max(max(quickCatPSTH))]);
      end
      caxis([min(min(quickCatPSTH)),max(max(quickCatPSTH))]);
      imagesc(xrange, yrange, quickCatPSTH); h = colorbar; ylabel(h,'Hz', 'FontSize',14);
      colormap(mystyle);
      ylimits= ylim();
      yRange = ylimits(2) - ylimits(1);  %todo: remove egregious use of this variable twice, but with different caps
      hold on
      draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
      draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
      set(gca,'YTick',linspace(ylimits(1)+yRange/(2*length(categoryList)),ylimits(2)-yRange/(2*length(categoryList)),length(categoryList)),'YTicklabel',categoryList,...
        'box','off','TickDir','out','FontSize',14,'TickLength',[.012 .012],'LineWidth',1.3);
      xlabel('Time from stimulus onset (ms)', 'FontSize',14);
      hold off
      if save_fig 
        savefig(strcat(figures_dir, sprintf('catPSTH_%s_Unit%d.fig',channel_names{channel_i},unit_i)));
        export_fig([figures_dir sprintf('catPSTH_%s_Unit%d.png',channel_names{channel_i},unit_i)],'-m1.2','-transparent','-opengl');    
      end
    end
  end
end

% firing rate RF color plot
rfGrid  = unique(taskData.pictureJumps,'rows');
gridX = unique(rfGrid(:,1));
gridsize = 2*mean(gridX(2:end,1)-gridX(1:end-1,1));
% these are the x and y values at which to interpolate RF values
xi = linspace(min(rfGrid(:,1)),max(rfGrid(:,1)),200);
yi = linspace(min(rfGrid(:,2)),max(rfGrid(:,2)),200);
if taskData.RFmap  %todo: put in check for analysis select conditional
  for channel_i = 1:length(spikesByImage{1})
    for unit_i = 1:length(spikesByImage{1}{channel_i}) %TODO: BUG! if no spikes on first stimulus, can skip unit entirely (12/12/16: still true?)
      meanRF = zeros(length(rfGrid),1);
      for image_i = 1:length(picFiles)
        imageRF = zeros(length(rfGrid),1);
        spikeLatencyRF = zeros(length(rfGrid),1);
        %also calc evoked peak time? power-weighted latency? peak of mean evoked correlogram?
        evokedPowerRF = zeros(length(rfGrid),1);
        for grid_i = 1:length(rfGrid)
          gridPointTrials = ismember(jumpsByImage{image_i},rfGrid(grid_i,:),'rows');
          if sum(gridPointTrials) == 0 
            imageRF(grid_i) = 0;
            disp('no data for this image and location');
            continue;
          end
          trialSpikes = spikesByImage{image_i}{channel_i}{unit_i}(gridPointTrials);
          totalSpikes = 0;
          spikeLatency = 0;
          evokedPotential = zeros(size(lfpByImage{image_i}(1,channel_i,1,samPerMS*(frCalcOn+psthPre):samPerMS*(frCalcOff+psthPre)))); 
          gridPointTrialInds = 1:length(gridPointTrials);
          gridPointTrialInds = gridPointTrialInds(gridPointTrials);
          for trial_i = 1:length(trialSpikes)
            totalSpikes = totalSpikes + sum(trialSpikes(trial_i).times > frCalcOn & trialSpikes(trial_i).times < frCalcOff);
            % bad latency measure; try weighting spike times by 1/ISI
            spikeLatency = spikeLatency + mean(trialSpikes(trial_i).times(trialSpikes(trial_i).times > frCalcOn & trialSpikes(trial_i).times < frCalcOff));
            evokedPotential = evokedPotential + lfpByImage{image_i}(1,channel_i,gridPointTrialInds(trial_i),samPerMS*(frCalcOn+psthPre):samPerMS*(frCalcOff+psthPre));
          end
          imageRF(grid_i) = (1000/(frCalcOff-frCalcOn))*totalSpikes/length(trialSpikes);
          spikeLatencyRF(grid_i) = 1/sum(gridPointTrials)*spikeLatency;
          evokedPotential = 1/sum(gridPointTrials)*(evokedPotential - mean(evokedPotential));
          evokedPowerRF(grid_i) = sum(abs(evokedPotential)); % should be squared rather than abs?
        end
        meanRF = meanRF + imageRF;
        display_map(rfGrid(:,1),rfGrid(:,2),imageRF,xi,yi,2.2857*gridsize,0,save_fig,sprintf('Channel %d, Unit %d, %s RF',channel_i,unit_i,picFilesPretty{image_i}),...
          [figures_dir sprintf('RF_%s_Unit%d_%s_Run%s.png',channel_names{channel_i},unit_i,picFilesPretty{image_i},runNum)]);
        if calcLatencyRF
          display_map(rfGrid(:,1),rfGrid(:,2),spikeLatencyRF,xi,yi,2.2857*gridsize,0,save_fig,sprintf('Channel %d, Unit %d, %s Latency RF',channel_i,unit_i,picFilesPretty{image_i}),...
            [figures_dir sprintf('LatencyRF_%s_Unit%d_%s_Run%s.png',channel_names{channel_i},unit_i,picFilesPretty{image_i},runNum)]);
        end
        if calcEvokedPowerRF
          display_map(rfGrid(:,1),rfGrid(:,2),evokedPowerRF,xi,yi,2.2857*gridsize,0,save_fig,sprintf('Channel %d, Unit %d, %s Evoked Power RF',channel_i,unit_i,picFilesPretty{image_i}),...
            [figures_dir sprintf('EvokedPowerRF_%s_Unit%d_%s_Run%s.png',channel_names{channel_i},unit_i,picFilesPretty{image_i},runNum)]);
        end
        % todo: add background subtracted version
        % todo: add subplots version
        % todo: coherency RFs (cc, cpt, ptpt)
        % todo: granger RF
        % todo: bandpassed power RFs
        % todo: stimulus decodability RF
      end
    end
  end
  display_map(rfGrid(:,1),rfGrid(:,2),meanRF,xi,yi,2.2857*gridsize,0,save_fig,sprintf('Channel %d, Unit %d, Mean RF',channel_i,unit_i),...
    [figures_dir sprintf('MeanRF_%s_Unit%d_Run%s.png',channel_names{channel_i},unit_i,runNum)]);
end  


%%%% make evoked potential plots by category

faceCatNum = find(strcmp(categoryList,'face'));
nonfaceCatNum = find(strcmp(categoryList,'nonface'));
categorySlimInds = zeros(length(categoryListSlim),1);
for i = 1:length(categoryListSlim)
  categorySlimInds(i) = find(strcmp(categoryList,categoryListSlim{i}));
end
colors = ['b','c','y','g','m','r','k'];
chColors = ['b','g','m'];
evokedFace = zeros(length(ns5channels),size(lfpArray,4));
evokedNonface = zeros(length(ns5channels),size(lfpArray,4));
%%%

for channel_i = 1:length(ns5channels)
  % TW=3 with T=.2, then W = 15 Hz (5 tapers)
  % TW=1.5 with T=.1, then W = 15 Hz (2 tapers)
  % TW = 1.5 with T=.2, then W = 7.5 Hz (2 tapers)
  chr_params.tapers = [1 1]; %[3 5] is chronux default; 
  chr_params.pad = 1;
  chr_params.fs = 1;
  chr_params.trialave = 1;
  chr_params.err = [1 .05];
  chr_params.fpass = [0 .1];
  movingWin = [50 5];
  specgramRowAve = 0;
  
  % create evoked-psth lineplot, face vs. non, one pane
  if channel_i == 1
    psthEvokedFig = figure();
    subplot(2,1,1);
    title('Face');
    xlabel('time (ms)');
    yyaxis right
    ylabel('power (normalized)'); % todo: units
    yyaxis left
    ylabel('firing rate (Hz)');
    hold on
    subplot(2,1,2);
    title('Non-face');
    xlabel('time (ms)');
    yyaxis right
    ylabel('power (normalized)'); % todo: units
    yyaxis left
    ylabel('firing rate (Hz)');
    hold on
    handlesForLegend1 = [];
    handlesForLegend2 = [];
    forLegend = {};
  end
  [faceV,faceT,faceE] = evoked(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',1000,[.001,.6],0,'n');
  [nfaceV,nfaceT,nfaceE] = evoked(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',1000,[.001,.6],0,'n');
% todo: fix length mismatch; evokedFace is 601, faceV is 598
%   disp(size(evokedFace(channel_i,:)));
%   disp(size(faceV));
%   evokedFace(channel_i,:) = faceV;
%   evokedNonface(channel_i,:) = nfaceV;
  % face vs. nonface evoked potential plot, one channel
  figure();
  plot(1000*faceT,faceV/max(nfaceV),1000*nfaceT, nfaceV/max(nfaceV), 'linewidth',3);
  hold on
  plot([0, psthImDur],[1.05*min(min(faceV),min(nfaceV))/max(nfaceV), 1.05*min(min(faceV),min(nfaceV))/max(nfaceV)],'color','k','linewidth',3);
  legend('face','nonface');
  title(sprintf('%s evoked potential',channel_names{channel_i}));
  xlabel('time (ms)', 'FontSize',18);
  ylabel('lfp (uV)', 'FontSize',18);
  set(gca,'fontsize',18);
  if save_fig    
    savefig(strcat(figures_dir, sprintf('Evoked_faceVnon_%s_Run%s.fig',channel_names{channel_i},runNum)));
    export_fig([figures_dir sprintf('Evoked_faceVnon_%s_Run%s.png',channel_names{channel_i},runNum)],'-m1.2','-transparent','-opengl');     % todo: add save to all these figs      
  end
  
  % contribute to evoked-psth lineplot, face vs. non, one pane
  figure(psthEvokedFig);
  subplot(2,1,1);
  yyaxis right
  lhLFP = plot(1000*faceT,faceV/max(faceV),'color',chColors(channel_i), 'linestyle','-','linewidth',2);  %todo: add errorbars?
  yyaxis left
  facePSTH = 1000*psth(spikesByCategory{faceCatNum}{channel_i}{end},smoothingWidth,'n',[1 psthPre+psthImDur+psthPost],0,1:psthPre+psthImDur+psthPost);
  lhMUA = plot(-1*psthPre:length(facePSTH)-psthPre-1,facePSTH,'color',chColors(channel_i),'linestyle','--','linewidth',2); 
  handlesForLegend1 = [handlesForLegend1,lhLFP,lhMUA];
  forLegend = [forLegend,strcat(channel_names{channel_i},' LFP'),strcat(channel_names{channel_i},' MUA')];
  if channel_i == length(ns5channels)
    legend(handlesForLegend1,forLegend);
    h = get(gca,'ylim');
    plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3,'linestyle','-');
  end
  subplot(2,1,2);
  yyaxis right
  lhLFP = plot(1000*nfaceT,nfaceV/max(nfaceV),'color',chColors(channel_i), 'linestyle', '-','linewidth',2);
  yyaxis left
  nfacePSTH = 1000*psth(spikesByCategory{nonfaceCatNum}{channel_i}{end},smoothingWidth,'n',[1 psthPre+psthImDur+psthPost],0,1:psthPre+psthImDur+psthPost);
  lhMUA = plot(-1*psthPre:length(nfacePSTH)-psthPre-1,nfacePSTH,'color',chColors(channel_i),'linestyle','--','linewidth',2);
  handlesForLegend2 = [handlesForLegend2, lhLFP, lhMUA];
  if channel_i == length(ns5channels)
    legend(handlesForLegend2,forLegend);
    h = get(gca,'ylim');
    plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3,'linestyle','-');
  end
  
  
  %lfp category psth
  figure();
  for i = 1:length(categoryListSlim)
    plot(evoked(squeeze(lfpByCategory{categorySlimInds(i)}(:,channel_i,:,:))',1000,[.001,.59],0,'n'), 'color',colors(i), 'linewidth',3);
    hold on
  end
  legend(categoryListSlim);
  h = get(gca,'ylim');
  plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
  hold off
  title(sprintf('%s evoked potentials',channel_names{channel_i}), 'FontSize',18);
  xlabel('time after stimulus (ms)', 'FontSize',18);
  ylabel('lfp (uV)', 'FontSize',18);
  set(gca,'fontsize',18);
  if save_fig    
    savefig(strcat(figures_dir,sprintf('Evoked_byCat_%s_Unit%d_Run%s.fig',channel_names{channel_i},unit_i,runNum))); 
    export_fig([figures_dir sprintf('Evoked_byCat_%s_Unit%d_Run%s.png',channel_names{channel_i},unit_i,runNum)],'-m1.2','-transparent','-opengl');     % todo: add save to all these figs      
  end
  
  % make lfp - psth subplot
  unit_i = length(spikesByImage{1}{channel_i}); % for now, only make this plot for MUA
  quickCatPSTH = zeros(length(categoryListSlim),psthPre+psthImDur+psthPost);
  for catsSlim_i = 1:length(categoryListSlim)
    cat_i = categorySlimInds(catsSlim_i);
    if ~psthEmptyByCategory{cat_i}{channel_i}{unit_i}        
      quickCatPSTH(catsSlim_i,:) = 1000*psth(spikesByCategory{cat_i}{channel_i}{unit_i},smoothingWidth,'n',[1 psthPre+psthImDur+psthPost],0,1:psthPre+psthImDur+psthPost);
    end
  end
  xrange= [psthPre psthImDur+psthPost];
  yrange= [1 length(categoryList)];
  figure();
  ah1 = subplot(2,1,1);
  caxis([min(min(quickCatPSTH)),max(max(quickCatPSTH))]);
  imagesc(xrange, yrange, quickCatPSTH); %h = colorbar; 
  h = colorbar();
  ylabel(h,'Hz', 'FontSize',14);
  colormap(mystyle);
  ylimits= ylim();
  yRange = ylimits(2) - ylimits(1);  %todo: remove egregious use of this variable twice, but with different caps
  hold on
  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
  set(gca,'YTick',linspace(ylimits(1)+yRange/(2*length(categoryListSlim)),ylimits(2)-yRange/(2*length(categoryListSlim)),length(categoryListSlim)),'YTicklabel',categoryListSlim,...
    'box','off','TickDir','out','FontSize',14,'TickLength',[.012 .012],'LineWidth',1.3);
  xlabel('Time from stimulus onset (ms)', 'FontSize',14);
  title(sprintf('%s evoked potentials and psth', channel_names{channel_i}), 'FontSize',18); 
  yyaxis right
  plot(-1*psthPre:length(facePSTH)-psthPre-1,sum(quickCatPSTH,1),'Color',[0.8,0.8,0.9],'LineWidth',4);
  ylabel('firing rate (Hz)');
  hold off
  ah2 = subplot(2,1,2);
  for i = 1:length(categoryListSlim)
    plot(evoked(squeeze(lfpByCategory{categorySlimInds(i)}(1,channel_i,:,:))',1000,[.001,.59],0,'n'), 'color',colors(i), 'linewidth',3);
    hold on
  end
  h = get(gca,'ylim');
  plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
  hold off
  legend(categoryListSlim,'Location','northeastoutside','FontSize',10);
  xlabel('time after stimulus (ms)', 'FontSize',14);
  ylabel('lfp (uV)', 'FontSize',14);
  set(gca,'fontsize',14);
  % align x axes with colorbar; see http://stackoverflow.com/questions/5259853/matlab-how-to-align-the-axes-of-subplots-when-one-of-them-contains-a-colorbar
  % note: drawnow line synchronizes the rendering thread with the main
  drawnow;
  pos1 = get(ah1,'Position');
  pos2 = get(ah2,'Position');
  pos2(1) = max(pos1(1),pos2(1)); % right limit
  pos1(1) = max(pos1(1),pos2(1));
  pos2(3) = min(pos1(3),pos2(3)); % left limit
  pos1(3) = min(pos1(3),pos2(3));
  set(ah2,'Position',pos2);
  set(ah1,'Position',pos1);
  if save_fig
    savefig(strcat(figures_dir,sprintf('PSTH_Evoked_%s_Run%s.fig',channel_names{channel_i},runNum)));
    export_fig([figures_dir sprintf('PSTH_Evoked_%s_Run%s.png',channel_names{channel_i},runNum)],'-m1.2','-transparent','-opengl');      
  end
  
  % single trial evoked lfps; probably better as subplot
  figure();
  subplot(2,1,1)
  hold on
  for i = 1:50
    plot(1000*faceT, squeeze(lfpByCategory{faceCatNum}(1,channel_i,i,1:length(faceT))));
  end;
  h = get(gca,'ylim');
  plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
  hold off
  title(sprintf('single trial face LFPs, %s', channel_names{channel_i}), 'FontSize',18);
  xlabel('time after stimulus (ms)', 'FontSize',18);
  ylabel('voltage (uV)', 'FontSize',18);
  set(gca,'fontsize',18);
  xlim([1000*min(faceT) 1000*max(faceT)]);
  
  subplot(2,1,2);
  hold on
  for i = 1:50
    plot(1000*nfaceT, squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,i,1:length(nfaceT))));
  end;
  h = get(gca,'ylim');
  plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
  hold off
  title(sprintf('single trial nonface LFPs, %s', channel_names{channel_i}), 'FontSize',18);
  xlabel('time after stimulus (ms)', 'FontSize',18);
  ylabel('voltage (uV)', 'FontSize',18);
  xlim([1000*min(nfaceT) 1000*max(nfaceT)]);
  if save_fig    
    savefig(strcat(figures_dir,sprintf('Evoked_singleTrials_%s_Run%s.fig',channel_names{channel_i},runNum)));
    export_fig([figures_dir sprintf('Evoked_singleTrials_%s_Run%s.png',channel_names{channel_i},runNum)],'-m1.2','-transparent','-opengl');     % todo: add save to all these figs      
  end
  % spectra
  [faceS,faceF, faceE] = mtspectrumc(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))', chr_params);
  [nfaceS,nfaceF, nfaceE] = mtspectrumc(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))', chr_params);
  figure();
  plot(1000*faceF,log(faceS),'linewidth',3,'marker','o');
  hold on
  plot(1000*nfaceF,log(nfaceS),'linewidth',3,'color','r','marker','o');
  hold off
  title(sprintf('%s evoked power spectrum',channel_names{channel_i}), 'FontSize',18);
  xlabel('frequency (Hz)', 'FontSize',18);
  ylabel('voltage, log(uV)', 'FontSize',18);
  legend('face','non');
  set(gca,'fontsize',18);
  if save_fig    
    savefig(strcat(figures_dir,sprintf('spectrum_faceVnon_%s_Run%s.fig',channel_names{channel_i},runNum)));
    export_fig([figures_dir sprintf('spectrum_faceVnon_%s_Run%s.png',channel_names{channel_i},runNum)],'-m1.2','-transparent','-opengl');     % todo: add save to all these figs      
  end
  figure();
  loglog(1000*faceF,faceS,'linewidth',3,'marker','o');
  hold on
  loglog(1000*nfaceF,nfaceS,'linewidth',3,'color','r','marker','o');
  hold off
  title(sprintf('%s evoked power spectrum, loglog',channel_names{channel_i}), 'FontSize',18);
  xlim([1000*min(faceF) 1000*max(faceF)]); %todo: fix error
  xlabel('frequency (Hz)', 'FontSize',18);
  ylabel('voltage (uV)', 'FontSize',18);
  legend('face','non');
  set(gca,'fontsize',18);
  if save_fig    
    savefig(strcat(figures_dir,sprintf('spectrum_log_faceVnon_%s_Run%s.fig',channel_names{channel_i},runNum)));
    export_fig([figures_dir sprintf('spectrum_log_faceVnon_%s_Run%s.png',channel_names{channel_i},runNum)],'-m1.2','-transparent','-opengl');     % todo: add save to all these figs      
  end
  
  
  %%%%%% TODO %%%%%%
  % todo: multichannel evoked lineplot, categories, subplots
  % todo: find local maxes and mins in evoked and psth
  % psth category lineplot (probably in subplot with evoked by category)

  % image-wise time-frequency plots, spikes and lfp
  if imageTF
    for image_i = 1:length(picFiles)
      [S,t,f,R]=mtspecgrampt(spikesByImage{image_i}{channel_i}{end},movingWin,chr_params); %last optional param not used: fscorr
      figure();
      imagesc(t,f,S); %TODO: fix to match cat version, which is correct!!
      set(gca,'Ydir','normal'); % also need 'Yscale,'log'??
      xlabel('Time'); % todo: units
      ylabel('Frequency (Hz)');
      c = colorbar();
      ylabel(c,'Power'); % todo: check units
      title(sprintf('%s MUA Time-Frequency, %s',channel_names{channel_i},picFilesPretty{image_i}));
      if save_fig
        savefig(strcat(figures_dir,sprintf('TF_MUA_%s_%s.fig',channel_names{channel_i},picFilesPretty{image_i})));
        export_fig([figures_dir sprintf('TF_MUA_%s_%s.png',channel_names{channel_i},picFilesPretty{image_i})],'-m1.2','-transparent','-opengl'); 
        close();
      end
      [S,t,f]=mtspecgramc(squeeze(lfpByImage{image_i}(1,channel_i,:,:))',movingWin,chr_params);
      figure();
      imagesc(t,f,S);
      set(gca,'Ydir','normal'); % also need 'Yscale,'log'??
      xlabel('Time'); % todo: units
      ylabel('Frequency (Hz)');
      c = colorbar();
      ylabel(c,'Power'); % todo: check units
      title(sprintf('%s LFP Time-Frequency, %s',channel_names{channel_i},picFilesPretty{image_i}));
      if save_fig
        savefig(strcat(figures_dir,sprintf('TF_LFP_%s_%s.fig',channel_names{channel_i},picFilesPretty{image_i})));
        export_fig([figures_dir sprintf('TF_LFP_%s_%s.png',channel_names{channel_i},picFilesPretty{image_i})],'-m1.2','-transparent','-opengl'); 
        close();
      end
    end
  end
  % category-wise time-frequency plots, spikes and lfp
  
  if catTF
    if channel_i == 1
      lfpTfFig = figure();
      muaTfFig = figure();
    end
    for cat_i = 1:length(categoryList)
      if cat_i ~= faceCatNum && cat_i ~= nonfaceCatNum
        continue
      end
      [S,t,f,R]=mtspecgrampt(spikesByCategory{cat_i}{channel_i}{end},movingWin,chr_params); %last optional param not used: fscorr
      totalPower = sum(S,2)';
      figure();
      if specgramRowAve
        for i = 1:size(S,2)
          S(:,i) = S(:,i)/mean(S(:,i)); 
        end
        imagesc(t,1000*f,S'); axis xy; c = colorbar();
        ylabel(c,'Row-Normalized Power');
      else
        S = 10*log10(S);
        imagesc(t,1000*f,S'); axis xy; c = colorbar();
        ylabel(c,'Power (dB)');
      end
      xlabel('Time (ms)'); 
      ylabel('Frequency (Hz)');
      hold on
      yyaxis right
      plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
      ylabel('Integrated Power');
      hold off
      title(sprintf('%s MUA Time-Frequency, %s',channel_names{channel_i},categoryList{cat_i}));
      if save_fig
        savefig(strcat(figures_dir,sprintf('TF_MUA_%s_%s.fig',channel_names{channel_i},categoryList{cat_i})));
        export_fig([figures_dir sprintf('TF_MUA_%s_%s.png',channel_names{channel_i},categoryList{cat_i})],'-m1.2','-transparent','-opengl'); 
        close();
      end
      % contribute to shared figure
      figure(muaTfFig);
      if cat_i == faceCatNum
        subplot(3,2,2*channel_i-1);
        forTitle = sprintf('Face MUA TF, %s',channel_names{channel_i});
      else
        subplot(3,2,2*channel_i);
        forTitle = sprintf('Nonface MUA TF, %s',channel_names{channel_i});
      end
      imagesc(t,1000*f,S'); axis xy; c = colorbar(); ylabel(c,'Power (dB)'); %todo: fix unit when row-normalizing
      xlabel('Time (ms)'); 
      ylabel('Frequency (Hz)');
      title(forTitle);
      
      [S,t,f]=mtspecgramc(squeeze(lfpByCategory{cat_i}(1,channel_i,:,:))',movingWin,chr_params);
      totalPower = sum(S,2)';
      
      figure();
      if specgramRowAve
        for i = 1:size(S,2)
          S(:,i) = S(:,i)/mean(S(:,i));
        end
        imagesc(t,1000*f,S'); axis xy; c = colorbar();
        ylabel(c,'Row-Normalized Power');
      else
        S = 10*log10(S);
        imagesc(t,1000*f,S'); axis xy; c = colorbar();
        ylabel(c,'Power (dB)');
      end
      xlabel('Time (ms)'); 
      ylabel('Frequency (Hz)');
      hold on
      yyaxis right
      plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
      ylabel('Integrated Power');
      hold off
      title(sprintf('%s LFP Time-Frequency, %s',channel_names{channel_i},categoryList{cat_i}));
      if save_fig
        savefig(strcat(figures_dir,sprintf('TF_LFP_%s_%s.fig',channel_names{channel_i},categoryList{cat_i})));
        export_fig([figures_dir sprintf('TF_LFP_%s_%s.png',channel_names{channel_i},categoryList{cat_i})],'-m1.2','-transparent','-opengl');
      end
      % contribute to shared figure
      figure(lfpTfFig);
      if cat_i == faceCatNum
        subplot(3,2,2*channel_i-1);
        forTitle = sprintf('Face LFP TF, %s',channel_names{channel_i});
      else
        subplot(3,2,2*channel_i);
        forTitle = sprintf('Nonface LFP TF, %s',channel_names{channel_i});
      end
      imagesc(t,1000*f,S'); axis xy; c = colorbar(); ylabel(c,'Power (dB)'); %todo: fix unit when row-normalizing
      xlabel('Time (ms)'); 
      ylabel('Frequency (Hz)');
      title(forTitle);
    end
    
    if save_fig && channel_i == length(ns5channels) && cat_i == max(faceCatNum, nonfaceCatNum)
      figure(muaTfFig);
      savefig(strcat(figures_dir,sprintf('TF_MUA.fig',channel_names{channel_i},categoryList{cat_i})));
      export_fig([figures_dir sprintf('TF_MUA.png',channel_names{channel_i},categoryList{cat_i})],'-m1.2','-transparent','-opengl');
      figure(lfpTfFig);
      savefig(lfpTfFig,strcat(figures_dir,sprintf('TF_LFP.fig',channel_names{channel_i},categoryList{cat_i})));
      export_fig([figures_dir sprintf('TF_LFP.png',channel_names{channel_i},categoryList{cat_i})],'-m1.2','-transparent','-opengl');
    end
  end
  
  drawnow;
  
  %%% across channels
  if crossTF
    for channel2_i = channel_i:length(ns5channels)
      % face vs. nonface spike field coherence, within and across channels
      % channel_i spike -- channel_i field
      useJacknife = 0;
      if useJacknife
        chr_params.err = [2 .05];
        % todo: eliminate 'faceMuaML' etc.
        [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
          spikesByCategory{faceCatNum}{channel_i}{end},chr_params);
        [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
          spikesByCategory{nonfaceCatNum}{channel_i}{end},chr_params);
        errs = vertcat(faceErrs,nfaceErrs);
      else
        [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
          spikesByCategory{faceCatNum}{channel_i}{end},chr_params);
        [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
          spikesByCategory{nonfaceCatNum}{channel_i}{end},chr_params);
        errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
      end
      figure();
      mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
      legend('face', 'nonface');
      xlabel('frequency (Hz)');
      ylabel('coherency');
      title(sprintf('%s spike - %s field coherence',channel_names{channel_i},channel_names{channel_i}),'FontSize',18);
      if save_fig
        export_fig([figures_dir sprintf('coh_MUA-LFP_%s_%s_FACEvsNON.png',channel_names{channel_i},channel_names{channel_i})],'-m1.2','-transparent','-opengl'); 
        close();
      end
      if channel2_i > channel_i 
        % channel2_i spike -- channel_i field
        useJacknife = 0;
        if useJacknife
          chr_params.err = [2 .05];
          [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
            spikesByCategory{faceCatNum}{channel2_i}{end},chr_params);
          [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
            spikesByCategory{nonfaceCatNum}{channel2_i}{end},chr_params);
          errs = vertcat(faceErrs,nfaceErrs);
        else
          [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
            spikesByCategory{faceCatNum}{channel2_i}{end},chr_params);
          [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
            spikesByCategory{nonfaceCatNum}{channel2_i}{end},chr_params);
          errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
        end
        figure();
        mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
        legend('face', 'nonface');
        xlabel('frequency (Hz)');
        ylabel('coherency');
        title(sprintf('%s spike - %s field coherence',channel_names{channel2_i},channel_names{channel_i}),'FontSize',18);
        if save_fig
          export_fig([figures_dir sprintf('coh_MUA-LFP_%s_%s_FACEvsNON.png',channel_names{channel2_i},channel_names{channel_i})],'-m1.2','-transparent','-opengl'); 
          close();
        end
        % channel_i spike -- channel2_i field
        useJacknife = 0;
        if useJacknife
          chr_params.err = [2 .05];
          [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
            spikesByCategory{faceCatNum}{channel_i}{end},chr_params);
          [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
            spikesByCategory{nonfaceCatNum}{channel_i}{end},chr_params);
          errs = vertcat(faceErrs,nfaceErrs);
        else
          [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
            spikesByCategory{faceCatNum}{channel_i}{end},chr_params);
          [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
            spikesByCategory{nonfaceCatNum}{channel_i}{end},chr_params);
          errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
        end
        figure();
        mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
        legend('face', 'nonface');
        xlabel('frequency (Hz)');
        ylabel('coherency');
        title(sprintf('%s spike - %s field coherence',channel_names{channel_i},channel_names{channel2_i}),'FontSize',18);
        if save_fig
          export_fig([figures_dir sprintf('coh_MUA-LFP_%s_%s_FACEvsNON.png',channel_names{channel_i},channel_names{channel2_i})],'-m1.2','-transparent','-opengl'); 
          close();
        end
        

        % field-field
        [Cface,phiface,S12,S1,S2,fface,confCface,phistdface]=coherencyc(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
          squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',chr_params);
        [Cnface,phinface,S12,S1,S2,fnface,confCnface,phistdnface]=coherencyc(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
          squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',chr_params);
        errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
        figure();
        mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
        legend('face', 'nonface');
        xlabel('frequency (Hz)');
        ylabel('coherency');
        title(sprintf('%s field - %s field coherence',channel_names{channel_i},channel_names{channel2_i}),'FontSize',18);
        if save_fig
          export_fig([figures_dir sprintf('coh_LFP-LFP_%s_%s_FACEvsNON.png',channel_names{channel_i},channel_names{channel2_i})],'-m1.2','-transparent','-opengl'); 
          close();
        end
        % field-field time-frequency coherency, face
        [Cface,phiface,S12,S1,S2,tface,fface,confCface,phistdface]=cohgramc(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
          squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',movingWin,chr_params);
        figure();
        imagesc(t,f,S'); axis xy
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'Coherency');
        title(sprintf('%s field - %s field coherence',channel_names{channel_i},channel_names{channel2_i}),'FontSize',18);
        if save_fig
          export_fig([figures_dir sprintf('coh_LFP-LFP_%s_%s_FACE.png',channel_names{channel_i},channel_names{channel2_i})],'-m1.2','-transparent','-opengl'); 
          close();
        end
        % lfp-lfp time-frequency coherency, nonface
        [Cnface,phinface,S12,S1,S2,tnface,fnface,confCnface,phistdnface]=cohgramc(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
          squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',movingWin,chr_params);
        figure();
        imagesc(t,f,S'); axis xy
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'Coherency'); 
        title(sprintf('%s field - %s field coherence',channel_names{channel_i},channel_names{channel2_i}),'FontSize',18);
        if save_fig
          export_fig([figures_dir sprintf('coh_LFP-LFP_%s_%s_NONFACE.png',channel_names{channel_i},channel_names{channel2_i})],'-m1.2','-transparent','-opengl'); 
          close();
        end
        
%         % spike-spike coherency, face vs. non
%         [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypt(spikesByCategory{faceCatNum}{channel_i}{end},...
%           spikesByCategory{faceCatNum}{channel2_i}{end},chr_params,0,1);
%         [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(spikesByCategory{faceCatNum}{channel_i}{end},...
%           spikesByCategory{nonfaceCatNum}{channel2_i}{end},chr_params,0,1);
%         errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
%         figure();
%         mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
%         legend('face', 'nonface');
%         xlabel('frequency (Hz)');
%         ylabel('coherency');
%         title(sprintf('%s Spike - %s Spike coherence',channel_names{channel_i},channel_names{channel2_i}),'FontSize',18);
%         if save_fig
%           export_fig([figures_dir sprintf('coh_MUA-MUA_%s_%s_FACEvsNON.png',channel_names{channel_i},channel_names{channel2_i})],'-m1.2','-transparent','-opengl'); 
%           close();
%         end
%         % spike-spike time-frequency coherency, face vs. non
        
        % todo: granger (full, passbands)
      end
    end
  end 
end
end

