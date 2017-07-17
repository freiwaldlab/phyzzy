function [ spikesByImage, spikesByCategory, lfpByImage, lfpByCategory, categoryList, picFiles ] = processRun( varargin )
%processRun is a top level function to analyze units, MUA, lfps, coupling, and RFs. 
%   - handles single-channel and multi-channel sessions
%   - relies on raw visiko (.log) and blackrock (.ns5) files
%   - aligns visiko events to blackrock clock (preprocessLogFile.m)
%   - excludes trials with broken fixation, fix spot flash, or (optionally)
%     juice delivery (via exludeStimuli)
%   - decimates and (optionally) filters raw lfp data (default, to 1 kHz)
%   Inputs:
%   - varargin can have the following forms:
%       - empty (default assignments: buildAnalysisParamFile, runAnalyses
%       - 'paramBuilder', paramBuilderFilename
%       - 'analyzer', analyzerFilename
%       - 'paramBuilder', paramBuilderFilename, 'analyzer', analyzerFilename
%   Notes:
%   Depends:
%   - contents of 'dependencies' folder (details coming)
%   - R2016a (or later) if joint psth-evoked potential plots desired
%   - Signal Processing Toolbox (for dpss taper calculation, LFP filters)

addpath(genpath('dependencies/genpath_exclude'));
addpath(genpath_exclude('dependencies',{'*mvgc_v1.0'})); %note: use this to exclude libraries that overwrite matlab builtin namespaces, until they're needed
% load analysis parameters
if nargin == 0
  analysisParamFilename = buildAnalysisParamFile();
else
  if strcmp(varargin{1},'preprocessed')
    load(varargin{2}); 
  else
    if strcmp(varargin{1},'paramBuilder')
      analysisParamFilename = feval(varargin{2});
    else
      analysisParamFilename = varargin{1};
    end
  end 
end
if nargin == 0 || ~strcmp(varargin{1},'preprocessed')
  load(analysisParamFilename);
  % extract parameters needed in this function from structures
  channelNames = ephysParams.channelNames;
  spikeChannels = ephysParams.spikeChannels;
  lfpChannels = ephysParams.lfpChannels;
  psthPre = psthParams.psthPre;
  psthImDur = psthParams.psthImDur;
  psthPost = psthParams.psthPost;
  smoothingWidth = psthParams.smoothingWidth;
  preAlign = spikeAlignParams.preAlign;
  postAlign = spikeAlignParams.postAlign;
  spikeAlignParams.preAlign;
  movingWin = tfParams.movingWin;
  specgramRowAve = tfParams.specgramRowAve;
  samPerMS = ephysParams.samPerMS;
  % set verbose level
  if strcmp(verbosity,'VERBOSE')
    Output.level(Output.DISP_VERBOSE);
  end
  if strcmp(verbosity,'DEBUG')
    Output.level(Output.DISP_DEBUG);
  end

  %%%%%%%%%%%%%%%%%
  analogInData = preprocessAnalogIn(analogInFilename, analogInParams); %todo: implement return variables
  [spikesByChannel, taskTriggers] = preprocessSpikes(spikeFilename, ephysParams);
  lfpData = preprocessLFP(lfpFilename, ephysParams);
  [ taskData, stimTiming ] = preprocessLogFile(taskFilename, taskTriggers, stimSyncParams); % load visual stimulus data and transform its timestamps to ephys clock reference
  if stimTiming.shortest == stimTiming.longest
    psthParams.psthImDur = stimTiming.shortest;
    psthImDur = psthParams.psthImDur;
    spikeAlignParams.postAlign = psthImDur+psthPost+3*smoothingWidth;
    lfpAlignParams.msPostAlign = psthImDur+psthPost+tfParams.movingWin(1)/2;
    if length(lfpAlignParams.DCSUB_SAM) > 1 && any(lfpAlignParams.DCSUB_SAM(2,:))  % todo: handle case where lfp analysis window doesn't cover full ISI
      lfpAlignParams.DCSUB_SAM(2,:) = (psthImDur+stimTiming.ISI)+lfpAlignParams.DCSUB_SAM(2,:);
    end
    preAlign = spikeAlignParams.preAlign;
    postAlign = spikeAlignParams.postAlign;
  else
    assert(psthImDur < stimTiming.longest, 'psthImDur is longer than longest trial; nothing to analyze');
    assert(psthImDur > 0, 'psthImDur = 0; nothing to analyze. Likely cause: unexpected variable stimulus length');
    fprintf('Variable stimulus length run. Excluding trials shorter than %d\n',psthImDur);
  end
  if frCalcOff < frCalcOn
    frCalcOff = psthImDur+frCalcOn;
  end

  taskDataAll = taskData;
  % exclude stimuli for fixation out, flash on, frame dropped, (accel high, juice on)
  % params are ( taskData, fixPre, fixPost, flashPre, flashPost, varargin )
  taskData = excludeStimuli( taskData, excludeStimParams);

  %sort trials by image and image category
  tmp = load(picParamsFilename); %loads variables paramArray, categoryLabels,pictureLabels
  picCategories = tmp.paramArray;
  categoryList = tmp.categoryLabels;
  pictureLabels = tmp.pictureLabels;
  picFiles = {};
  for pic_i = 1:length(picCategories)
    picFiles = vertcat(picFiles,picCategories{pic_i}{1}); 
  end
  onsetsByImage = cell(length(picFiles),1);
  offsetsByImage = cell(length(picFiles),1);
  picsNotPresented = zeros(length(picFiles),1);
  for i = 1:length(picFiles)
    onsetsByImage{i} = taskData.pictureStartTimes(strcmp(taskData.pictureFilenames,picFiles{i}));
    offsetsByImage{i} = taskData.pictureEndTimes(strcmp(taskData.pictureFilenames,picFiles{i}));
    picsNotPresented(i) = isempty(onsetsByImage{i});
  end
  % todo: add similar defense for categories
  % todo: need defense against image with onset but no offset? 
  % todo: add similar defense for rf map locations here?
  if ~taskData.RFmap
    disp('No presentations of the following images survived exclusion:');
    disp(picFiles(picsNotPresented == 1));
  end
  onsetsByImage = onsetsByImage(picsNotPresented == 0);
  offsetsByImage = offsetsByImage(picsNotPresented == 0);
  picFiles = picFiles(picsNotPresented == 0);
  pictureLabels = pictureLabels(picsNotPresented == 0);
  picCategories = picCategories(picsNotPresented == 0);

  onsetsByCategory = cell(length(categoryList));
  offsetsByCategory = cell(length(categoryList));
  for cat_i = 1:length(categoryList)
    catOnsets = [];
    catOffsets = [];
    for image_i = 1:length(picFiles)
      if any(strcmp(picCategories{image_i},categoryList{cat_i}))
        catOnsets = vertcat(catOnsets,onsetsByImage{image_i});
        catOffsets = vertcat(catOffsets, offsetsByImage{image_i});
      end
    end
    onsetsByCategory{cat_i} = catOnsets;
    offsetsByCategory{cat_i} = catOffsets;
    Output.DEBUG('numel cat onsets');
    Output.DEBUG(numel(catOnsets));
  end

  if taskData.RFmap
    jumpsByImage = cell(length(picFiles),1);
    for i = 1:length(picFiles)
      jumpsByImage{i} = taskData.pictureJumps(strcmp(taskData.pictureFilenames,picFiles{i}),:);
    end
  else
    jumpsByImage = [];
  end
  %align spikes by trial, and sort by image and category
  spikeAlignParams.refOffset = 0;
  [spikesByImage, psthEmptyByImage] = alignSpikes( spikesByChannel, onsetsByImage, spikeChannels, spikeAlignParams );
  [spikesByCategory, psthEmptyByCategory] = alignSpikes( spikesByChannel, onsetsByCategory, spikeChannels, spikeAlignParams );

  % align spikes again, but this time reference time to lfp sample number (required for chronux TF, even spike-spike)
  spikeAlignParamsTF.preAlign = lfpAlignParams.msPreAlign;
  spikeAlignParamsTF.postAlign = lfpAlignParams.msPostAlign;
  spikeAlignParamsTF.refOffset = -lfpAlignParams.msPreAlign;
  spikesByImageForTF = alignSpikes( spikesByChannel, onsetsByImage, spikeChannels, spikeAlignParamsTF );
  spikesByCategoryForTF = alignSpikes( spikesByChannel, onsetsByCategory, spikeChannels, spikeAlignParamsTF );
  %  align LFP data  by trial, sort by image and category, and possibly remove DC and linear components
  lfpByImage = alignLFP(lfpData, onsetsByImage, lfpChannels, lfpAlignParams);
  lfpByCategory = alignLFP(lfpData, onsetsByCategory, lfpChannels, lfpAlignParams);

  for cat_i = 1:length(categoryList)
    Output.VERBOSE(categoryList{cat_i});
    Output.VERBOSE(size(lfpByCategory{cat_i}));
  end

  channelUnitNames = cell(length(spikeChannels),1);
  for channel_i = 1:length(spikesByImage{1})
    channelUnitNames{channel_i} = cell(length(spikesByImage{1}{channel_i}),1);
    channelUnitNames{channel_i}{1} = 'Unsorted';
    channelUnitNames{channel_i}{end} = 'MUA';
    for unit_i = 2:length(spikesByImage{1}{channel_i})-1
      channelUnitNames{channel_i}{unit_i} = sprintf('Unit %d',unit_i-1); 
    end
  end
  if savePreprocessed
    save(preprocessedDataFilename,'analysisParamFilename', 'spikesByChannel', 'lfpData', 'analogInData', 'taskData', 'taskDataAll', 'psthImDur', 'preAlign', 'postAlign',...
      'categoryList', 'pictureLabels', 'jumpsByImage', 'spikesByImage', 'psthEmptyByImage', 'spikesByCategory', 'psthEmptyByCategory',...
      'spikesByImageForTF', 'spikesByCategoryForTF', 'lfpByImage', 'lfpByCategory', 'channelUnitNames', 'stimTiming', 'picCategories', 'onsetsByImage', 'onsetsByCategory')
  end
end
if nargin == 0 || (nargin == 2 && strcmp(varargin{1},'paramBuilder')) || (nargin == 2 && strcmp(varargin{1},'preprocessed'))
  runAnalyses( analysisParamFilename, spikesByChannel, lfpData, analogInData, taskData, taskDataAll, psthImDur, preAlign, postAlign,...
    categoryList, pictureLabels, jumpsByImage, spikesByImage, psthEmptyByImage, spikesByCategory, psthEmptyByCategory,...
    spikesByImageForTF, spikesByCategoryForTF, lfpByImage, lfpByCategory, channelUnitNames, stimTiming, picCategories, onsetsByImage, onsetsByCategory);
else
  feval(varargin{end},analysisParamFilename, spikesByChannel, lfpData, analogInData, taskData, taskDataAll, psthImDur, preAlign, postAlign,...
    categoryList, pictureLabels, jumpsByImage, spikesByImage, psthEmptyByImage, spikesByCategory, psthEmptyByCategory,...
    spikesByImageForTF, spikesByCategoryForTF, lfpByImage, lfpByCategory, channelUnitNames, stimTiming, picCategories, onsetsByImage, onsetsByCategory);
end
end

