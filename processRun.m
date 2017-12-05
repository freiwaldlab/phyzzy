function [ spikesByImage, spikesByCategory, lfpByImage, lfpByCategory, categoryList, picFiles ] = processRun( varargin )
%processRun loads, preprocesses, and aligns task events, lfps, muas, units, and analog signals (eg eye, accelerometer,photodiodes)
%   - handles single-channel and multi-channel sessions
%   - relies only on raw task log and physiology fiels: currently visiko (.log) and blackrock (.ns5,.ns2)
%   - aligns task events to ephys system clock (preprocessLogFile.m)
%   - excludes trials with broken fixation, fix spot flash, and (optionally)
%     juice delivery and high head acceleration(via exludeStimuli)
%   - decimates and (optionally) filters raw lfp data (default, to 1 kHz)
%   - calibrates eye signals (via preprocessEyeSignals)
%   - calibrates accelerometer signals (via processAccelSignals)
%   - calculates spike isolation measures (via preprocessSpikes)
%   - shows task event summary (via excludeTrials)
%   - optionally calls an analysis routine; default is runAnalyses.m
%   Inputs:
%   - varargin can have the following forms:
%       - empty (default assignments: buildAnalysisParamFile, runAnalyses
%       - 'paramBuilder', paramBuilderFilename
%       - 'analyzer', analyzerFilename
%       - 'paramBuilder', paramBuilderFilename, 'analyzer', analyzerFilename
%   Notes:=
%   Depends:
%   - contents of 'dependencies' folder (details coming)
%   - R2016a (or later) if joint psth-evoked potential plots desired
%   - Signal Processing Toolbox (for dpss taper calculation, LFP filters)

addpath(genpath('dependencies/genpath_exclude'));
addpath(genpath_exclude('dependencies',{'*mvgc_v1.0'})); %note: use this to exclude libraries that overwrite matlab builtin namespaces, until they're needed
% load analysis parameters
if nargin == 0
  analysisParamFilename = buildAnalysisParamFileSofiFamFacesVoices();
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
  checkAnalysisParamFile(analysisParamFilename);
  load(analysisParamFilename);
  % extract parameters needed in this function from structures
  channelNames = ephysParams.channelNames;
  spikeChannels = ephysParams.spikeChannels;
  lfpChannels = ephysParams.lfpChannels;
  analogInChannels = analogInParams.analogInChannels;
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
  analogInData = preprocessAnalogIn(analogInFilename, analogInParams); 
  [spikesByChannel, taskTriggers, channelUnitNames] = preprocessSpikes(spikeFilename, ephysParams);
  %lfpData = preprocessLFP(lfpFilename, ephysParams);
  [ lfpData, diodeTrace, audioTrace ] = preprocessNS5data( lfpFilename, ephysParams, photodiodeParams, audioParams );
  diodeTriggers = preprocessPhotodiodeStrobe(photodiodeFilename, photodiodeParams);
  
  if ~stimSyncParams.useAudio
      [ taskData, stimTiming ] = preprocessLogFile(taskFilename, taskTriggers, diodeTriggers, stimSyncParams); % load visual stimulus data and transform its timestamps to ephys clock reference
  else
      [ taskData, stimTiming ] = preprocessLogFileAudio(taskFilename, audioTrace, stimSyncParams); % load visual stimulus data and transform its timestamps to ephys clock reference      
  end
  
  analogInData = preprocessEyeSignals(analogInData,taskData,eyeCalParams);
  analogInData = preprocessAccelSignals(analogInData, accelParams); 
  taskDataAll = taskData;
  %taskData = excludeTrials( taskData, excludeStimParams); %exclude trials for lost fixation etc. 
  
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

  %sort trials by image and image category
  tmp = load(stimParamsFilename); %loads variables paramArray, categoryLabels,pictureLabels
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
    onsetsByImage{i} = taskData.stimStartTimes(strcmp(taskData.stimFilenames,picFiles{i}));
    offsetsByImage{i} = taskData.stimEndTimes(strcmp(taskData.stimFilenames,picFiles{i}));
    picsNotPresented(i) = isempty(onsetsByImage{i});
  end

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
  catsNotPresented = zeros(length(categoryList),1);
  for cat_i = 1:length(categoryList)
    catOnsets = [];
    catOffsets = [];
    for image_i = 1:length(picFiles)
      if any(strcmp(picCategories{image_i},categoryList{cat_i}))
        catOnsets = horzcat(catOnsets,onsetsByImage{image_i});
        catOffsets = horzcat(catOffsets, offsetsByImage{image_i});
      end
    end
    onsetsByCategory{cat_i} = catOnsets;
    offsetsByCategory{cat_i} = catOffsets;
    catsNotPresented(cat_i) = isempty(onsetsByCategory{cat_i});
    Output.DEBUG('numel cat onsets');
    Output.DEBUG(numel(catOnsets));
  end
  
  onsetsByCategory = onsetsByCategory(catsNotPresented == 0);
  offsetsByCategory = offsetsByCategory(catsNotPresented == 0);
  categoryList = categoryList(catsNotPresented == 0);

  if taskData.RFmap
    jumpsByImage = cell(length(picFiles),1);
    for i = 1:length(picFiles)
      jumpsByImage{i} = taskData.stimJumps(strcmp(taskData.stimFilenames,picFiles{i}),:);
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
  analogInByImage = alignAnalogIn(analogInData, onsetsByImage, analogInChannels, lfpAlignParams);
  analogInByCategory = alignAnalogIn(analogInData, onsetsByCategory, analogInChannels, lfpAlignParams);
  
  for cat_i = 1:length(categoryList)
    Output.VERBOSE(categoryList{cat_i});
    Output.VERBOSE(size(lfpByCategory{cat_i}));
  end

  if savePreprocessed
    save(preprocessedDataFilename,'analysisParamFilename', 'spikesByChannel', 'lfpData', 'analogInData', 'taskData', 'taskDataAll', 'psthImDur', 'preAlign', 'postAlign',...
      'categoryList', 'pictureLabels', 'jumpsByImage', 'spikesByImage', 'psthEmptyByImage', 'spikesByCategory', 'psthEmptyByCategory',...
      'spikesByImageForTF', 'spikesByCategoryForTF', 'lfpByImage', 'lfpByCategory', 'analogInByImage','analogInByCategory','channelUnitNames', ...
      'stimTiming', 'picCategories', 'onsetsByImage', 'onsetsByCategory')
  end
end

runAnalysisInputs.analysisParamFilename = analysisParamFilename;
runAnalysisInputs.spikesByChannel = spikesByChannel; 
runAnalysisInputs.lfpData = lfpData; 
runAnalysisInputs.analogInData = analogInData; 
runAnalysisInputs.taskData = taskData;
runAnalysisInputs.taskDataAll = taskDataAll; 
runAnalysisInputs.psthImDur = psthImDur; 
runAnalysisInputs.preAlign = preAlign; 
runAnalysisInputs.postAlign = postAlign;
runAnalysisInputs.categoryList = categoryList; 
runAnalysisInputs.pictureLabels = pictureLabels; 
runAnalysisInputs.jumpsByImage = jumpsByImage; 
runAnalysisInputs.spikesByImage = spikesByImage; 
runAnalysisInputs.psthEmptyByImage = psthEmptyByImage;  
runAnalysisInputs.spikesByCategory = spikesByCategory; 
runAnalysisInputs.psthEmptyByCategory = psthEmptyByCategory; 
runAnalysisInputs.spikesByImageForTF = spikesByImageForTF;  
runAnalysisInputs.spikesByCategoryForTF = spikesByCategoryForTF;  
runAnalysisInputs.lfpByImage = lfpByImage;  
runAnalysisInputs.lfpByCategory = lfpByCategory;  
runAnalysisInputs.analogInByImage = analogInByImage; 
runAnalysisInputs.analogInByCategory = analogInByCategory;  
runAnalysisInputs.channelUnitNames = channelUnitNames;  
runAnalysisInputs.stimTiming = stimTiming;  
runAnalysisInputs.picCategories = picCategories;  
runAnalysisInputs.onsetsByImage = onsetsByImage;  
runAnalysisInputs.onsetsByCategory = onsetsByCategory;



if nargin == 0 || (nargin == 2 && strcmp(varargin{1},'paramBuilder')) || (nargin == 2 && strcmp(varargin{1},'preprocessed'))
  runAnalyses(runAnalysisInputs);
else
  feval(varargin{end},runAnalysisInputs);
end
end

