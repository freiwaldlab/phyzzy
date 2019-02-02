function [runAnalysisInputs, analysisOutFilename] = processRun( varargin )
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
%   - varargin consists of name-value argument pairs:
%       - 'paramBuilder', paramBuilderFunctionName
%       - 'paramFile', paramFilename (must be .mat)
%       - 'analyzer', analyzerFunctionName
%       - 'preprocessed', preProcessedDataFilename (must be .mat), or '-d'
%             for default, in which case preprocessedDataFilename will be
%             extracted from the dateSubj+runNum in buildAnalysisParamFile,
%             paramBuilder, or paramFile. Unless paramFile is supplied,
%             will overwrite the calcSwitch, plotSwitch, and analysisGroups
%             fields of the saved analysisParamFile, based on their current
%             values in the param file builder. 
%       - 'preprocessedCC', preProcessedDataFilename (must be .mat), or '-d'
%             identical behavior to preprocessed, except that no variables
%             in analysisParamFile will be updated or overwritten. This is
%             useful to reproduce an analysis, or to run identically specified analyses
%             with a different or modified analyzer.
%       - 'keepItemsNotPresented', logical, default 0. If 0, events that
%           appear in the stimulus param files, but that have no valid trials,
%           do not appear in output data structures. 
%       Note: functionNames are strings, and do not include the trailing '.m'; see 'feval' docs
%   Notes:
%   Depends:
%   - contents of 'dependencies' folder (details coming)
%   - R2016a (or later) if joint psth-evoked potential plots desired
%   - Signal Processing Toolbox (for dpss taper calculation, LFP filters)

addpath(genpath('dependencies'));
rmpath(genpath('dependencies/mvgc_v1.0')); %note: Not sure if this is appropriate replacement for genpath_exclude. previous line caused issues in parallel runs.
addpath('buildAnalysisParamFileLib');

%% load analysis parameters

%Default values
usePreprocessed = 0;
preprocessedCC = 0;
defaultPreprocessed = 0;
keepItemsNotPresented = 0;

%Parse inputs
assert(mod(nargin,2) == 0, 'processRun takes arguments as name-value pairs; odd number of arguments provided');
for argPair_i=1:nargin/2
  argName = varargin{1+2*(argPair_i-1)};
  argVal = varargin{2+2*(argPair_i-1)};
  switch argName
    case 'paramBuilder'
      paramBuilder = argVal;
    case 'analyzer'
      analyzer = argVal;
    case 'paramFile'
      analysisParamFilename = argVal;
    case 'keepItemsNotPresented'
      keepItemsNotPresented = argVal;
    case 'preprocessed'
      if strcmp(argVal,'-d')
        defaultPreprocessed = 1;
      else
        preprocessedDataFilename = argVal;
      end
      usePreprocessed = 1;
    case 'preprocessedCC'
      if strcmp(argVal,'-d')
        defaultPreprocessed = 1;
      else
        preprocessedDataFilename = argVal;
      end
      usePreprocessed = 1;
      preprocessedCC = 1;
    otherwise
      error('Invalid field %s provided to processRun',argName);
  end
end

%% ParamBuilder and preprocessed file handling

if usePreprocessed
  if defaultPreprocessed
    if preprocessedCC 
      if ~exist('analysisParamFilename','var')
        if exist('paramBuilder','var')
          analysisParamFilename = feval(paramBuilder,'noSave');
        else
          analysisParamFilename = buildAnalysisParamFile('noSave');
        end
      end
    else
      if ~exist('analysisParamFilename','var')
        if exist('paramBuilder','var')
          analysisParamFilename = feval(paramBuilder,'saveNoPreprocParams');
        else
          analysisParamFilename = buildAnalysisParamFile('saveNoPreprocParams');
        end
      end
    end
    load(analysisParamFilename,'preprocessedDataFilename');
  end
  load(preprocessedDataFilename);
else
  if ~exist('analysisParamFilename','var')
    if exist('paramBuilder','var')
      analysisParamFilename = feval(paramBuilder);
    else
      analysisParamFilename = buildAnalysisParamFile();
    end
  end
end

if ~usePreprocessed
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

  %% Preprocess inputs
  
  analogInData =                                          preprocessAnalogIn(analogInFilename, analogInParams); 
  [spikesByChannel, taskTriggers, channelUnitNames] =     preprocessSpikes(spikeFilename, ephysParams);
  lineNoiseTriggers =                                     preprocessStrobe(lineNoiseTriggerFilename, lineNoiseTriggerParams);
  lfpData =                                               preprocessLFP(lfpFilename, ephysParams, lineNoiseTriggers);
  diodeTriggers =                                         preprocessStrobe(photodiodeFilename, photodiodeParams);
  [taskData, stimTiming ] =                               preprocessLogFile(taskFilename, taskTriggers, diodeTriggers, stimSyncParams); % load visual stimulus data and transform its timestamps to ephys clock reference
%   [pre, post] =                                           spikeBackground( spikesByChannel, taskData, spikeChannels, params )
  analogInData =                                          preprocessEyeSignals(analogInData,taskData,eyeCalParams);
  analogInData =                                          preprocessAccelSignals(analogInData, accelParams); 
  
  
  % determine the duration of ephys data collection, in ms
  if ~isempty(lfpData)
    ephysDuration = size(lfpData,2); 
  elseif ~isempty(analogInData)
    ephysDuration = size(analogInData,2); 
  elseif ~isempty(spikesByChannel)
    ephysDuration = 0;
    for channel_i = 1:length(spikesByChannel)
      ephysDuration = max(ephysDuration, max(spikesByChannel{channel_i}.times)); 
    end
    lastTaskTriggerTime = 1000*max(taskTriggers.TimeStampSec); 
    if lastTaskTriggerTime > ephysDuration
      ephysDuration = lastTaskTriggerTime;
      warning('Using final digital trigger time to define end of ephys data collection, because it followed the final spike, and because no LFP or analog in data available. May bias results');
    else
      warning('Using final spike time to define end of ephys data collection, because it followed the final digital trigger, and because no LFP or analog in data available. May bias results');
    end
  end
  excludeStimParams.ephysDuration = ephysDuration;
  
  taskDataAll = taskData;
  taskData = excludeTrials(taskData, excludeStimParams); %exclude trials for lost fixation etc. 
  

%% find Stimulus length

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

%% Sort trials by image and image category, align data appropriately.

  tmp = load(stimParamsFilename); %loads variables paramArray, categoryLabels,eventLabels
  eventCategories = tmp.paramArray;
  categoryList = tmp.categoryLabels;
  try
    eventLabels = tmp.eventLabels;
  catch
    eventLabels = tmp.pictureLabels; %implements back-compatibility with previous variable names
  end
  eventIDs = cell(size(eventCategories,1),1);
  for event_i = 1:length(eventCategories)
    eventIDs{event_i} = eventCategories{event_i}{1}; %per the stimParamFile spec, this is the event ID
  end
  onsetsByEvent = cell(length(eventIDs),1);
  offsetsByEvent = cell(length(eventIDs),1);
  eventsNotObserved = zeros(length(eventIDs),1);
  trialIDsByEvent = cell(length(eventIDs),1);
  for event_i = 1:length(eventIDs)
    onsetsByEvent{event_i} = taskData.taskEventStartTimes(strcmp(taskData.taskEventIDs,eventIDs{event_i}));
    offsetsByEvent{event_i} = taskData.taskEventEndTimes(strcmp(taskData.taskEventIDs,eventIDs{event_i}));
    eventsNotObserved(event_i) = isempty(onsetsByEvent{event_i});
    trialIDsByEvent{event_i} = find(strcmp(taskData.taskEventIDs,eventIDs{event_i}));
  end

  assert(~(sum(eventsNotObserved) == length(eventsNotObserved)),'No Events observed. Confirm correct stimulus Param file is in use.');
  if ~(sum(eventsNotObserved == 0) == length(unique(taskData.taskEventIDs)))
    uni = unique(taskData.taskEventIDs);
    for ii = 1:length(uni)
      chck = sum(strcmp(eventIDs,uni(ii)));
      if ~chck
        fprintf('Event %s missing from stimParamFile \n', uni{ii});
      end
    end
    error('Not all Events observed are represented in stimParamFile.')
  end
  % todo: need defense against image with onset but no offset? 
  % todo: add similar defense for rf map locations here?
  if ~taskData.RFmap
    disp('No presentations of the following images survived exclusion:');
    %disp(eventIDs(eventsNotObserved == 1));
  end
  if ~keepItemsNotPresented
    onsetsByEvent = onsetsByEvent(eventsNotObserved == 0);
    offsetsByEvent = offsetsByEvent(eventsNotObserved == 0);
    trialIDsByEvent = trialIDsByEvent(eventsNotObserved == 0);
    eventIDs = eventIDs(eventsNotObserved == 0);
    eventLabels = eventLabels(eventsNotObserved == 0);
    eventCategories = eventCategories(eventsNotObserved == 0);
  end
  
  onsetsByCategory = cell(length(categoryList),1);
  offsetsByCategory = cell(length(categoryList),1);
  catsNotObserved = zeros(length(categoryList),1);
  trialIDsByCategory = cell(length(categoryList),1);
  for cat_i = 1:length(categoryList)
    catOnsets = [];
    catOffsets = [];
    catTrialIDs = [];
    for event_i = 1:length(eventIDs)
      if any(strcmp(eventCategories{event_i},categoryList{cat_i}))
        catOnsets = vertcat(catOnsets,onsetsByEvent{event_i});
        catOffsets = vertcat(catOffsets,offsetsByEvent{event_i});
        catTrialIDs = vertcat(catTrialIDs,trialIDsByEvent{event_i});
      end
    end
    onsetsByCategory{cat_i} = catOnsets;
    offsetsByCategory{cat_i} = catOffsets;
    catsNotObserved(cat_i) = isempty(onsetsByCategory{cat_i});
    trialIDsByCategory{cat_i} = catTrialIDs;
    Output.DEBUG('numel cat onsets');
    Output.DEBUG(numel(catOnsets));
  end
  
  if ~keepItemsNotPresented
    onsetsByCategory = onsetsByCategory(catsNotObserved == 0);
    offsetsByCategory = offsetsByCategory(catsNotObserved == 0);
    categoryList = categoryList(catsNotObserved == 0);
    trialIDsByCategory = trialIDsByCategory(catsNotObserved == 0);
  end
  
  if taskData.RFmap
    jumpsByImage = cell(length(eventIDs),1);
    for i = 1:length(eventIDs)
      jumpsByImage{i} = taskData.stimJumps(strcmp(taskData.taskEventIDs,eventIDs{i}),:);
    end
  else
    jumpsByImage = [];
  end
  
  % align spikes by trial, and sort by image and category
  spikeAlignParamsToCoverMovingWin.preAlign = lfpAlignParams.msPreAlign; % need to include spikes throughout the window chronux will use to calculate spectra
  spikeAlignParamsToCoverMovingWin.postAlign = lfpAlignParams.msPostAlign;
  spikeAlignParamsToCoverMovingWin.refOffset = 0;  
  spikeAlignParams.refOffset = 0;
  [spikesByEvent, psthEmptyByEvent] = alignSpikes( spikesByChannel, onsetsByEvent, spikeChannels, spikeAlignParamsToCoverMovingWin);
  [spikesByCategory, psthEmptyByCategory] = alignSpikes( spikesByChannel, onsetsByCategory, spikeChannels, spikeAlignParamsToCoverMovingWin);

  % align spikes again, but this time reference time to lfp sample number (required for chronux TF, even spike-spike)
  spikeAlignParamsTF.preAlign = lfpAlignParams.msPreAlign;
  spikeAlignParamsTF.postAlign = lfpAlignParams.msPostAlign;
  spikeAlignParamsTF.refOffset = -lfpAlignParams.msPreAlign;
  spikesByEventForTF = alignSpikes( spikesByChannel, onsetsByEvent, spikeChannels, spikeAlignParamsTF );
  spikesByCategoryForTF = alignSpikes( spikesByChannel, onsetsByCategory, spikeChannels, spikeAlignParamsTF );
  % align LFP data  by trial, sort by image and category, and possibly remove DC and linear components
  lfpByEvent = alignLFP(lfpData, onsetsByEvent, lfpChannels, lfpAlignParams);
  lfpByCategory = alignLFP(lfpData, onsetsByCategory, lfpChannels, lfpAlignParams);
  analogInByEvent = alignAnalogIn(analogInData, onsetsByEvent, analogInChannels, lfpAlignParams);
  analogInByCategory = alignAnalogIn(analogInData, onsetsByCategory, analogInChannels, lfpAlignParams);
  
  for cat_i = 1:length(categoryList)
    Output.VERBOSE(categoryList{cat_i});
    Output.VERBOSE(size(lfpByCategory{cat_i}));
  end

  if savePreprocessed
    save(preprocessedDataFilename,'analysisParamFilename', 'spikesByChannel', 'lfpData', 'analogInData', 'taskData', 'taskDataAll', 'psthImDur', 'preAlign', 'postAlign',...
      'categoryList', 'eventLabels', 'eventIDs', 'jumpsByImage', 'spikesByEvent', 'psthEmptyByEvent', 'spikesByCategory', 'psthEmptyByCategory',...
      'spikesByEventForTF', 'spikesByCategoryForTF', 'lfpByEvent', 'lfpByCategory', 'analogInByEvent','analogInByCategory','channelUnitNames', ...
      'stimTiming', 'eventCategories', 'onsetsByEvent', 'offsetsByEvent', 'onsetsByCategory', 'offsetsByCategory', 'trialIDsByEvent','trialIDsByCategory', 'excludeStimParams');
  end
end

%% Package Analysis inputs, and run Analysis

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
runAnalysisInputs.eventLabels = eventLabels; 
runAnalysisInputs.eventIDs = eventIDs;
runAnalysisInputs.jumpsByImage = jumpsByImage; 
runAnalysisInputs.spikesByEvent = spikesByEvent; 
runAnalysisInputs.psthEmptyByEvent = psthEmptyByEvent;  
runAnalysisInputs.spikesByCategory = spikesByCategory; 
runAnalysisInputs.psthEmptyByCategory = psthEmptyByCategory; 
runAnalysisInputs.spikesByEventForTF = spikesByEventForTF;  
runAnalysisInputs.spikesByCategoryForTF = spikesByCategoryForTF;  
runAnalysisInputs.lfpByEvent = lfpByEvent;  
runAnalysisInputs.lfpByCategory = lfpByCategory;  
runAnalysisInputs.analogInByEvent = analogInByEvent; 
runAnalysisInputs.analogInByCategory = analogInByCategory;  
runAnalysisInputs.channelUnitNames = channelUnitNames;
runAnalysisInputs.stimTiming = stimTiming;  
runAnalysisInputs.eventCategories = eventCategories;  
runAnalysisInputs.onsetsByEvent = onsetsByEvent;
runAnalysisInputs.offsetsByEvent = offsetsByEvent; 
runAnalysisInputs.trialIDsByEvent = trialIDsByEvent;
runAnalysisInputs.onsetsByCategory = onsetsByCategory;
runAnalysisInputs.offsetsByCategory = offsetsByCategory;
runAnalysisInputs.trialIDsByCategory = trialIDsByCategory;
runAnalysisInputs.excludeStimParams = excludeStimParams;

if ~exist('analyzer','var')
  [analysisOutFilename] = runAnalyses(runAnalysisInputs);
else
  feval(analyzer,runAnalysisInputs);
end
end

