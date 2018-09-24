function [ spikesByEvent, spikesByCategory, lfpByEvent, lfpByCategory, categoryList, eventIDs ] = processRun( varargin )
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
%       - empty (default assignments: buildAnalysisParamFile, runAnalyses)
%       - 'paramBuilder', paramBuilderFunctionName
%       - 'analyzer', analyzerFunctionName
%       - 'paramBuilder', paramBuilderFunctionName, 'analyzer', analyzerFunctionName
%       Note: functionNames are strings, and do not include the trailing '.m'; see 'feval' docs
%   Notes:
%   Depends:
%   - contents of 'dependencies' folder (details coming)
%   - R2016a (or later) if joint psth-evoked potential plots desired
%   - Signal Processing Toolbox (for dpss taper calculation, LFP filters)

addpath(genpath('dependencies/genpath_exclude'));
addpath(genpath_exclude('dependencies',{'*mvgc_v1.0'})); %note: use this to exclude libraries that overwrite matlab builtin namespaces, until they're needed
addpath('buildAnalysisParamFileLib');
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

  %%%%%%%%%%%%%%%%%
  analogInData = preprocessAnalogIn(analogInFilename, analogInParams); 
  [spikesByChannel, taskTriggers, channelUnitNames] = preprocessSpikes(spikeFilename, ephysParams);
  lineNoiseTriggers = preprocessStrobe(lineNoiseTriggerFilename, lineNoiseTriggerParams);
  lfpData = preprocessLFP(lfpFilename, ephysParams, lineNoiseTriggers);
  diodeTriggers = preprocessStrobe(photodiodeFilename, photodiodeParams);
  [ taskData, stimTiming ] = preprocessLogFile(taskFilename, taskTriggers, diodeTriggers, stimSyncParams); % load visual stimulus data and transform its timestamps to ephys clock reference
%   [ pre, post ] = spikeBackground( spikesByChannel, taskData, spikeChannels, params )
  analogInData = preprocessEyeSignals(analogInData,taskData,eyeCalParams);
  analogInData = preprocessAccelSignals(analogInData, accelParams); 
  taskDataAll = taskData;
  taskData = excludeTrials( taskData, excludeStimParams); %exclude trials for lost fixation etc. 
  
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

  % todo: need defense against image with onset but no offset? 
  % todo: add similar defense for rf map locations here?
  if ~taskData.RFmap
    disp('No presentations of the following images survived exclusion:');
    disp(eventIDs(eventsNotObserved == 1));
  end
  onsetsByEvent = onsetsByEvent(eventsNotObserved == 0);
  offsetsByEvent = offsetsByEvent(eventsNotObserved == 0);
  trialIDsByEvent = trialIDsByEvent(eventsNotObserved == 0);
  eventIDs = eventIDs(eventsNotObserved == 0);
  eventLabels = eventLabels(eventsNotObserved == 0);
  eventCategories = eventCategories(eventsNotObserved == 0);
  
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
  
  onsetsByCategory = onsetsByCategory(catsNotObserved == 0);
  offsetsByCategory = offsetsByCategory(catsNotObserved == 0);
  categoryList = categoryList(catsNotObserved == 0);
  trialIDsByCategory = trialIDsByCategory(catsNotObserved == 0);
  
  if taskData.RFmap
    jumpsByImage = cell(length(eventIDs),1);
    for i = 1:length(eventIDs)
      jumpsByImage{i} = taskData.stimJumps(strcmp(taskData.taskEventIDs,eventIDs{i}),:);
    end
  else
    jumpsByImage = [];
  end
  %align spikes by trial, and sort by image and category
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
  %  align LFP data  by trial, sort by image and category, and possibly remove DC and linear components
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
      'categoryList', 'eventLabels', 'jumpsByImage', 'spikesByEvent', 'psthEmptyByEvent', 'spikesByCategory', 'psthEmptyByCategory',...
      'spikesByEventForTF', 'spikesByCategoryForTF', 'lfpByEvent', 'lfpByCategory', 'analogInByEvent','analogInByCategory','channelUnitNames', ...
      'stimTiming', 'eventCategories', 'onsetsByEvent', 'onsetsByCategory');
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
runAnalysisInputs.eventLabels = eventLabels; 
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
runAnalysisInputs.trialIDsByEvent = trialIDsByEvent;
runAnalysisInputs.onsetsByCategory = onsetsByCategory;
runAnalysisInputs.trialIDsByCategory = trialIDsByCategory;


if nargin == 0 || (nargin == 2 && strcmp(varargin{1},'paramBuilder')) || (nargin == 2 && strcmp(varargin{1},'preprocessed'))
  runAnalyses(runAnalysisInputs);
else
  feval(varargin{end},runAnalysisInputs);
end
end

