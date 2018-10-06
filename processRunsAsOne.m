function [ runAnalysisInputs ] = processRunsAsOne(runList, varargin)
%processRunsAsOne merges the data structures that processRun passes to
%runAnalyses from multiple runs that should be treated as one.
%   - assumes that all stimuli are described in one stimulusParamFile
%   - adds'RUN_BREAK' event at beginning of second and following runs
%   - assumes that unit numbering is constant; unsorts any unit numbers
%     that appear in some but not all trials.
%   - runList: runs to analyze, as cell array of runNum strings, e.g. {'002';'003'}
%   - varargin: name-value pairs:
%               'analysisParamSource': full path to a mat file, or handle to analysisParamFileBuilder;
%                                       default is to use  buildAnalysisParamFile for all params except dateSubj and runNum
%               'unitsIdentical': logical, if true, don't check unit assignment consistency
%               'savePreprocessed': logical, if true, save the preprocessed combined output data structures
%   Notes and todo:
%     - currently overwrites the analysisParams file for each constituent run
%     - cannot currently combine an RF map run with a non-RF map run

for argPair_i=1:length(varargin)/2
  argName = varargin{1+2*(argPair_i-1)};
  argVal = varargin{2+2*(argPair_i-1)};
  if strcmp(argName, 'analysisParamSource')
    if isstring(argVal)
      analysisParamFilename = argVal;
    else
      analysisParamFilename = argVal();
    end
  elseif strcmp(argName, 'unitsIdentical')
    unitsIdentical = argVal;
  end
end
if ~exist('analysisParamFilename', 'var')
    analysisParamFilename = buildAnalysisParamFile();
end
if ~exist('unitsIdentical','var')
  unitsIdentical = 0;
end
load(analysisParamFilename);

% resolve unit mapping:
%   - assume all units that appear with same number are the same
%   - unsort units that appear in some but not all runs
if ~unitsIdentical
  unitsToInclude = cell(length(ephysParams.spikeChannels),1);
  for run_i = 1:length(runList)
    spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runList{run_i});
    NEV = openNEV(spikeFilename,'read','nosave','noparse'); %note: add param 'report' for verbose
    for channel_i = 1:length(ephysParams.spikeChannels)
      unitNumbers = unique(NEV.Data.Spikes.Unit(NEV.Data.Spikes.Electrode == ephysParams.spikeChannels(channel_i)));
      if channel_i == 1
        unitsToInclude{channel_i} = unitNumbers;
      else
        unitsToInclude{channel_i} = intersect(unitsToInclude{channel_i},unitNumbers);
        ephysParams.unitsToUnsort{channel_i} = union(ephysParams.unitsToUnsort{channel_i}, setdiff(unitsToInclude{channel_i}, unitNumbers));
      end
    end  
  end
end

%%% set paths and directories for combined run%%%
runNum = '';
for run_i = 1:length(runList)
  runNum = strcat(runNum, runList{run_i});
end
finalOutDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
finalAnalysisParamFilename = strcat(finalOutDir,analysisParamFilenameStem);
finalPreprocessedDataFilename = strcat(finalOutDir,preprocessedDataFilenameStem);                     
%
if ~exist(finalOutDir,'dir')
  mkdir(finalOutDir);
end
outDir = finalOutDir;                                       %#ok
analysisParamFilename = finalAnalysisParamFilename;         %#ok
preprocessedDataFilename = finalPreprocessedDataFilename;   
save(finalAnalysisParamFilename);

for run_i = 1:length(runList)
  runNum = runList{run_i};
  analysisParamFilename = strcat(finalOutDir,'run',runNum,analysisParamFilenameStem);
  save(analysisParamFilename);
  runAnalysisInputsTmp = processRun('paramFile',analysisParamFilename,'analyzer', @(x) 1,'keepItemsNotPresented',1);
  if run_i == 1
    runAnalysisInputs = runAnalysisInputsTmp;
    runAnalysisInputs.analysisParamFilename = finalAnalysisParamFilename;
    accumulatedTimeOffset = runAnalysisInputs.excludeStimParams.ephysDuration;
    accumulatedTrialOffset = length(runAnalysisInputs.taskData.taskEventStartTimes);
    continue
  end
  % first, defend against spurrious computation of timeseries across run boundaries
  dataReqPre = psthParams.psthPre + max(3*psthParams.smoothingWidth, tfParams.movingWin(1)/2);
  if exist('spikeBackgroundParams','var') && spikeBackgroundParams.trialWise
    dataReqPre = max(dataReqPre, spikeBackgroundParams.compWin(1));
  end
  if exist('lfpBackgroundParams','var') && lfpBackgroundParams.trialWise
    dataReqPre = max(dataReqPre, lfpBackgroundParams.compWin(1));
  end
  dataReqPost = psthParams.psthPost + max(3*psthParams.smoothingWidth, tfParams.movingWin(1)/2);

  if ~(isfield(runAnalysisInputsTmp.excludeStimParams,'ephysDataPre') && ...
      runAnalysisInputsTmp.excludeStimParams.ephysDataPre > dataReqPre && ...
      isfield(runAnalysisInputsTmp.excludeStimParams,'ephysDataPost') && ...
      runAnalysisInputsTmp.excludeStimParams.ephysDataPost > dataReqPost)
    assert(false,'processRunAsOne currently requires that trials be excluded per-run for insufficient ephys data pre- and post');
  end
  runAnalysisInputs.lfpData = cat(2,runAnalysisInputs.lfpData,runAnalysisInputsTmp.lfpData); 
  runAnalysisInputs.analogInData = cat(2,runAnalysisInputs.analogInData,runAnalysisInputsTmp.analogInData); 
  % merge taskData
  runAnalysisInputs.taskData.taskEventStartTimes = vertcat(runAnalysisInputs.taskData.taskEventStartTimes,accumulatedTimeOffset,...
    runAnalysisInputsTmp.taskData.taskEventStartTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskData.taskEventEndTimes = vertcat(runAnalysisInputs.taskData.taskEventEndTimes,accumulatedTimeOffset+1,...
    runAnalysisInputsTmp.taskData.taskEventEndTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskData.fixationInTimes = vertcat(runAnalysisInputs.taskData.fixationInTimes,...
    runAnalysisInputsTmp.taskData.fixationInTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskData.fixationOutTimes = vertcat(runAnalysisInputs.taskData.fixationOutTimes,...
    runAnalysisInputsTmp.taskData.fixationOutTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskData.juiceOnTimes = vertcat(runAnalysisInputs.taskData.juiceOnTimes,...
    runAnalysisInputsTmp.taskData.juiceOnTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskData.juiceOffTimes = vertcat(runAnalysisInputs.taskData.juiceOffTimes,...
    runAnalysisInputsTmp.taskData.juiceOffTimes + accumulatedTimeOffset);
  if runAnalysisInputsTmp.taskData.RFmap
    assert(isfield(runAnalysisInputs.taskData.stimJumps),'combination of RF map and non-RF-map runs not currently implemented');
    runAnalysisInputs.taskData.stimJumps = vertcat(runAnalysisInputs.taskData.stimJumps,[0,0],runAnalysisInputsTmp.taskData.stimJumps);
    runAnalysisInputs.taskData.gridPointsDegX = sort(union(runAnalysisInputs.taskData.gridPointsDegX,runAnalysisInputsTmp.taskData.gridPointsDegX),'ascend');
    runAnalysisInputs.taskData.gridPointsDegY = sort(union(runAnalysisInputs.taskData.gridPointsDegY,runAnalysisInputsTmp.taskData.gridPointsDegY),'ascend');
  end
  runAnalysisInputs.taskData.taskEventIDs = vertcat(runAnalysisInputs.taskData.taskEventIDs,'RUN_BREAK',...
    runAnalysisInputsTmp.taskData.taskEventIDs);
  
  runAnalysisInputs.taskData.stimFramesLost = vertcat(runAnalysisInputs.taskData.stimFramesLost,0,runAnalysisInputsTmp.taskData.stimFramesLost);
  runAnalysisInputs.taskData.eventTimeAdjustments = vertcat(runAnalysisInputs.taskData.eventTimeAdjustments,0,...
    runAnalysisInputsTmp.taskData.eventTimeAdjustments);
  % merge taskDataAll
  runAnalysisInputs.taskDataAll.taskEventStartTimes = vertcat(runAnalysisInputs.taskDataAll.taskEventStartTimes,accumulatedTimeOffset,...
    runAnalysisInputsTmp.taskDataAll.taskEventStartTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskDataAll.taskEventEndTimes = vertcat(runAnalysisInputs.taskDataAll.taskEventEndTimes,accumulatedTimeOffset+1,...
    runAnalysisInputsTmp.taskDataAll.taskEventEndTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskDataAll.fixationInTimes = vertcat(runAnalysisInputs.taskDataAll.fixationInTimes,...
    runAnalysisInputsTmp.taskDataAll.fixationInTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskDataAll.fixationOutTimes = vertcat(runAnalysisInputs.taskDataAll.fixationOutTimes,...
    runAnalysisInputsTmp.taskDataAll.fixationOutTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskDataAll.juiceOnTimes = vertcat(runAnalysisInputs.taskDataAll.juiceOnTimes,...
    runAnalysisInputsTmp.taskDataAll.juiceOnTimes + accumulatedTimeOffset);
  runAnalysisInputs.taskDataAll.juiceOffTimes = vertcat(runAnalysisInputs.taskDataAll.juiceOffTimes,...
    runAnalysisInputsTmp.taskDataAll.juiceOffTimes + accumulatedTimeOffset);
  if runAnalysisInputsTmp.taskData.RFmap
    assert(isfield(runAnalysisInputs.taskDataAll.stimJumps),'combination of RF map and non-RF-map runs not currently implemented');
    runAnalysisInputs.taskDataAll.stimJumps = vertcat(runAnalysisInputs.taskDataAll.stimJumps,[0,0],runAnalysisInputsTmp.taskDataAll.stimJumps);
    runAnalysisInputs.taskDataAll.gridPointsDegX = sort(union(runAnalysisInputs.taskDataAll.gridPointsDegX,runAnalysisInputsTmp.taskDataAll.gridPointsDegX),'ascend');
    runAnalysisInputs.taskDataAll.gridPointsDegY = sort(union(runAnalysisInputs.taskDataAll.gridPointsDegY,runAnalysisInputsTmp.taskDataAll.gridPointsDegY),'ascend');
  end
  runAnalysisInputs.taskDataAll.taskEventIDs = vertcat(runAnalysisInputs.taskDataAll.taskEventIDs,'RUN_BREAK',...
    runAnalysisInputsTmp.taskDataAll.taskEventIDs);
  runAnalysisInputs.taskDataAll.stimFramesLost = vertcat(runAnalysisInputs.taskDataAll.stimFramesLost,0,runAnalysisInputsTmp.taskDataAll.stimFramesLost);
  runAnalysisInputs.taskDataAll.eventTimeAdjustments = vertcat(runAnalysisInputs.taskDataAll.eventTimeAdjustments,0,...
    runAnalysisInputsTmp.taskDataAll.eventTimeAdjustments);

  % merge 'byEvent' data structures for spikes, lfps, analogIns, jumps, onsets, and trialIDs
  for event_i = 1:length(runAnalysisInputs.eventLabels)
    for channel_i = 1:length(runAnalysisInputsTmp.spikesByEvent{event_i})
      for unit_i = 1:length(runAnalysisInputsTmp.spikesByEvent{1}{channel_i})
        runAnalysisInputs.spikesByEvent{event_i}{channel_i}{unit_i} = vertcat(runAnalysisInputs.spikesByEvent{event_i}{channel_i}{unit_i},...
          runAnalysisInputsTmp.spikesByEvent{event_i}{channel_i}{unit_i}); 
        runAnalysisInputs.psthEmptyByEvent{event_i}{channel_i}{unit_i} = vertcat(runAnalysisInputs.psthEmptyByEvent{event_i}{channel_i}{unit_i},...
          runAnalysisInputsTmp.psthEmptyByEvent{event_i}{channel_i}{unit_i});
        runAnalysisInputs.spikesByEventForTF{event_i}{channel_i}{unit_i} = vertcat(runAnalysisInputs.spikesByEventForTF{event_i}{channel_i}{unit_i},...
          runAnalysisInputsTmp.spikesByEventForTF{event_i}{channel_i}{unit_i});
      end      
    end
    runAnalysisInputs.lfpByEvent{event_i} = cat(3, runAnalysisInputs.lfpByEvent{event_i}, runAnalysisInputsTmp.lfpByEvent{event_i});
    runAnalysisInputs.analogInByEvent{event_i} = cat(3, runAnalysisInputs.analogInByEvent{event_i}, runAnalysisInputsTmp.analogInByEvent{event_i});
    if runAnalysisInputsTmp.taskData.RFmap
      runAnalysisInputs.jumpsByImage{event_i} = vertcat(runAnalysisInputs.jumpsByImage{event_i}, runAnalysisInputsTmp.jumpsByImage{event_i});
    end
    runAnalysisInputs.onsetsByEvent{event_i} = vertcat(runAnalysisInputs.onsetsByEvent{event_i},runAnalysisInputsTmp.onsetsByEvent{event_i});
    runAnalysisInputs.trialIDsByEvent{event_i} = vertcat(runAnalysisInputs.trialIDsByEvent{event_i},...
      runAnalysisInputsTmp.trialIDsByEvent{event_i}+accumulatedTrialOffset);
  end
  % merge 'byCategory' data structures for spikes, lfps, analogIns, jumps, onsets, and trialIDs
  for cat_i = 1:length(runAnalysisInputs.categoryList)
    for channel_i = 1:length(runAnalysisInputsTmp.spikesByCategory{cat_i})
      for unit_i = 1:length(runAnalysisInputsTmp.spikesByCategory{1}{channel_i})
        runAnalysisInputs.spikesByCategory{cat_i}{channel_i}{unit_i} = vertcat(runAnalysisInputs.spikesByCategory{cat_i}{channel_i}{unit_i},...
          runAnalysisInputsTmp.spikesByCategory{cat_i}{channel_i}{unit_i}); 
        runAnalysisInputs.psthEmptyByCategory{cat_i}{channel_i}{unit_i} = vertcat(runAnalysisInputs.psthEmptyByCategory{cat_i}{channel_i}{unit_i},...
          runAnalysisInputsTmp.psthEmptyByCategory{cat_i}{channel_i}{unit_i});
        runAnalysisInputs.spikesByCategoryForTF{cat_i}{channel_i}{unit_i} = vertcat(runAnalysisInputs.spikesByCategoryForTF{cat_i}{channel_i}{unit_i},...
          runAnalysisInputsTmp.spikesByCategoryForTF{cat_i}{channel_i}{unit_i});
      end
    end
    runAnalysisInputs.lfpByCategory{cat_i} = cat(3, runAnalysisInputs.lfpByCategory{cat_i}, runAnalysisInputsTmp.lfpByCategory{cat_i});
    runAnalysisInputs.analogInByCategory{cat_i} = cat(3, runAnalysisInputs.analogInByCategory{cat_i}, runAnalysisInputsTmp.analogInByCategory{cat_i});
    
    runAnalysisInputs.onsetsByCategory{cat_i} = vertcat(runAnalysisInputs.onsetsByCategory{cat_i},runAnalysisInputsTmp.onsetsByCategory{cat_i});
    runAnalysisInputs.trialIDsByCategory{cat_i} = vertcat(runAnalysisInputs.trialIDsByCategory{cat_i},...
      runAnalysisInputsTmp.trialIDsByCategory{cat_i} + accumulatedTrialOffset);
  end
  runAnalysisInputs.stimTiming.shortest = min(runAnalysisInputs.stimTiming.shortest,runAnalysisInputsTmp.stimTiming.shortest);
  runAnalysisInputs.stimTiming.longest = min(runAnalysisInputs.stimTiming.longest,runAnalysisInputsTmp.stimTiming.longest);

  accumulatedTimeOffset = accumulatedTimeOffset + runAnalysisInputsTmp.excludeStimParams.ephysDuration;
  accumulatedTrialOffset = accumulatedTrialOffset + length(runAnalysisInputsTmp.taskData.taskEventStartTimes);
end

%remove lines for events not observed from data structures
eventsNotObserved = zeros(length(runAnalysisInputs.eventLabels),1);
for event_i = 1:length(runAnalysisInputs.eventLabels)
  eventsNotObserved(event_i) = isempty(runAnalysisInputs.onsetsByEvent{event_i});
end
runAnalysisInputs.onsetsByEvent = runAnalysisInputs.onsetsByEvent(eventsNotObserved == 0);
runAnalysisInputs.offsetsByEvent = runAnalysisInputs.offsetsByEvent(eventsNotObserved == 0);
runAnalysisInputs.trialIDsByEvent = runAnalysisInputs.trialIDsByEvent(eventsNotObserved == 0);
runAnalysisInputs.eventLabels = runAnalysisInputs.eventLabels(eventsNotObserved == 0);
runAnalysisInputs.eventCategories = runAnalysisInputs.eventCategories(eventsNotObserved == 0);
runAnalysisInputs.spikesByEvent = runAnalysisInputs.spikesByEvent(eventsNotObserved == 0);
runAnalysisInputs.psthEmptyByEvent = runAnalysisInputs.psthEmptyByEvent(eventsNotObserved == 0);
runAnalysisInputs.spikesByEventForTF = runAnalysisInputs.spikesByEventForTF(eventsNotObserved == 0);
runAnalysisInputs.lfpByEvent = runAnalysisInputs.lfpByEvent(eventsNotObserved == 0);
runAnalysisInputs.analogInByEvent = runAnalysisInputs.analogInByEvent(eventsNotObserved == 0);
if runAnalysisInputsTmp.taskData.RFmap
  runAnalysisInputs.jumpsByImage = runAnalysisInputs.jumpsByImage(eventsNotObserved == 0);
end
%remove lines for categories not observed from data structures
catsNotObserved = zeros(length(runAnalysisInputs.categoryList),1);
for cat_i = 1:length(runAnalysisInputs.categoryList)
  catsNotObserved(cat_i) = isempty(runAnalysisInputs.onsetsByCategory{cat_i});
end

runAnalysisInputs.onsetsByCategory = runAnalysisInputs.onsetsByCategory(catsNotObserved == 0);
runAnalysisInputs.offsetsByCategory = runAnalysisInputs.offsetsByCategory(catsNotObserved == 0);
runAnalysisInputs.categoryList = runAnalysisInputs.categoryList(catsNotObserved == 0);
runAnalysisInputs.trialIDsByCategory = runAnalysisInputs.trialIDsByCategory(catsNotObserved == 0);
runAnalysisInputs.spikesByCategory = runAnalysisInputs.spikesByCategory(catsNotObserved == 0);
runAnalysisInputs.psthEmptyByCategory = runAnalysisInputs.psthEmptyByCategory(catsNotObserved == 0);
runAnalysisInputs.spikesByCategoryForTF = runAnalysisInputs.spikesByCategoryForTF(catsNotObserved == 0);
runAnalysisInputs.lfpByCategory = runAnalysisInputs.lfpByCategory(catsNotObserved == 0);
runAnalysisInputs.analogInByCategory = runAnalysisInputs.analogInByCategory(catsNotObserved == 0);

if savePreprocessed
  analysisParamFilename = runAnalysisInputs. analysisParamFilename; %#ok
  spikesByChannel = runAnalysisInputs. spikesByChannel;             %#ok
  lfpData = runAnalysisInputs. lfpData;                             %#ok
  analogInData = runAnalysisInputs. analogInData;                   %#ok
  taskData = runAnalysisInputs. taskData;                           %#ok
  taskDataAll = runAnalysisInputs. taskDataAll;                     %#ok
  psthImDur = runAnalysisInputs. psthImDur;                         %#ok
  preAlign = runAnalysisInputs. preAlign;                           %#ok
  postAlign = runAnalysisInputs. postAlign;                         %#ok
  categoryList = runAnalysisInputs. categoryList;                   %#ok
  eventLabels = runAnalysisInputs. eventLabels;                     %#ok
  jumpsByImage = runAnalysisInputs. jumpsByImage;                   %#ok
  spikesByEvent = runAnalysisInputs. spikesByEvent;                 %#ok
  psthEmptyByEvent = runAnalysisInputs. psthEmptyByEvent;           %#ok
  spikesByCategory = runAnalysisInputs. spikesByCategory;           %#ok
  psthEmptyByCategory = runAnalysisInputs. psthEmptyByCategory;     %#ok
  spikesByEventForTF = runAnalysisInputs. spikesByEventForTF;       %#ok
  spikesByCategoryForTF = runAnalysisInputs. spikesByCategoryForTF; %#ok
  lfpByEvent = runAnalysisInputs. lfpByEvent;                       %#ok
  lfpByCategory = runAnalysisInputs. lfpByCategory;                 %#ok
  analogInByEvent = runAnalysisInputs. analogInByEvent;             %#ok
  analogInByCategory = runAnalysisInputs. analogInByCategory;       %#ok
  channelUnitNames = runAnalysisInputs. channelUnitNames;           %#ok
  stimTiming = runAnalysisInputs. stimTiming;                       %#ok
  eventCategories = runAnalysisInputs. eventCategories;             %#ok
  onsetsByEvent = runAnalysisInputs. onsetsByEvent;                 %#ok
  onsetsByCategory = runAnalysisInputs. onsetsByCategory;           %#ok
  trialIDsByEvent = runAnalysisInputs. trialIDsByEvent;             %#ok
  trialIDsByCategory = runAnalysisInputs. trialIDsByCategory;       %#ok
  save(finalPreprocessedDataFilename,'analysisParamFilename', 'spikesByChannel', 'lfpData', 'analogInData', 'taskData', 'taskDataAll', 'psthImDur', 'preAlign', 'postAlign',...
    'categoryList', 'eventLabels', 'jumpsByImage', 'spikesByEvent', 'psthEmptyByEvent', 'spikesByCategory', 'psthEmptyByCategory',...
    'spikesByEventForTF', 'spikesByCategoryForTF', 'lfpByEvent', 'lfpByCategory', 'analogInByEvent','analogInByCategory','channelUnitNames', ...
    'stimTiming', 'eventCategories', 'onsetsByEvent', 'onsetsByCategory','trialIDsByEvent','trialIDsByCategory');
end

if ~exist('analyzer','var')
  runAnalyses(runAnalysisInputs);
else
  feval(analyzer,runAnalysisInputs);
end
end

