function [] = processRunsAsOne(runList, varargin)
%UNTITLED Summary of this function goes here
%   - hwo to handle across-run dependencies? like fixation/trial before?
%     - put in dummy blank period betweek? 
%       - could also put in 'runBreak' event --> right idea
%     - alternative is to compute everything first, then combine outputs
%       - issues: need to keep building new stuff for new analyses
%
%   - runList: runs to analyze, as cell array of runNum strings, e.g. {'002';'003'}
%   - varargin: name-value pairs:
%               'analysisParamSource': full path to a mat file, or handle to analysisParamFileBuilder;
%                   default is to use  buildAnalysisParamFile for all params except dateSubj and runNum
%               'unitsIdentical': logical, if true, don't check unit assignment consistency

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
    for channel_i = 1:length(ephysParams.spikeChannels)
      NEV = openNEV(spikeFilename,'read','nosave','noparse'); %note: add param 'report' for verbose
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
runNums = '';
for run_i = 1:length(runList)
  runNums = strcat(runNums, runList{run_i});
end
finalOutDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNums);
finalAnalysisParamFilename = strcat(finalOutDir,analysisParamFilenameStem);
finalPreprocessedDataFilename = strcat(finalOutDir,preprocessedDataFilenameStem);                     %#ok
%
if ~exist(finalOutDir,'dir')
  mkdir(finalOutDir);
end
save(finalAnalysisParamFilename);


for run_i = 1:length(runList)
  runNum = runList{run_i};
  tmpAnalysisParamFilename = save(strcat(finalOutDir,analysisParamFilenameStem,'run',runNum));
  runAnalysisInputsTmp = processRun('paramFile',tmpAnalysisParamFilename,'analyzer', @(x) 1);
  if run_i == 1
    runAnalysisInputs = runAnalysisInputsTmp;
    accumulatedTimeOffset = runAnalysisInputsTmp.excludeStimParams.ephysDuration;
    accumulatedTrialOffset = length(runAnalysisInputsTmp.taskData.taskEventStartTimes);
  else
    % first, defend against spurrious computation of timeseries across run boundaries
    % pre
    dataReqPre = psthParams.psthPre + max(3*psthParams.smoothingWidth, tfParams.movingWin(1)/2);
    if exist('spikeBackgroundParams','var') && spikeBackgroundParams.trialWise
      dataReqPre = max(dataReqPre, spikeBackgroundParams.compWin(1));
    end
    if exist('lfpBackgroundParams','var') && lfpBackgroundParams.trialWise
      dataReqPre = max(dataReqPre, lfpBackgroundParams.compWin(1));
    end
    % post
    dataReqPost = psthParams.psthPost + max(3*psthParams.smoothingWidth, tfParams.movingWin(1)/2);
  
    if ~(isfield(runAnalysisInputsTmp.excludeStimParams,'ephysDataPre') && ...
        runAnalysisInputsTmp.excludeStimParams,'ephysDataPre') > dataReqPre && ...
        isfield(runAnalysisInputsTmp.excludeStimParams,'ephysDataPost') && ...
        runAnalysisInputsTmp.excludeStimParams,'ephysDataPost') > dataReqPost)
      assert(0,'processRunAsOne currently requires that trials be excluded per-run for insufficient ephys data pre- and post');
%       to do: implement removing trials at the beginning and end of runs if it hasn't been handled per-run 
%       taskDataValid = runAnalysisInputsTmp.taskData.taskEventStartTimes;
%       taskDataValid.taskEventIDs = taskData.taskEventIDs(trialValid);
%       taskDataValid.stimJumps = taskData.stimJumps(trialValid,:);
%       taskDataValid.stimFramesLost = taskData.stimFramesLost(trialValid);
%       taskDataValid.taskEventStartTimes = taskData.taskEventStartTimes(trialValid);
%       taskDataValid.taskEventEndTimes = taskData.taskEventEndTimes(trialValid); 
    end
    runAnalysisInputs.lfpData = cat(2,runAnalysisInputs.lfpData,runAnalysisInputsTmp.lfpData); 
    runAnalysisInputs.analogInData = cat(2,runAnalysisInputs.analogInData,runAnalysisInputsTmp.analogInData); 
    
    runAnalysisInputs.taskData.taskEventStartTimes = runAnalysisInputsTmp.taskData.taskEventStartTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskData.taskEventEndTimes = runAnalysisInputsTmp.taskData.taskEventEndTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskData.fixationInTimes = runAnalysisInputsTmp.taskData.fixationInTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskData.fixationOutTimes = runAnalysisInputsTmp.taskData.fixationOutTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskData.juiceOnTimes = runAnalysisInputsTmp.taskData.juiceOnTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskData.juiceOffTimes = runAnalysisInputsTmp.taskData.juiceOffTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskData.gridPointsDegX = runAnalysisInputsTmp.taskData.gridPointsDegX;
    runAnalysisInputs.taskData.gridPointsDegY = runAnalysisInputsTmp.taskData.gridPointsDegY;
    runAnalysisInputs.taskData.taskEventIDs = runAnalysisInputsTmp.taskData.taskEventIDs;
    runAnalysisInputs.taskData.stimJumps = runAnalysisInputsTmp.taskData.stimJumps;
    runAnalysisInputs.taskData.stimFramesLost = runAnalysisInputsTmp.taskData.stimFramesLost;
    runAnalysisInputs.taskData.eventTimeAdjustments = runAnalysisInputsTmp.taskData.eventTimeAdjustments;
    
    runAnalysisInputs.taskDataAll.taskEventStartTimes = runAnalysisInputsTmp.taskDataAll.taskEventStartTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskDataAll.taskEventEndTimes = runAnalysisInputsTmp.taskDataAll.taskEventEndTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskDataAll.fixationInTimes = runAnalysisInputsTmp.taskDataAll.fixationInTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskDataAll.fixationOutTimes = runAnalysisInputsTmp.taskDataAll.fixationOutTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskDataAll.juiceOnTimes = runAnalysisInputsTmp.taskDataAll.juiceOnTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskDataAll.juiceOffTimes = runAnalysisInputsTmp.taskDataAll.juiceOffTimes + accumulatedTimeOffset;
    runAnalysisInputs.taskDataAll.gridPointsDegX = runAnalysisInputsTmp.taskDataAll.gridPointsDegX;
    runAnalysisInputs.taskDataAll.gridPointsDegY = runAnalysisInputsTmp.taskDataAll.gridPointsDegY;
    runAnalysisInputs.taskDataAll.taskEventIDs = runAnalysisInputsTmp.taskDataAll.taskEventIDs;
    runAnalysisInputs.taskDataAll.stimJumps = runAnalysisInputsTmp.taskDataAll.stimJumps;
    runAnalysisInputs.taskDataAll.stimFramesLost = runAnalysisInputsTmp.taskDataAll.stimFramesLost;
    runAnalysisInputs.taskDataAll.eventTimeAdjustments = runAnalysisInputsTmp.taskDataAll.eventTimeAdjustments;
     
    
    runAnalysisInputs.categoryList = union(runAnalysisInputs.categoryList,runAnalysisInputsTmp.categoryList); 
    runAnalysisInputs.eventLabels = union(runAnalysisInputs.eventLabels,runAnalysisInputsTmp.eventLabels); 
    
    for event_i = 1:length(runAnalysisInputs.eventLabels)
      runAnalysisInputs.jumpsByImage{event_i} = cat(1, runAnalysisInputs.jumpsByImage{event_i}, runAnalysisInputsTmp.jumpsByImage{event_i});
      %todo
      runAnalysisInputs.spikesByEvent = spikesByEvent; 
      runAnalysisInputs.psthEmptyByEvent = psthEmptyByEvent;
      runAnalysisInputs.spikesByEventForTF = spikesByEventForTF;
      runAnalysisInputs.lfpByEvent = lfpByEvent;
      runAnalysisInputs.analogInByEvent = analogInByEvent;
      runAnalysisInputs.onsetsByEvent = onsetsByEvent; 
      runAnalysisInputs.trialIDsByEvent = trialIDsByEvent;
      
    end
    
    for cat_i = 1:length(runAnalysisInputs.categoryList)
      %todo
      runAnalysisInputs.spikesByCategory = spikesByCategory; 
      runAnalysisInputs.psthEmptyByCategory = psthEmptyByCategory;
      runAnalysisInputs.spikesByCategoryForTF = spikesByCategoryForTF;
      runAnalysisInputs.lfpByCategory = lfpByCategory; 
      runAnalysisInputs.analogInByCategory = analogInByCategory;
      runAnalysisInputs.onsetsByCategory = onsetsByCategory;
      runAnalysisInputs.trialIDsByCategory = trialIDsByCategory;
      
    end
      

    runAnalysisInputs.stimTiming.shortest = min(runAnalysisInputs.stimTiming.shortest,runAnalysisInputsTmp.stimTiming.shortest);
    runAnalysisInputs.stimTiming.longest = min(runAnalysisInputs.stimTiming.longest,runAnalysisInputsTmp.stimTiming.longest);
    
    %this one is tricky; need to make sure it reflects any indexing updates
    %if event lists are different, e.g. if different stimulus subsets
    %presented
    runAnalysisInputs.eventCategories = eventCategories;  
    
    
    
    accumulatedTimeOffset = accumulatedTimeOffset + runAnalysisInputsTmp.excludeStimParams.ephysDuration;
  end
end


end

