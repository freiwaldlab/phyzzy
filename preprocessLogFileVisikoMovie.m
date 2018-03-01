function [ taskData, stimTiming ] = preprocessLogFileVisikoMovie(logFilename, taskTriggers, diodeTriggers, params )
%Special preprocessor for Visiko movie log files. 
%Sync log file's timestamps to blackrock clock, and store: 
%  - in presentation order, stimulus filenames, start/end times and jump RF mapping positions
%  - task event times: fixation in/out, fixspot flash start/end, juice delivery start/end 
%  - currently, uses digital IO packets for synchronization. In future
%    implementation, will optionally refine this sync with photodiode traces
%   Inputs:
%   - logFilename: an xml file generated by Visiko (or, in other
%     implementations, some other stimulation software)
%   - taskTriggers: digital IO data from blackrock NEV file
%   - params: currently must be struct of form params.usePhotodiode = 0
%   Outputs:
%   - taskData: struct with fields:
%       - taskEventIDs: nTrials x 1 cell array of alignment point identifiers, e.g. stimulus filenames. 
%         Used to retrieve information from stimulusParamFile.
%       - taskEventStartTimes: nTrials x 1 array of task event (e.g. stimulus) start times in ms
%       - taskEventEndTimes: nTrials x 1 array of task event (e.g. stimulus) end times in ms
%       _ stimParams: struct containing information about stimulus size etc. Not currently used elsewhere, so may be left empty
%       - RFmap: 1 for runs where stimulus position varies, else 0
%       - stimJumps: if RFmap, nTrialsx2 array or stim x and y positions, in degrees of visual angle, else may be empty
%           - note: position relative to center of jump grid, not absolute position
%       - gridPointsDegX: if RFmap, 1xN array of unique stimulus jump X locations, else may be empty 
%       - gridPointsDegY: if RFmap, 1xN array of unique stimulus jump Y locations, else may be empty
%       - fields likely required for excludeStimuli and runSummary plots (may be left empty if not needed):
%         - stimFramesLost: nTrials x 1, number of frames lost during trial
%         - fixationInTimes: nTrials x 1, times in ms when fixation epochs begin
%         - fixationOutTimes: nTrials x 1, times in ms when fixation epochs end
%         - juiceOnTimes: nTrials x 1, times in ms when juice delivery begins
%         - juiceOffTimes: nTrials x 1,times in ms when juice delivery ends
%         - fixSpotFlashStartTimes: nTrials x 1, times in ms when flash begins
%         - fixSpotFlashEndTimes: nTrials x 1, times in ms when flash ends
%   Dependencies:
%   - Statistics and Machine Learning Toolbox (for synchronizaton)
%   - xml2struct (from Matlab fileExchange)

if params.usePhotodiode
  error('photodiode synchronization not enabled');
end

disp('parsing serial IO packets');
packetTimes = taskTriggers.TimeStampSec;
packetData = dec2bin(taskTriggers.UnparsedData);
% note: bit (end-5) is 1 when fixating; changes in the last three bits
% (possibly more, haven't tested) signify the start of a new video

% note: the following assumes that the first trigger is not the start of a
% stimulation object: todo: make conditional, or include in params, or
% check log file
objectChangeIndexing = (abs(diff(packetData(:,end-2))) |  abs(diff(packetData(:,end-1))) | abs(diff(packetData(:,end))));
taskEventStartTimesBlk = 1000*packetTimes(vertcat(false,objectChangeIndexing) ~= 0);  %Blk affix signifies Blackrock reference frame
disp('number of stim and task-end triggers received by blackrock');
disp(length(taskEventStartTimesBlk));
clear packetData; 

% parse log file
disp('Loading visiko log file and converting to matlab structure');
assert(logical(exist(logFilename,'file')),'The stimulus log file you requested does not exist.');
if isfield(params,'tryPreparsedLogFile')
  try
    load(strcat(logFilename,'.mat'));
  catch
    logStruct = xml2struct(logFilename);
  end
else
  logStruct = xml2struct(logFilename);
end
assert(isfield(logStruct.VISIKOLOG,'EndStimulation'),'Error: stimulation end not included in log file');
assert(strcmp(logStruct.VISIKOLOG.Attributes.tasktype,'Bitmap Continuous'),'Error: unknown Visiko task type. Must be Bitmap Continuous');
if isa(logStruct.VISIKOLOG.DOCDATA,'cell')
  s = input(sprintf('Visiko parameters changed during task; %d parameter sets found. Enter the number to analyze, or n to quit: ',length(logStruct.VISIKOLOG.DOCDATA)),'s');
  if strcmp(s,'n')
    return;
  end
  logStruct.VISIKOLOG.DOCDATA = logStruct.VISIKOLOG.DOCDATA{str2double(s)};
end

%stimParams = logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT;
taskEventIDs = {};
% note: 'Objects' in the log file are either 'Picture' fields, which we
% want, or 'PictureCompleted' fields, which we don't care about. We
% initialize arrays to length(Objects), then, since many are non-negative, 
% we initialize to -1 then use >= 0 logical indexing to remove interlopers  
stimFramesLost = -1*ones(length(logStruct.VISIKOLOG.Object),1);
stimJumps = zeros(length(logStruct.VISIKOLOG.Object),2); %note: jumps can be negative, so we will use startTime for logical indexing
taskEventStartTimesLog = -1*ones(length(logStruct.VISIKOLOG.Object),1); % Log affix signifies stimulation computer reference frame
taskEventEndTimesLog = -1*ones(length(logStruct.VISIKOLOG.Object),1);
fixationInTimesLog = -1*ones(length(logStruct.VISIKOLOG.Trigger),1);
fixationOutTimesLog = -1*ones(length(logStruct.VISIKOLOG.Trigger),1);
juiceOnTimesLog = -1*ones(length(logStruct.VISIKOLOG.Trigger),1);
juiceOffTimesLog = -1*ones(length(logStruct.VISIKOLOG.Trigger),1);
if isfield(logStruct.VISIKOLOG,'FixspotFlash') 
  fixSpotFlashStartTimesLog = zeros(length(logStruct.VISIKOLOG.FixspotFlash),1);
  fixSpotFlashEndTimesLog = zeros(length(logStruct.VISIKOLOG.FixspotFlash),1);
else
  fixSpotFlashStartTimesLog = 0;
  fixSpotFlashEndTimesLog = 0;
end

for i = 1:length(logStruct.VISIKOLOG.Object)-1
  stimulusStruct = logStruct.VISIKOLOG.Object{i};
  tmp = regexp(stimulusStruct.Video.Filename.Text,'\','split');
  taskEventIDs = vertcat(taskEventIDs, strtrim(tmp{end})); %note: for some reason, filenames have trailing whitespace; trim it off
  stimFramesLost(i) = str2double(stimulusStruct.Frameslost.Text);
  stimJumps(i,:) = [0 0]; %stimJumps not implemented for videos; leave here for now in case expected somewhere else in code
  taskEventStartTimesLog(i) = str2double(stimulusStruct.Start.Text);
  taskEventEndTimesLog(i) = str2double(stimulusStruct.End.Text);
  % note: video struct also has fields FrameRate and FramePresentationTiming, if we want to be fancy
end
Output.VERBOSE(sprintf('number of stimulus trials: %s',length(taskEventIDs)));
% the following two commands account for two effects: first, we ignored the
% final, incomplete video above. Second, blackrock sometimes, but not
% always, receives a stimulation-end trigger from visiko, in which case we
% need to trim that trigger as well as the final video-start triger
assert(length(taskEventIDs) == length(taskEventStartTimesBlk)-1 || ...
  length(taskEventIDs) == length(taskEventStartTimesBlk)-2 ,'Inconsistent number of object starts in blackrock vs. Visiko');
taskEventStartTimesBlk = taskEventStartTimesBlk(1:length(taskEventIDs));
stimJumps = stimJumps(taskEventStartTimesLog >= 0,:); %note: jumps can be negative, so use startTime for logical indexing
stimFramesLost = stimFramesLost(stimFramesLost >= 0);
taskEventStartTimesLog = taskEventStartTimesLog(taskEventStartTimesLog >= 0);
taskEventEndTimesLog = taskEventEndTimesLog(taskEventEndTimesLog >= 0);

for i = 1:length(logStruct.VISIKOLOG.Trigger)
  triggerStruct = logStruct.VISIKOLOG.Trigger{i};
  switch triggerStruct.Name.Text
    case 'Fixation_in'
      fixationInTimesLog(i) = str2double(triggerStruct.Time.Text);
    case 'Fixation_out'
      fixationOutTimesLog(i) = str2double(triggerStruct.Time.Text);
    case 'RewardOn'
      juiceOnTimesLog(i)= str2double(triggerStruct.Time.Text);
    case 'RewardOff'
      juiceOffTimesLog(i)= str2double(triggerStruct.Time.Text);
    otherwise
      disp(strcat('unknown trigger name, ',TriggerStruct.Name.Text)); 
  end
end
fixationInTimesLog = fixationInTimesLog(fixationInTimesLog >= 0);
fixationOutTimesLog = fixationOutTimesLog(fixationOutTimesLog >= 0);
juiceOnTimesLog = juiceOnTimesLog(juiceOnTimesLog >= 0);
juiceOffTimesLog = juiceOffTimesLog(juiceOffTimesLog >= 0);

if isfield(logStruct.VISIKOLOG,'FixspotFlash')
  if iscell(logStruct.VISIKOLOG.FixspotFlash)
    for i = 1:length(logStruct.VISIKOLOG.FixspotFlash)
      flashStruct = logStruct.VISIKOLOG.FixspotFlash{i};
      fixSpotFlashStartTimesLog(i) = str2double(flashStruct.Time.Text);
      % note: flash duration is measured in sec in the log file; convert to msec
      fixSpotFlashEndTimesLog(i) = fixSpotFlashStartTimesLog(i) + 1000*str2double(flashStruct.Duration.Text);
    end
  else
    fixSpotFlashStartTimesLog(1) = str2double(logStruct.VISIKOLOG.FixspotFlash.Time.Text);
    fixSpotFlashEndTimesLog(1) = fixSpotFlashStartTimesLog(1) + 1000*str2double(logStruct.VISIKOLOG.FixspotFlash.Duration.Text);
  end
end

% now, calculate stimulation log to ephys clock conversion
if ~isfield(params,'syncMethod') || any(params.syncMethod == {'digitalTrigger','digitalTriggerNearestFrame'})
  % first, if nev has one more start trigger than log, throw out final nev
  % trigger (this is a known visiko bug, according to Michael Borisov's code)
  assert(length(taskEventStartTimesBlk) - length(taskEventStartTimesLog) <= 1, 'Error: Start triggers missing from log file'); %note: redundant, since we already checked they're equal
  taskEventStartTimesBlk = taskEventStartTimesBlk(1:length(taskEventStartTimesLog));
  %note: don't use first trigger in fit; sometimes off (known visiko bug)
  logVsBlkModel = fitlm(taskEventStartTimesBlk(2:end), taskEventStartTimesLog(2:end));
  disp(logVsBlkModel);
  m = logVsBlkModel.Coefficients.Estimate(2);
  y0 = logVsBlkModel.Coefficients.Estimate(1);
  % for debugging
  taskEventStartTimesFit = (1/m)*(taskEventStartTimesLog - y0);
  disp(strcat('Max magnitude fit residual, msec: ',num2str(max(abs(taskEventStartTimesBlk-taskEventStartTimesFit)))));
  % end for debugging
  taskEventEndTimesBlk = (1/m)*(taskEventEndTimesLog - y0);
  fixationInTimesBlk = (1/m)*(fixationInTimesLog - y0);
  fixationOutTimesBlk = (1/m)*(fixationOutTimesLog - y0);
  juiceOnTimesBlk = (1/m)*(juiceOnTimesLog - y0);
  juiceOffTimesBlk = (1/m)*(juiceOffTimesLog - y0);
  fixSpotFlashStartTimesBlk = (1/m)*(fixSpotFlashStartTimesLog - y0);
  fixSpotFlashEndTimesBlk = (1/m)*(fixSpotFlashEndTimesLog - y0);
end
if isfield(params,'syncMethod') && params.syncMethod == 'digitalTriggerNearestFrame'
  assert(~isempty(diodeTriggers.frameTimes),'You requested nearest frame correction, but did not supply a preprocesed photodiode trace.');
  for stim_i = 1:length(taskEventStartTimesBlk)
    [~,i] = min(abs(diodeTriggers.frameTimes-taskEventStartTimes(stim_i))); %or should we force it to be > 0?
    taskEventStartTimesBlk(stim_i) = frameTimes(i);
  end
end
if isfield(params,'syncMethod') && params.syncMethod == 'highFrames'
  error('highFrames sync method not yet implemented.');
end

% now, find within-video triggers
if isfield(params,'subTriggerArrayFilenames')
  subTriggerArrayFilenames = params.subTriggerArrayFilenames;
  for subTrigArray_i = 1:length(params.subTriggerArrayFilenames)
    tmp = load(subTriggerArrayFilenames{subTrigArray_i},'subTriggerArray');
    subTriggerArray = tmp.subTriggerArray;
    if params.keepTriggersAndSubTriggers
      warning('Trigger list will not be sorted by time; that will be implemetent in the future');
      newTaskEventIDs = taskEventIDs;
      newTaskEventStartTimesBlk = taskEventStartTimesBlk;
      newTaskEventEndTimesBlk = taskEventEndTimesBlk;
      newStimFramesLost = stimFramesLost;
    else
      newTaskEventIDs = {};
      newTaskEventStartTimesBlk = [];
      newTaskEventEndTimesBlk = [];
      newStimFramesLost = [];
    end
    for presentedVideo_i = 1:length(taskEventIDs)
      for trigArrayVideo_i = 1:length(subTriggerArray)
        if strcmp(taskEventIDs{presentedVideo_i},subTriggerArray{trigArrayVideo_i}{1})
          for subTrigger_i = 2:length(subTriggerArray{trigArrayVideo_i})
            newTaskEventIDs = vertcat(newTaskEventIDs,{subTriggerArray{trigArrayVideo_i}{subTrigger_i}{1}});
            newTaskEventStartTimesBlk = vertcat(newTaskEventStartTimesBlk,taskEventStartTimesBlk(presentedVideo_i) + subTriggerArray{trigArrayVideo_i}{subTrigger_i}{2});
            newStimFramesLost = vertcat(newStimFramesLost,stimFramesLost(presentedVideo_i));
            if length(subTriggerArray{trigArrayVideo_i}{subTrigger_i}) > 2
              newTaskEventEndTimesBlk = vertcat(newTaskEventEndTimesBlk,newTaskEventStartTimesBlk(end)+subTriggerArray{trigArrayVideo_i}{subTrigger_i}{3});
            else
              newTaskEventEndTimesBlk = vertcat(newTaskEventEndTimesBlk,newTaskEventStartTimesBlk(end)+1); %duration = 1 ms if no end time supplied
            end
          end
        end
      end
    end
    taskEventIDs = newTaskEventIDs;
    taskEventStartTimesBlk = newTaskEventStartTimesBlk;
    taskEventEndTimesBlk = newTaskEventEndTimesBlk;
    stimFramesLost = newStimFramesLost;
  end
end

%todo: re-sort task events for case where keepTriggersAndSubtriggers == 1
stimTiming.shortest = 1000*min(diff(sort(taskEventEndTimesBlk - taskEventStartTimesBlk, 'ascend'))); 
stimTiming.longest = 1000*max(diff(sort(taskEventEndTimesBlk - taskEventStartTimesBlk, 'ascend'))); 
if length(taskEventStartTimesBlk) > 1
  stimTiming.ISI = 1000*min(diff(sort(taskEventStartTimesBlk(2:end) - taskEventEndTimesBlk(1:end-1), 'ascend')));
else
  stimTiming.ISI = 0;
end

% finally, build the output structure
taskData.taskEventIDs = taskEventIDs;
taskData.stimJumps = stimJumps;
taskData.stimFramesLost = stimFramesLost;
taskData.taskEventStartTimes = taskEventStartTimesBlk;
taskData.taskEventEndTimes = taskEventEndTimesBlk;
taskData.fixationInTimes = fixationInTimesBlk;
taskData.fixationOutTimes = fixationOutTimesBlk;
taskData.juiceOnTimes = juiceOnTimesBlk;
taskData.juiceOffTimes = juiceOffTimesBlk;
taskData.fixSpotFlashStartTimes = fixSpotFlashStartTimesBlk;
taskData.fixSpotFlashEndTimes = fixSpotFlashEndTimesBlk;
taskData.stimParams = 0;
taskData.RFmap = 0;

%disp(taskEventIDs);

end
%
