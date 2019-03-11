function [ taskDataValid ] = excludeTrialsMonkeyLogic( taskData, params )
%excludeStimuli removes stim based on
%   - fix spot flashes within flashPre before or flashPost after (msec)
%   - broken fixation within fixPre before or fixPost after (msec)
%   - juice delivery within juicePre before or juicePost after (msec)
%   - acceleration above thresh. within accelPre before or accelPost after
%   - photodiode vs. digital trigger alignment worse than threshold (worse
%     defined as difference from median value)
%
%   params has fields: 
%   - fixPre, fixPost, flashPre, flashPost (in ms)
%   and optional fields:
%   - juicePre, juicePost (in ms)
%   - accel1, accel2 (structs with fields data (1d timeseries, 1 ks/sec) and threshold);
%   - maxEventTimeAdjustmentDeviation (in ms)
%   - minStimDuration (ms)
%   todo: exclude stimuli shorter than minStimDuration (for arrythmic runs)

if isfield(params,'needExcludeTrials') && ~params.needExcludeTrials
  taskDataValid = taskData;
  return
end

fixPre = params.fixPre;
fixPost = params.fixPost; 
flashPre = params.flashPre;  
flashPost = params.flashPost;
if isfield(params,'ephysDuration')
  ephysDuration = params.ephysDuration;
else
  ephysDuration = max(taskData.taskEventStartTimes);
  disp('No ephysDuration supplied to excludeParams; using final task event start time as duration lower bound.');
end
if isfield(params, 'juicePre') && isfield(params, 'juicePost')
  juicePre = params.juicePre;
  juicePost = params.juicePost;
end
if isfield(params, 'accel1')
  accel1 = params.accel1;
end
if isfield(params, 'accel2')
  accel2 = params.accel2;
end
if isfield(params,'minStimDur')
  minStimDur = params.minStimDur;
end
if isfield(params, 'maxEventTimeAdjustmentDeviation') && isfield(taskData,'eventTimeAdjustments')
  maxEventTimeAdjustmentDeviation = params.maxEventTimeAdjustmentDeviation;
  eventTimeAdjustmentDeviations = abs(taskData.eventTimeAdjustments - median(taskData.eventTimeAdjustments));
end
if isfield(params, 'ephysDataPre')
 ephysDataPre = params.ephysDataPre;
else
  ephysDataPre = fixPre;
end
if isfield(params, 'ephysDataPost')
 ephysDataPost = params.ephysDataPost;
else
  ephysDataPost = fixPost;
end

trialValid = zeros(length(taskData.taskEventIDs),1);
%initialize array of 0s, add numbers to the same Ntrials*1 array. Non-zeros
%will be excluded at the end.

%exclude failed trials
if isfield(params, 'excludeFailed') && params.excludeFailed
  trialValid = taskData.errorArray + trialValid;
end

%exclude trials where too many frames are dropped;
if isfield(params, 'frameDropThreshold')
  trialValid = trialValid + floor(taskData.stimFramesLost/params.frameDropThreshold);
end

trialValid = (trialValid == 0);

%
fprintf('Percent of trials excluded: %f\n', sum(trialValid)/length(trialValid)*100)
fprintf('Percent of trials excluded: %d\n', sum(trialValid))
trialValid = logical(trialValid);
taskDataValid = struct; %taskDataValid = taskData;
taskDataValid.translationTable = taskData.translationTable;
taskDataValid.frameMotionData = taskData.frameMotionData;
taskDataValid.taskEventIDs = taskData.taskEventIDs(trialValid);
taskDataValid.stimFramesLost = taskData.stimFramesLost(trialValid);
taskDataValid.taskEventStartTimes = taskData.taskEventStartTimes(trialValid);
taskDataValid.taskEventEndTimes = taskData.taskEventEndTimes(trialValid);
taskDataValid.taskEventStartTimes = taskData.taskEventStartTimes(trialValid);
taskDataValid.juiceOnTimes = taskData.juiceOnTimes(trialValid);
taskDataValid.juiceOffTimes = taskData.juiceOffTimes(trialValid);
taskDataValid.RFmap = taskData.RFmap;
taskDataValid.RFmap = taskData.RFmap;


if params.DEBUG
  figure();
  hold on
  plot(taskDataValid.taskEventStartTimes,ones(size(taskDataValid.taskEventStartTimes)),'color','red','marker','o', 'linestyle','none');
  plot(taskData.juiceOnTimes, 3*ones(size(taskData.juiceOnTimes)),'color','red','marker','o', 'linestyle','none');
  plot(taskData.juiceOffTimes, 3*ones(size(taskData.juiceOffTimes)),'color','green','marker','o', 'linestyle','none');
  plot([0,0],[-20 20],'marker','none');
  hold off
end
end



