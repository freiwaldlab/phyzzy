function [ taskDataValid ] = excludeTrials( taskData, params )
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

if isfield(params,'excludeProcessor')
  [taskDataValid] = params.excludeProcessor( taskData, params );
  return
end

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
for trial_i = 1:length(taskData.taskEventIDs)
  if taskData.stimFramesLost(trial_i) > 0
    continue
  end
  lastFixIn = max(taskData.fixationInTimes(taskData.fixationInTimes < taskData.taskEventStartTimes(trial_i) - fixPre));
  lastFixOut = max(taskData.fixationOutTimes(taskData.fixationOutTimes < taskData.taskEventEndTimes(trial_i) + fixPost));
  if isempty(lastFixIn)
    continue
  end
  if ~isempty(lastFixOut) && lastFixOut > lastFixIn %note: && short circuits, so won't throw error
    continue
  end
  lastFlashEnd = max(taskData.fixSpotFlashEndTimes(taskData.fixSpotFlashEndTimes < taskData.taskEventStartTimes(trial_i) - flashPre));
  lastFlashStart = max(taskData.fixSpotFlashStartTimes(taskData.fixSpotFlashStartTimes < taskData.taskEventEndTimes(trial_i) + flashPost));
  if ~isempty(lastFlashStart) && isempty(lastFlashEnd)
    continue
  end
  if ~isempty(lastFlashStart) && lastFlashStart > lastFlashEnd %note: short circuit
    continue
  end
  if exist('juicePre','var') && exist('juicePost','var')
    lastJuiceOff = max(taskData.juiceOffTimes(taskData.juiceOffTimes < taskData.taskEventStartTimes(trial_i) - juicePre));
    lastJuiceOn = max(taskData.juiceOnTimes(taskData.juiceOnTimes < taskData.taskEventEndTimes(trial_i) + juicePost));
    if ~isempty(lastJuiceOn) && isempty(lastJuiceOff)
      continue
    end
    if ~isempty(lastJuiceOn) && lastJuiceOn > lastJuiceOff %note: short circuit
      continue
    end
  end
  if exist('maxEventTimeAdjustmentDeviation','var') && exist('eventTimeAdjustmentDeviations','var')
    if eventTimeAdjustmentDeviations(trial_i) > maxEventTimeAdjustmentDeviation
      continue
    end
  end
  if exist('accel1','var')
    shakeOn = find(diff(accel1.data > accel1.threshold) > 0);
    if accel1.data(1) > accel1.threshold
      shakeOn = vertcat(1,shakeOn);
    end
    shakeOff = find(diff(accel1.data > accel1.threshold) < 0);
    lastShakeOff = max(shakeOff(shakeOff < taskData.taskEventStartTimes(trial_i) - accel1.pre));
    lastShakeOn = max(shakeOn(shakeOn < taskData.taskEventEndTimes(trial_i) + accel1.post));
    if ~isempty(lastShakeOn) && isempty(lastShakeOff)
      continue
    end
    if ~isempty(lastShakeOn) && lastShakeOn > lastShakeOff %note: short circuit
      continue
    end
  end
  if exist('accel2','var')
    shakeOn = find(diff(accel2.data > accel2.threshold) > 0);
    if accel2.data(1) > accel2.threshold
      shakeOn = vertcat(1,shakeOn);
    end
    shakeOff = find(diff(accel2.data > accel2.threshold) < 0);
    lastShakeOff = max(shakeOff(shakeOff < taskData.taskEventStartTimes(trial_i) - accel2.pre));
    lastShakeOn = max(shakeOn(shakeOn < taskData.taskEventEndTimes(trial_i) + accel2.post));
    if ~isempty(lastShakeOn) && isempty(lastShakeOff)
      continue
    end
    if ~isempty(lastShakeOn) && lastShakeOn > lastShakeOff %note: short circuit
      continue
    end
  end
  if (taskData.taskEventStartTimes(trial_i) - ephysDataPre) <= 0
    continue
  end
  if taskData.taskEventEndTimes(trial_i) + ephysDataPost > ephysDuration
    continue
  end
  trialValid(trial_i) = 1;
end
%
disp(strcat('Percent of trials excluded: ', num2str(100-round(100*sum(trialValid)/length(trialValid)))));
disp(strcat('total remaining trials: ', num2str(sum(trialValid))));
trialValid = logical(trialValid);
taskDataValid = taskData;
taskDataValid.taskEventIDs = taskData.taskEventIDs(trialValid);
taskDataValid.stimJumps = taskData.stimJumps(trialValid,:);
taskDataValid.stimFramesLost = taskData.stimFramesLost(trialValid);
taskDataValid.taskEventStartTimes = taskData.taskEventStartTimes(trialValid);
taskDataValid.taskEventEndTimes = taskData.taskEventEndTimes(trialValid);

if params.DEBUG
  figure();
  hold on
  plot(taskDataValid.taskEventStartTimes,ones(size(taskDataValid.taskEventStartTimes)),'color','red','marker','o', 'linestyle','none');
  plot(taskData.fixationInTimes, 2*ones(size(taskData.fixationInTimes)),'color','green','marker','o','linestyle','none');
  plot(taskData.fixationOutTimes,2*ones(size(taskData.fixationOutTimes)),'color','red','marker','o', 'linestyle','none');
  plot(taskData.juiceOnTimes, 3*ones(size(taskData.juiceOnTimes)),'color','red','marker','o', 'linestyle','none');
  plot(taskData.juiceOffTimes, 3*ones(size(taskData.juiceOffTimes)),'color','green','marker','o', 'linestyle','none');
  plot(taskData.fixSpotFlashStartTimes, 4*ones(size(taskData.fixSpotFlashStartTimes)),'color','red','marker','o', 'linestyle','none');
  plot(taskData.fixSpotFlashEndTimes, 4*ones(size(taskData.fixSpotFlashEndTimes)),'color','green','marker','o', 'linestyle','none');
  plot([0,0],[-20 20],'marker','none');
  hold off
end
end



