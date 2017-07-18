function [ taskDataValid ] = excludeTrials( taskData, params )
%excludeStimuli removes stim based on
%   - fix spot flashes within flashPre before or flashPost after (msec)
%   - broken fixation within fixPre before or fixPost after (msec)
%   - juice delivery within juicePre before or juicePost after (msec)
%   - acceleration above thresh. within accelPre before or accelPost after
%
%   params has fields: 
%   - fixPre, fixPost, flashPre, flashPost (in ms)
%   and optional fields:
%   - juicePre, juicePost (in ms)
%   - accel1, accel2 (structs with fields data (1d timeseries, 1 ks/sec) and threshold);
%   - minStimDuration (ms)
%   todo: exclude stimuli shorter than minStimDuration (for arrythmic runs)

fixPre = params.fixPre;
fixPost = params.fixPost; 
flashPre = params.flashPre;  
flashPost = params.flashPost;
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

pictureValid = zeros(length(taskData.stimFilenames),1);
for i = 1:length(taskData.stimFilenames)
  if taskData.stimFramesLost(i) > 0
    continue
  end
  lastFixIn = max(taskData.fixationInTimes(taskData.fixationInTimes < taskData.stimStartTimes(i) - fixPre));
  lastFixOut = max(taskData.fixationOutTimes(taskData.fixationOutTimes < taskData.stimEndTimes(i) + fixPost));
  if isempty(lastFixIn)
    continue
  end
  if ~isempty(lastFixOut) && lastFixOut > lastFixIn %note: && short circuits, so won't throw error
    continue
  end
  lastFlashEnd = max(taskData.fixSpotFlashEndTimes(taskData.fixSpotFlashEndTimes < taskData.stimStartTimes(i) - flashPre));
  lastFlashStart = max(taskData.fixSpotFlashStartTimes(taskData.fixSpotFlashStartTimes < taskData.stimEndTimes(i) + flashPost));
  if ~isempty(lastFlashStart) && isempty(lastFlashEnd)
    continue
  end
  if ~isempty(lastFlashStart) && lastFlashStart > lastFlashEnd %note: short circuit
    continue
  end
  if exist('juicePre','var') && exist('juicePost','var')
    lastJuiceOff = max(taskData.juiceOffTimes(taskData.juiceOffTimes < taskData.stimStartTimes(i) - juicePre));
    lastJuiceOn = max(taskData.juiceOnTimes(taskData.juiceOnTimes < taskData.stimEndTimes(i) + juicePost));
    if ~isempty(lastJuiceOn) && isempty(lastJuiceOff)
      continue
    end
    if ~isempty(lastJuiceOn) && lastJuiceOn > lastJuiceOff %note: short circuit
      continue
    end
  end
  if exist('accel1','var')
    shakeOn = find(diff(accel1.data > accel1.threshold) > 0);
    if accel1.data(1) > accel1.threshold
      shakeOn = vertcat(1,shakeOn);
    end
    shakeOff = find(diff(accel1.data > accel1.threshold) < 0);
    lastShakeOff = max(shakeOff(shakeOff < taskData.stimStartTimes(i) - accel1.pre));
    lastShakeOn = max(shakeOn(shakeOn < taskData.stimEndTimes(i) + accel1.post));
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
    lastShakeOff = max(shakeOff(shakeOff < taskData.stimStartTimes(i) - accel2.pre));
    lastShakeOn = max(shakeOn(shakeOn < taskData.stimEndTimes(i) + accel2.post));
    if ~isempty(lastShakeOn) && isempty(lastShakeOff)
      continue
    end
    if ~isempty(lastShakeOn) && lastShakeOn > lastShakeOff %note: short circuit
      continue
    end
  end
  pictureValid(i) = 1;
end
%
disp(strcat('Percent of stimuli excluded: ', num2str(100-round(100*sum(pictureValid)/length(pictureValid)))));
disp(strcat('total remaining stimuli: ', num2str(sum(pictureValid))));
pictureValid = logical(pictureValid);
taskDataValid = taskData;
taskDataValid.stimFilenames = taskData.stimFilenames(pictureValid);
taskDataValid.stimJumps = taskData.stimJumps(pictureValid,:);
taskDataValid.stimFramesLost = taskData.stimFramesLost(pictureValid);
taskDataValid.stimStartTimes = taskData.stimStartTimes(pictureValid);
taskDataValid.stimEndTimes = taskData.stimEndTimes(pictureValid);

if params.DEBUG
  figure();
  hold on
  plot(taskDataValid.stimStartTimes,ones(size(taskDataValid.stimStartTimes)),'color','red','marker','o', 'linestyle','none');
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



