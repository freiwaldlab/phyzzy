function [ taskData, stimTiming ] = preprocessLogFileMonkeyLogic(logfile, taskTriggers, diodeTriggers, params)
%% Reads in Black
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

%% Process MonkeyLogic File
%Parse the Log file
disp('Loading MonkeyLogic log file');
assert(logical(exist(logfile,'file')),'The logfile you requested does not exist.');
[data,MLConfig,TrialRecord] = mlread(logfile);

%Find stimulus timing range - in ML, we have no range within valid trials.
stimTiming.shortest = 2800; %Making this 0 for now
stimTiming.longest = 2800; %setting this to 2800 for now, but maybe save as an editable variable.
stimTiming.ISI = MLConfig.InterTrialInterval;

%Construct the translation table from the logfile
%Find the Conditions that need to be connected to codes
allConditions = unique(TrialRecord.ConditionsPlayed); %Pull unique members of this list
translationTable = cell(length(allConditions),1); %Initialize translation table
trial_ind = 1;

%Check for empty members of the translation table. 
while sum(find(cellfun('isempty', translationTable))) ~= 0
    trialHolder = data(trial_ind);
    stimName = trialHolder.TaskObject.Attribute(2).Name; %Pull the string containing the stimulus name.
    translationTable{trialHolder.Condition} = stimName(~isspace(stimName));
    trial_ind = trial_ind + 1;
end

%Pull Behavioral codes for other events
behavioralCodes = TrialRecord.TaskInfo.BehavioralCodes;

fixCueMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Fix Cue'));
trialStartMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Start trial'));
stimStartMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Stimuli On'));
stimEndMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Stimuli Off'));
rewardMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Reward'));
frameSkipMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Frame skipped'));
trialEndMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'End trial'));
manualRewardMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Manual reward'));

%MonkeyLogic's output has an absolute trial start (time within entire run)
%and additional time stamps within trials, relative to that start. Collect
%all the start times, then create subsequent ones through the addition of
%the right amount. Step below removes failed trials.
trueTrialArray = (TrialRecord.TrialErrors == 0);
mklTrialStarts = [data(trueTrialArray).AbsoluteTrialStartTime];
tmpStruct = [data(trueTrialArray).BehavioralCodes];
rwdTimes = [data(trueTrialArray).RewardRecord];

[taskEventStartTimesLog, taskEventEndTimesLog, juiceOnTimesLog, juiceOffTimesLog] = deal(zeros(sum(trueTrialArray), 1));

taskEventIDsLog = TrialRecord.ConditionsPlayed(trueTrialArray)';
  
for ii = 1:length(mklTrialStarts)
  taskEventStartTimesLog(ii) = mklTrialStarts(ii) + tmpStruct(ii).CodeTimes(tmpStruct(ii).CodeNumbers == stimStartMarker);
  taskEventEndTimesLog(ii) = mklTrialStarts(ii) + tmpStruct(ii).CodeTimes(tmpStruct(ii).CodeNumbers == stimEndMarker);
  juiceOnTimesLog(ii) = mklTrialStarts(ii) + rwdTimes(ii).StartTimes;
  juiceOffTimesLog(ii) = mklTrialStarts(ii) + rwdTimes(ii).EndTimes;
end

% Process the Eye data into a structure which can be used later to line up
% with Blackrock signals.
tmpEye = [data(trueTrialArray).AnalogData];
tmpEye = rmfield(tmpEye, {'SampleInterval', 'EyeExtra','Joystick','Mouse','PhotoDiode','General','Button'});

% Behavioral summary of performance during recording
behaviorsummaryPhyzzy(logfile)
[~, filename, ~] = fileparts(logfile);
savefig(sprintf('%sBehavSum_%s',params.outDir,filename));

%% Process the NEV derived data (Blackrock's Eventmarkers)
disp('parsing serisal IO packets');
packetTimes = double(taskTriggers.TimeStampSec)*1000; %convert seconds to milliseconds
packetData = double(taskTriggers.UnparsedData);

if ~isempty(packetData) % Means Blackrock/MKL Communication was intact, correct
  %Code to get rid of any markers prior to the first trial beginning (means
  % blackrock turned on AFTER beginning MKL.
  lol = find(packetData == trialStartMarker);
  packetData = packetData(lol(1):end);
  packetTimes = packetTimes(lol(1):end);
  
  %Comb through the Blackrock data and pull codes/times associated w/ valid
  %trials.
  [taskEventStartTimesBlk, taskEventEndTimesBlk, juiceOnTimesBlk, juiceOffTimesBlk, taskEventIDsBlk] = deal(zeros(sum(packetData == 40), 1));
  trueTrialcount = 1;
  
  for ii = 1:length(packetData)
    if packetData(ii) > 100
      stimCondTemp = packetData(ii);
    elseif packetData(ii) == stimStartMarker
      stimStartTemp = packetTimes(ii);
    elseif packetData(ii) == stimEndMarker
      stimEndTemp = packetTimes(ii);
    elseif packetData(ii) == rewardMarker %This assumes the "juice end time" is right after this marker.
      taskEventIDsBlk(trueTrialcount) = stimCondTemp;
      taskEventStartTimesBlk(trueTrialcount) = stimStartTemp;
      taskEventEndTimesBlk(trueTrialcount) = stimEndTemp;
      juiceOnTimesBlk(trueTrialcount) = packetTimes(ii);
      juiceOffTimesBlk(trueTrialcount) = packetTimes(ii + 1);
      trueTrialcount = trueTrialcount + 1;
    end
  end
  
  %for every trial, find the conditions number and convert it to the
  %appropriate filename for the stimulus, based on the translation table created earlier.
  taskEventIDsBlk = taskEventIDsBlk - 100; %This is a number I set in my timing script.
  taskEventIDs = cell(1, length(taskEventIDsBlk));
  for ii = 1:length(taskEventIDsBlk)
    taskEventIDs(ii) = translationTable(taskEventIDsBlk(ii)); %Reference the translation table w/ the condition ID as an index
  end
  
  
else % Means Blackrock/MKL Communication was not correctly connected.
  error('The Blackrock Digital inputs are empty. Digital inputs may have not been plugged in.');
  %disp('The Blackrock Digital inputs are empty. Digital inputs may have not been plugged in. Using MKL Time stamps.');
  [taskEventStartTimesBlk, taskEventEndTimesBlk, juiceOnTimesBlk, juiceOffTimesBlk, taskEventIDsBlk] = deal(zeros(sum(TrialRecord.TrialErrors == 0), 1));
  
  %Since there are no behavioral timestamps, you need to line these codes
  %up with something. The strobe is the only non-digital signal which is directed and
  %simply "lined up" to MonkeyLogic's control of trials. Assumptions must
  %be made. as of Oct '18, the grey background is "mid", the strobe in
  %white is high, and the task has the strobe @ low, or black. Every
  %trial has the fixation period, most proceed past it, few don't. Proper
  %indexing of the strobe will rely on capturing the fixation period.
  
  %Grab the array of start times from the data structure.
  %[dataConcat, ~, ~] = mlconcatenate(logfile);
  
  %Find the time stamp of the first trial's start time.
  firstTrialFix = data(1).BehavioralCodes.CodeTimes(find(data(1).BehavioralCodes.CodeNumbers == fixCueMarker));
  
  % This value should line up with the first switching on of the strobe,
  % and if recording of the strobe began after the MKL Exp was started
  % (but not initiated, before Trial 1), screen is grey, strobe is white
  % for fix, so 'midtoHigh' first marker is what we want to match this
  % to.
  strobeToMKLOffset = diodeTriggers.midToHigh - firstTrialFix;
  
  %Add this offset to all the startTimes.
  %taskEventStartTimesBlk = vertcat(dataConcat.AbsoluteTrialStartTime);
  
  
  %taskEventStartTimesBlk = vertcat(data(1:end).AbsoluteTrialStartTime) + strobeToMKLOffset;
  
  %Shift all of those start times by the amount between the start of the
  %blackrock (First transition from mid/high). This may not work in
  %scenarios where the Blackrock was turned on well before the trial
  %began (and if anything came on the screen to mess with it). (Assuming
  %these timestamps are in the same units).
  
  
  %Idea 1 - use the transition from Grey to White (Before task to First
  %Fix window) as a source of Lag, and add it to every time stamp once it
  %is collected.
end

%We now have the start of trials, Blackrock.Eventmarkers. We want to move
%that to Blackrock.Strobe, which we believe is more accurate. This will be done by finding the transitions for the
%strobe from High (white) to low (black) closest to the time stamp.

%What is the discrepency?
offsets = zeros(sum(packetData == rewardMarker),1);

%Using Photodiode instead of linear model 
if params.usePhotodiode
  if length(diodeTriggers.low) < length(diodeTriggers.mid)
    % For each of the event start times, find the closest transition in
    % the photodiode and overwrite the time stamp with this time.
    for ii = 1:length(taskEventStartTimesBlk)
      [offsets(ii), ind_truestart] = min(abs(diodeTriggers.highToMid-taskEventStartTimesBlk(ii)));
      taskEventStartTimesBlk(ii) = diodeTriggers.highToMid(ind_truestart);
    end
    %do the same for the end times, except referencing the nearest Low to high
    %transition.
    for ii = 1:length(taskEventEndTimesBlk)
      [offsets(ii), ind_truestart] = min(abs(diodeTriggers.midToHigh-taskEventEndTimesBlk(ii)));
      taskEventEndTimesBlk(ii) = diodeTriggers.midToHigh(ind_truestart);
    end
  else
    for ii = 1:length(taskEventStartTimesBlk)
      [offsets(ii), ind_truestart] = min(abs(diodeTriggers.highToLow-taskEventStartTimesBlk(ii)));
      taskEventStartTimesBlk(ii) = diodeTriggers.highToLow(ind_truestart);
    end
    %do the same for the end times, except referencing the nearest Low to high
    %transition.
    for ii = 1:length(taskEventEndTimesBlk)
      [offsets(ii), ind_truestart] = min(abs(diodeTriggers.lowToHigh-taskEventEndTimesBlk(ii)));
      taskEventEndTimesBlk(ii) = diodeTriggers.lowToHigh(ind_truestart);
    end
  end
end

%% Now, calculate stimulation log to ephys clock conversion 
  %Prelim Model - All Trial start times
  trialStartTimesBlk = packetTimes(packetData == 9);
  trialStartTimesBlk = trialStartTimesBlk(trueTrialArray);
  trialStartTimesLog = [data(trueTrialArray).AbsoluteTrialStartTime]';
  tmpStruct = [data(trueTrialArray).BehavioralCodes]';
  for ii = 1:length(trialStartTimesLog)
    trialStartTimesLog(ii) = trialStartTimesLog(ii) + tmpStruct(ii).CodeTimes(1);
  end
  
  logVsBlkModelTrial = fitlm(trialStartTimesBlk, trialStartTimesLog);
  disp(logVsBlkModelTrial);

  %Make the Model by comparing Blackrock event start times to monkeylogic.
  logVsBlkModel = fitlm(taskEventStartTimesBlk, taskEventStartTimesLog);
  disp(logVsBlkModel)
  
  m = logVsBlkModel.Coefficients.Estimate(2);
  y0 = logVsBlkModel.Coefficients.Estimate(1);
  taskEventStartTimesFit = (1/m)*(taskEventStartTimesLog - y0);
  
  disp(strcat('Max magnitude fit residual, msec: ',num2str(max(abs(taskEventStartTimesBlk-taskEventStartTimesFit)))));
  %Shift over all the events, using the calculated line.
  eventTimeAdjustments = taskEventStartTimesFit-taskEventStartTimesBlk;
  
%   figure();
%   hist(eventTimeAdjustments,40);
%   title('Sync adjustments from photodiode');
%   xlabel('offset (adjusted - original) (ms)');
%   ylabel('count');
%   disp(median(eventTimeAdjustments));
%   [~,adjSortInds] = sort(abs(eventTimeAdjustments-median(eventTimeAdjustments)), 'descend');
%   disp('worst alignments, log file times');
%   disp(taskEventStartTimesLog(adjSortInds(1:min(5,length(adjSortInds)))));
%   disp('worst alignments, adjusted times');
%   disp(taskEventStartTimesBlk(adjSortInds(1:min(5,length(adjSortInds)))));
%   disp('worst alignments, adjustment values (ms)');
%   disp(eventTimeAdjustments(adjSortInds(1:min(5,length(adjSortInds)))));
  
  %Assuming you're happy with this model, shift values only in the log to
  %Blackrock time.
  juiceOnTimesBlk = (1/m)*(juiceOnTimesLog - y0);
  juiceOffTimesBlk = (1/m)*(juiceOffTimesLog - y0);
  %Not doing this for event start times, since I trust the photodiode more
  %than the eventmarker on either machine. 
  
  fprintf('average offset %s ms\n', num2str(mean(offsets)));
  fprintf('range of offset %d ms - %d ms \n', [min(offsets), max(offsets)])
%% Output
%Adding random numbers to these - they aren't relevant for my current task,
%nor are they directly recorded by MKL.
fixSpotFlashStartTimesBlk = taskEventStartTimesBlk(1);
fixSpotFlashEndTimesBlk = taskEventEndTimesBlk(1);
fixationInTimesBlk = taskEventStartTimesBlk(1);
fixationOutTimesBlk = taskEventEndTimesBlk(1);

% finally, build the output structure
taskData.taskEventIDs = taskEventIDs';
%taskData.stimJumps = stimJumps;
%taskData.stimFramesLost = stimFramesLost;
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
taskData.eyeData = tmpEye;
taskData.eyeCal.origin = MLConfig.EyeTransform{2}.origin;
taskData.eyeCal.gain = MLConfig.EyeTransform{2}.gain;
taskData.mklTrialStarts = mklTrialStarts;
taskData.NEVTrialTimes = trialStartTimesBlk;


end
%
