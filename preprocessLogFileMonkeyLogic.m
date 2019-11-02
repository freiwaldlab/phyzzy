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
assert(length(unique(TrialRecord.ConditionsPlayed)) < ceil((length(TrialRecord.ConditionsPlayed))/2) , 'MonkeyLogic file reports each condition was not repeated at least 2 times')

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

%Clean up the stim table in the logfile. 
for ii = 1:length(TrialRecord.TaskInfo.Stimuli)
  logFileStimTable{ii,1} = TrialRecord.TaskInfo.Stimuli{ii}(~isspace(TrialRecord.TaskInfo.Stimuli{ii}));
end

%Check to see if what you collected by cycling through the stim matches
%this table.
assert(sum(strcmp(logFileStimTable,translationTable)) == length(translationTable), 'Stim table from MonkeyLogic doesnt match collected table');

%Pull Behavioral codes for other events
behavioralCodes = TrialRecord.TaskInfo.BehavioralCodes;

fixCueMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Fix Cue'));
trialStartMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Start trial'));
stimStartMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Stimuli On'));
stimEndMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Stimuli Off')); %Also means fixation off.
rewardMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Reward'));
frameSkipMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Frame skipped'));
trialEndMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'End trial'));
manualRewardMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Manual reward'));
fixFailMarker = 4;
stimFailMarker = 3;

%MonkeyLogic's output has an absolute trial start (time within entire run)
%and additional time stamps within trials, relative to that start. Collect
%all the start times, then create subsequent ones through the addition of
%the right amount. Step below removes failed trials.
errorArray = TrialRecord.TrialErrors';
correctTrialArray = (errorArray == 0);
mklTrialStarts = [data.AbsoluteTrialStartTime];
tmpStruct = [data.BehavioralCodes];
rwdTimes = [data.RewardRecord];

%initialize everything
[taskEventStartTimesLog, taskEventEndTimesLog, juiceOnTimesLog, juiceOffTimesLog, stimFramesLost] = deal(zeros(sum(correctTrialArray), 1));

%Collect trialEventIDs.
taskEventIDsLog = TrialRecord.ConditionsPlayed';
  
%Note - the below setup pulls trials which fail during the stimuli, but not
%during the fixation period. for those trials, both start and end are NaN,
%for stim failures, the numbers are meaningful.
for ii = 1:length(mklTrialStarts)
  if ~isempty(tmpStruct(ii).CodeTimes(tmpStruct(ii).CodeNumbers == stimStartMarker))
    taskEventStartTimesLog(ii) = mklTrialStarts(ii) + tmpStruct(ii).CodeTimes(tmpStruct(ii).CodeNumbers == stimStartMarker);
    taskEventEndTimesLog(ii) = mklTrialStarts(ii) + tmpStruct(ii).CodeTimes(tmpStruct(ii).CodeNumbers == stimEndMarker);
  else
    taskEventStartTimesLog(ii) = nan;
    taskEventEndTimesLog(ii) = nan;
  end
  if ~isempty(rwdTimes(ii).StartTimes)
    juiceOnTimesLog(ii) = mklTrialStarts(ii) + rwdTimes(ii).StartTimes(end);
    juiceOffTimesLog(ii) = mklTrialStarts(ii) + rwdTimes(ii).EndTimes(end);
  else
    juiceOnTimesLog(ii) = nan;
    juiceOffTimesLog(ii) = nan;
  end
end

% Process the Eye data into a structure which can be used later to line up
% with Blackrock signals.
tmpEye = [data.AnalogData];
tmpEye = rmfield(tmpEye, {'SampleInterval', 'EyeExtra','Joystick','Mouse','PhotoDiode','General','Button'});

%Package the information necessary for proper eye/stimuli relationship.
tmp = split(MLConfig.Resolution,' ');
screenStats.screenSize = {str2double(tmp{1}), str2double(tmp{3})};
screenStats.PixelsPerDegree = abs(MLConfig.PixelsPerDegree);
screenStats.DiagonalSize = MLConfig.DiagonalSize;

%Translations - Function changes the name of stimuli from an Old convention
translationTable = updateTranslationTable(translationTable);

%Create Frame lost array
for ii = 1:length(data)
    stimFramesLost(ii) = sum(data(ii).BehavioralCodes.CodeNumbers == frameSkipMarker);
end

% Behavioral summary of performance during recording
taskData = struct();
taskData = behaviorsummaryPhyzzy(logfile, taskData);
[~, filename, ~] = fileparts(logfile);
savefig(sprintf('%sBehavSum_%s',params.outDir,filename));

%% Process the NEV derived data (Blackrock's Eventmarkers)
disp('parsing serisal IO packets');
packetTimes = double(taskTriggers.TimeStampSec)*1000; %convert seconds to milliseconds
packetData = double(taskTriggers.UnparsedData);

if ~isempty(packetData) % Means Blackrock/MKL Communication was intact, correct
  %Code to get rid of any markers prior to the first trial beginning, or after the final trial end.
  %a defense against double 9's
  trueStart = find(packetData == trialStartMarker);
  trueStartTrials = (diff(trueStart) > 1);
  trueStart = trueStart(trueStartTrials);
  trueStart = trueStart(1);
  trueEnd = find(packetData == trialEndMarker,1,'last');
  
  packetData = packetData(trueStart:trueEnd);
  packetTimes = packetTimes(trueStart:trueEnd);
  
%   if trueEnd ~= length(packetData) || trueStart ~= 1
%       %This happens if Blackrock is on less time than Monkeylogic -
%       %fragments of complete trial marker sets arrive at Blackrock. to keep
%       %things complete, we need to see which trials are cut, assuming
%       %Blackrock 
%       origStartInd = find(packetData == trialStartMarker);
%       origEndInd = find(packetData == trialEndMarker);
%       
%       packetData = packetData(trueStart:trueEnd);
%       packetTimes = packetTimes(trueStart:trueEnd);
%       
%       cutStartInd = find(packetData == trialStartMarker);
%       cutEndInd = find(packetData == trialEndMarker);
% 
%       indStartCutTrial = setdiff(origStartInd, cutStartInd);
%       indEndCutTrial = setdiff(origEndInd, cutEndInd);
%       indCutTrial = [indStartCutTrial indEndCutTrial];
%       %Use this index to modify relevant vectors from Monkeylogic.
%       keptTrialInd = ones(length(origStartInd),1);
%       for ii = 1:length(indCutTrial)
%         keptTrialInd(origStartInd == indCutTrial(ii)) = 0;
%       end
%   end
  
  %to later figure out how to shape data, we need to know what we saw and
  %got rid of.
  
  %Below are things which should happen every trial  
  trialStartInds = find(packetData == trialStartMarker); 
  trialEndInds = find(packetData == trialEndMarker); 
  condNumbers = packetData(packetData > 100)-100; 
  fixStartInds = find(packetData == fixCueMarker);
  stimFixEndInds = find(packetData == stimEndMarker);
  
  assert(length(trialStartInds) == length(condNumbers), 'Every trial doesnt have a condition number - this may be due to changes in when the marker is sent.')
  
  %Now, use the correct trial array from MKL to pick out only the correct
  %trials from the whole array. Assign them them to the correct vectors. 
  [taskEventStartTimesBlk, taskEventEndTimesBlk, juiceOnTimesBlk, juiceOffTimesBlk, taskEventIDsBlk] = deal(nan(length(trialStartInds), 1));
  
  %Construct Error array from Blackrock files only, useful in the case that
  %incomplete trials lead to unpaired start and stop triggers.
  errorArray = zeros(sum(packetData == trialStartMarker), 1);
  for ii = 1:length(trialStartInds) %for every trial
    %Go through packet data, starting at that stim, and see if you hit a
    %40, 3, or 4 first.
    found = 0;
    stepsAhead = 1;
    while ~found
      switch packetData(trialStartInds(ii)+stepsAhead)
        case fixFailMarker
          errorArray(ii) = fixFailMarker;
          found = 1;
        case stimFailMarker
          errorArray(ii) = stimFailMarker;
          found = 1;
        case rewardMarker
          errorArray(ii) = 0;
          found = 1;
        otherwise
          stepsAhead = stepsAhead + 1;
      end
    end
  end
  
  %packetData is directly referenced for events which only happen
  %sometimes.
  taskEventIDsBlk = condNumbers;
  taskEventStartTimesBlk(errorArray ~= 4) = packetTimes(packetData == stimStartMarker); %Stores every trial where stim started (no fail during fix);
  taskEventEndTimesBlk(errorArray ~= 4) = packetTimes(stimFixEndInds(errorArray ~= 4));
  juiceOnTimesBlk(errorArray == 0) = packetTimes(packetData == rewardMarker);
  juiceOffTimesBlk(errorArray == 0) = packetTimes(trialEndInds(errorArray == 0));
  
  %use the condition number to fill out the taskEventsIDs using the
  %translation table.
  taskEventIDs = translationTable(condNumbers);
  
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

%% Run some checks comparing data collected from Blackrock, MonkeyLogic.
if (length(taskEventIDsLog) ~= length(taskEventIDsBlk))
  warning('Blackrock and MonkeyLogic report different numbers of correct trials. Attempting to find best alignment')
  %Odd situation, likely due to hitting runAnalyses prior to the MKL being
  %stopped. Seems like MKL doesn't have the complete record, Blackrock
  %does. In these cases, use the shorter of the two.
  lag = finddelay(taskEventIDsBlk, taskEventIDsLog);
  %Switch below written for a single case in July, unlikely to happen
  %again, needs more extensive testing.
  %If find delay has an answer - this might invalidate the length checking method, but it seems to not always work.
  %This means shifting things forward improves the pairwise match - means
  %something extra at the start.
  if (lag > 0)
    taskEventIDsLog = taskEventIDsLog((lag+1):end);
    taskEventStartTimesLog = taskEventStartTimesLog((lag+1):end);
    taskEventEndTimesLog = taskEventEndTimesLog((lag+1):end);
    juiceOnTimesLog = juiceOnTimesLog((lag+1):end);
    juiceOffTimesLog = juiceOffTimesLog((lag+1):end);
    tmpEye = tmpEye((lag+1):end);
  elseif (lag < 0)
    %Means shifting back improves match,
    lag = abs(lag);
    taskEventIDsBlk = taskEventIDsBlk((lag+1):end);
    taskEventStartTimesBlk = taskEventStartTimesBlk((lag+1):end);
    taskEventEndTimesBlk = taskEventEndTimesBlk((lag+1):end);
    juiceOnTimesBlk = juiceOnTimesBlk((lag+1):end);
    juiceOffTimesBlk = juiceOffTimesBlk((lag+1):end);
    taskEventIDs = taskEventIDs((lag+1):end);
  end
  %Recalculate lengths
  LogLen = length(taskEventIDsLog);
  BlkLen = length(taskEventIDsBlk);
  if LogLen < BlkLen %monkeyLogic stopped storing things early.
    %Chop the Blackrock events down to the size of what Mkl has stored
    taskEventIDsBlk = taskEventIDsBlk(1:LogLen);
    if sum(taskEventIDsBlk == taskEventIDsLog) == length(taskEventIDsLog) %If They are now a match
      taskEventStartTimesBlk = taskEventStartTimesBlk(1:LogLen);
      taskEventEndTimesBlk = taskEventEndTimesBlk(1:LogLen);
      stimFramesLost = stimFramesLost(1:LogLen);
      errorArray = errorArray(1:LogLen); %Unique here because it is a stored variable from MKL.
      juiceOnTimesBlk = juiceOnTimesBlk(1:LogLen);
      juiceOffTimesBlk = juiceOffTimesBlk(1:LogLen);
      taskEventIDs = taskEventIDs(1:LogLen);
    else
      error("Simple alignment failed - Blackrock and MonkeyLogic logs still don't match")
    end
  else %Blackrock stopped early.
    %Blackrock shut off before monkeyLogic was done.
    taskEventIDsLog = taskEventIDsLog(1:BlkLen);
    if sum(taskEventIDsBlk == taskEventIDsLog) == length(taskEventIDsLog) %If They are now a match
      %do the same chopping to monkeyLogic related structures
      %Collect trialEventIDs.
      taskEventStartTimesLog = taskEventStartTimesLog(1:BlkLen);
      taskEventEndTimesLog = taskEventEndTimesLog(1:BlkLen);
      stimFramesLost = stimFramesLost(1:BlkLen);
      juiceOnTimesLog = juiceOnTimesLog(1:BlkLen);
      juiceOffTimesLog = juiceOffTimesLog(1:BlkLen);
      tmpEye = tmpEye(1:BlkLen);
    else
      error("Simple alignment failed - Blackrock and MonkeyLogic logs still don't match")
    end
  end
end
assert(sum(taskEventIDsLog == taskEventIDsBlk) == length(taskEventIDsLog), 'Blackrock and MonkeyLogic do not have matching Event numbers.')

%% The strobe is the most reliable marker that the stimuli has begun - use this for event start times.
%We now have the start of trials, Blackrock.Eventmarkers. We want to move that to Blackrock.Strobe, which we 
%believe is more accurate. This will be done by finding the transitions for
%the strobe from High (white) to low (black) closest to the time stamp.

%What is the discrepency?
offsets = zeros(sum(packetData == rewardMarker),1);

%Using Photodiode to define true start times.
if params.usePhotodiode
  taskEventStartTimesBlkPreStrobe = taskEventStartTimesBlk; %Saving these for Eye signal related processing.
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

%% In cases where the stimuli are .avi's, we should check for the appropriate frameMotion data (representing)
if any(strfind(translationTable{1},'.avi'))
  %assumes filename below sitting in directory with stimuli
  frameMotionFile = dir([params.stimDir '/**/frameMotion_complete.mat']);
  load([frameMotionFile(1).folder filesep frameMotionFile(1).name],'frameMotionData'); 
  %go through the translation table, comparing to frameMotionData.stimVid
  frameMotionNames = [{frameMotionData(:).stimVid}'];
  tmpFrameMotionData = struct('stimVid',[],'objNames',[],'objShapes',[],'objRadii',[],'vidB_XShift',[],'objLoc',[],'fps',[],'width',[],'height',[]);
  for table_i = 1:length(translationTable)
    dataInd = find(strcmp(frameMotionNames, translationTable{table_i}));
    if ~isempty(dataInd)
      tmpFrameMotionData(table_i) = frameMotionData(dataInd);
    else
      tmpFrameMotionData(table_i) = frameMotionData(1);
      tmpFrameMotionData(table_i).objNames = [];
      tmpFrameMotionData(table_i).objShapes = [];
      tmpFrameMotionData(table_i).objRadii = [];
      tmpFrameMotionData(table_i).objLoc = [];
      tmpFrameMotionData(table_i).vidB_XShift = [];
      tmpFrameMotionData(table_i).stimVid = translationTable{table_i};
    end
  end
end
%% Now, calculate stimulation log to ephys clock conversion 
  %Make the Model by comparing Blackrock event start times to monkeylogic.
  logVsBlkModel = fitlm(taskEventStartTimesLog, taskEventStartTimesBlkPreStrobe);
  disp(logVsBlkModel)

  taskEventStartTimesFit = predict(logVsBlkModel, taskEventStartTimesLog);

  %Shift over all the events, using the calculated line.
  eventTimeAdjustments = taskEventStartTimesFit-taskEventStartTimesBlkPreStrobe;
  disp(strcat('Max magnitude fit residual, msec: ',num2str(max(abs(eventTimeAdjustments)))));
  
  figure('Name','Sync adjustments from model fitting','NumberTitle','off')
  hist(eventTimeAdjustments,40);
  title('Sync adjustments from model fitting');
  xlabel('offset (adjusted - original) (ms)');
  ylabel('count');
  disp(median(eventTimeAdjustments));
  [~,adjSortInds] = sort(abs(eventTimeAdjustments-median(eventTimeAdjustments)), 'descend');
  disp('worst alignments, log file times');
  disp(taskEventStartTimesLog(adjSortInds(1:min(5,length(adjSortInds)))));
  disp('worst alignments, adjusted times');
  disp(taskEventStartTimesBlk(adjSortInds(1:min(5,length(adjSortInds)))));
  disp('worst alignments, adjustment values (ms)');
  disp(eventTimeAdjustments(adjSortInds(1:min(5,length(adjSortInds)))));
  
  %Assuming you're happy with this model, shift values only in the log to
  %Blackrock time.
  juiceOnTimesBlk = predict(logVsBlkModel, juiceOnTimesLog);
  juiceOffTimesBlk = predict(logVsBlkModel, juiceOffTimesLog);
  %Not doing this for event start times, since I trust the photodiode more
  %than the eventmarker as a metric of true stim start time.
  
  fprintf('average offset %s ms\n', num2str(mean(offsets)));
  fprintf('range of offset %d ms - %d ms \n', [min(offsets), max(offsets)])
%% Output
%Adding random numbers to these - they aren't relevant for my current task,
%nor are they directly recorded by MKL.
fixSpotFlashStartTimesBlk = taskEventStartTimesBlk(1);
fixSpotFlashEndTimesBlk = taskEventEndTimesBlk(1);
fixationInTimesBlk = taskEventStartTimesBlk(1);
fixationOutTimesBlk = taskEventEndTimesBlk(1);

% finally, build the output structure - This Structure is rebuilt
% selectively in the 
taskData.taskDataSummary.TrialRecord = TrialRecord;
taskData.errorArray = errorArray;
taskData.taskEventIDs = taskEventIDs;
taskData.translationTable = translationTable;
taskData.frameMotionData = tmpFrameMotionData';
%taskData.stimJumps = stimJumps;
taskData.stimFramesLost = stimFramesLost;
taskData.taskEventStartTimes = taskEventStartTimesBlk;
%taskData.taskEventStartTimesFit = taskEventStartTimesFit;
taskData.taskEventEndTimes = taskEventEndTimesBlk;
taskData.trialStartTimesMkl = mklTrialStarts';
taskData.logVsBlkModel = logVsBlkModel;
taskData.fixationInTimes = fixationInTimesBlk;
taskData.fixationOutTimes = fixationOutTimesBlk;
taskData.juiceOnTimes = juiceOnTimesBlk;
taskData.juiceOffTimes = juiceOffTimesBlk;
taskData.fixSpotFlashStartTimes = fixSpotFlashStartTimesBlk;
taskData.fixSpotFlashEndTimes = fixSpotFlashEndTimesBlk;
taskData.stimParams = 0;
taskData.RFmap = 0;
taskData.eyeData = tmpEye;
taskData.screenStats = screenStats;
taskData.eyeCal.PixelsPerDegree = MLConfig.PixelsPerDegree;
taskData.eyeCal.origin = MLConfig.EyeTransform{1,2}.origin;
taskData.eyeCal.gain = MLConfig.EyeTransform{1,2}.gain;

end
%

function translationTable = updateTranslationTable(translationTable)
%Function takes the translation table drawn from the MonkeyLogic file, and
%replaces entries with Old names with their new names.

RosettaStone = {...
{'monkeyObserving_1191.avi', 'monkeyObserving_1171.avi'}
{'monkeyObserving_1192.avi', 'monkeyObserving_1172.avi'}
{'monkeyObserving_1193.avi', 'monkeyObserving_1173.avi'}
{'monkeyObserving_1194.avi', 'monkeyObserving_1174.avi'}
{'monkeyObserving_1195.avi', 'monkeyObserving_1175.avi'}

{'monkeyGoalDir_1201.avi', 'monkeyGoalDir_1101.avi'}
{'monkeyGoalDir_1202.avi', 'monkeyGoalDir_1102.avi'}
{'monkeyGoalDir_1203.avi', 'monkeyGoalDir_1103.avi'}
{'monkeyGoalDir_1204.avi', 'monkeyGoalDir_1104.avi'}
{'monkeyGoalDir_1205.avi', 'monkeyGoalDir_1105.avi'}

{'humanGoalDir_1201.avi', 'humanGoalDir_1101.avi'}
{'humanGoalDir_1202.avi', 'humanGoalDir_1102.avi'}
{'humanGoalDir_1203.avi', 'humanGoalDir_1103.avi'}
{'humanGoalDir_1204.avi', 'humanGoalDir_1104.avi'}
{'humanGoalDir_1205.avi', 'humanGoalDir_1105.avi'}

{'Cut_obj_interact1_1.avi', 'objects_2101.avi'}; ...
{'Cut_obj_interact1_2.avi', 'objects_2102.avi'}; ...
{'Cut_Order1Interaction_1.avi', 'monkeyChasing_1111.avi'}; ...
{'Cut_Order1Interaction_2.avi', 'monkeyFighting_1121.avi'}; ...
{'Cut_Order1Interaction_3.avi', 'monkeyGrooming_1141.avi'}; ...
{'Cut_Order1Interaction_4.avi', 'monkeyGrooming_1142.avi'}; ...
{'Cut_Order1Interaction_5.avi', 'monkeyMounting_1131.avi'}; ...
{'Cut_Order2Interaction_1.avi', 'monkeyFighting_1124.avi'}; ...
{'Cut_Order2Interaction_2.avi', 'monkeyGrooming_1146.avi'}; ...
{'Cut_Order2Interaction_3.avi', 'monkeyMounting_1133.avi'}; ...
{'Cut_Order2Interaction_4.avi', 'monkeyGrooming_1147.avi'}; ...
{'Cut_Order2Interaction_5.avi', 'monkeyMounting_1136.avi'}; ...
{'Cut_Order2LandscapesFull_1.avi', 'landscape_4001.avi'}; ...
{'Cut_Order2LandscapesFull_2.avi', 'landscape_4002.avi'}; ...

{'Dephased_Order1Interactionchasing.avi', 'scramble_3001.avi'}; ...
{'Dephased_Order1Interactionfighting.avi', 'scramble_3002.avi'}; ...

{'Dephased_Order1Interaction_1.avi', 'scramble_3001.avi'}; ...
{'Dephased_Order1Interaction_2.avi', 'scramble_3002.avi'}; ...
{'Dephased_Order1Interaction_3.avi', 'scramble_3003.avi'}; ...
{'Dephased_Order1Interaction_4.avi', 'scramble_3004.avi'}; ...
{'Dephased_Order1Interaction_5.avi', 'scramble_3005.avi'}; ...

{'Dephased_Order2Interaction_1.avi', 'scramble_3006.avi'}; ...
{'Dephased_Order2Interaction_2.avi', 'scramble_3007.avi'}; ...
{'Dephased_Order2Interaction_3.avi', 'scramble_3008.avi'}; ...
{'Dephased_Order2Interaction_4.avi', 'scramble_3009.avi'}; ...
{'Dephased_Order2Interaction_5.avi', 'scramble_30010.avi'}; ...

{'Cut_human_interaction_1.avi', 'humanChasing_1111.avi'}; ...
{'Cut_human_interaction_2.avi', 'humanFollowing_1161.avi'}; ...
{'Cut_human_interaction_3.avi', 'humanGrooming_1141.avi'}; ...
{'Cut_human_interaction_4.avi', 'humanFighting_1121.avi'}; ...
{'Cut_human_interaction_5.avi', 'humanMounting_1131.avi'};...

{'Order1GoalLeft_1_Order1GoalRight_1.avi', 'monkeyGoalDir_1101.avi'};...
{'Order1GoalLeft_2_Order1GoalRight_2.avi', 'monkeyGoalDir_1102.avi'};...
{'Order1GoalLeft_3_Order1GoalRight_3.avi', 'monkeyGoalDir_1103.avi'};...
{'Order1GoalLeft_4_Order1GoalRight_4.avi', 'monkeyGoalDir_1104.avi'};...
{'Order1GoalLeft_5_Order1GoalRight_5.avi', 'monkeyGoalDir_1105.avi'};...

{'Order1MonkLeft_1_Order1MonkRight_1.avi', 'monkeyIdle_1301.avi'};...
{'Order1MonkLeft_2_Order1MonkRight_2.avi', 'monkeyIdle_1302.avi'};...
{'Order1MonkLeft_3_Order1MonkRight_3.avi', 'monkeyIdle_1303.avi'};...
{'Order1MonkLeft_4_Order1MonkRight_4.avi', 'monkeyIdle_1304.avi'};...
{'Order1MonkLeft_5_Order1MonkRight_5.avi', 'monkeyIdle_1305.avi'};...

{'human_alone1_1_1_human_alone2_1_1.avi', 'humanIdle_1301.avi'};...
{'human_alone1_1_2_human_alone2_1_2.avi', 'humanIdle_1302.avi'};...
{'human_alone1_1_3_human_alone2_1_3.avi', 'humanIdle_1303.avi'};...
{'human_alone1_1_4_human_alone2_1_4.avi', 'humanIdle_1304.avi'};...
{'human_alone1_1_5_human_alone2_1_5.avi', 'humanIdle_1305.avi'};...

{'human_goal1_1_4_human_goal2_1_4.avi', 'humanGoalDir_1101.avi'};
{'human_goal1_1_5_human_goal2_1_5.avi', 'humanGoalDir_1102.avi'};
{'human_goal1_1_1_human_goal2_1_1.avi', 'humanGoalDir_1103.avi'};
{'human_goal1_1_2_human_goal2_1_2.avi', 'humanGoalDir_1104.avi'};
{'human_goal1_1_3_human_goal2_1_3.avi', 'humanGoalDir_1105.avi'};

{'Cut_Order1Chasing1.avi', 'monkeyChasing_1113.avi'}; ...
{'Cut_Order1Chasing2.avi', 'monkeyChasing_1112.avi'}; ...
{'Cut_Order1Chasing3.avi', 'monkeyChasing_1111.avi'}; ...
{'Cut_Order1Chasing4.avi', 'monkeyChasing_1114.avi'}; ...
{'Cut_Order1Chasing5.avi', 'monkeyChasing_1115.avi'}; ...

{'Cut_Order1Fighting1.avi', 'monkeyFighting_1123.avi'}; ...
{'Cut_Order1Fighting2.avi', 'monkeyFighting_1122.avi'}; ...
{'Cut_Order1Fighting3.avi', 'monkeyFighting_1121.avi'}; ...
{'Cut_Order1Fighting4.avi', 'monkeyFighting_1124.avi'}; ...
{'Cut_Order1Fighting5.avi', 'monkeyFighting_1125.avi'}; ...

{'Cut_Order1Mounting1.avi', 'monkeyMounting_1131.avi'}; ...
{'Cut_Order1Mounting2.avi', 'monkeyMounting_1132.avi'}; ...
{'Cut_Order1Mounting3.avi', 'monkeyMounting_1133.avi'}; ...
{'Cut_Order1Mounting4.avi', 'monkeyMounting_1134.avi'}; ...
{'Cut_Order1Mounting5.avi', 'monkeyMounting_1135.avi'}; ...

};

oldHalf = cell(size(RosettaStone));
newHalf = cell(size(RosettaStone));

for event_i = 1:length(RosettaStone)
  oldHalf{event_i} = RosettaStone{event_i}{1};
  newHalf{event_i} = RosettaStone{event_i}{2};
end
 
% {'Cut_Order1Fighting1.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Fighting2.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Fighting3.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Fighting4.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Fighting5.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Grooming5.avi', 'socialInteraction', 'grooming'}
% 
% {'Dephased_Order1Interactionchasing.avi', 'scramble', 'highMotion'}; ...
% {'Dephased_Order1Interactionfighting.avi', 'scramble', 'highMotion'}; ...
% };

for ii = 1:length(translationTable)
  if any(strcmp(translationTable{ii}, oldHalf))
    translationTable{ii} = newHalf{strcmp(translationTable{ii}, oldHalf)};
  end
end

end

