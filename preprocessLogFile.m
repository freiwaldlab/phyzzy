function [ taskData, stimTiming ] = preprocessLogFile(logFilename, taskTriggers, diodeTriggers, params )
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

disp(logFilename);

if isfield(params,'logProcessor')
  [taskData, stimTiming] = params.logProcessor(logFilename, taskTriggers, diodeTriggers, params);
  return
end

if params.usePhotodiode
  error('photodiode synchronization not enabled');
end

disp('parsing serial IO packets');
packetTimes = taskTriggers.TimeStampSec;
packetData = dec2bin(taskTriggers.UnparsedData);
% note: bit (end-5) is 1 when fixating; bit (end) toggles on each stimulus start
stimBits = packetData(:,end);
% note: the following two lines assume that the first packet is stimulus,
% not fixation in
taskEventStartTimesBlk = 1000*packetTimes(vertcat(1,diff(stimBits)) ~= 0);  %Blk affix signifies Blackrock reference frame
disp('number of stim triggers received by blackrock');
disp(length(taskEventStartTimesBlk));
clear packetData; 

% parse log file
disp('Loading visiko log file and converting to matlab structure');
assert(logical(exist(logFilename,'file')),'The stimulus log file you requested does not exist.');
logStruct = xml2struct(logFilename);
assert(isfield(logStruct.VISIKOLOG,'EndStimulation'),'Error: stimulation end not included in log file');
assert(strcmp(logStruct.VISIKOLOG.Attributes.tasktype,'Bitmap Continuous'),'Error: unknown Visiko task type. Must be Bitmap Continuous');
if isa(logStruct.VISIKOLOG.DOCDATA,'cell')
  s = input(sprintf('Visiko parameters changed during task; %d parameter sets found. Enter the number to analyze, or n to quit: ',length(logStruct.VISIKOLOG.DOCDATA)),'s');
  if strcmp(s,'n')
    return;
  end
  logStruct.VISIKOLOG.DOCDATA = logStruct.VISIKOLOG.DOCDATA{str2double(s)};
end
stimTiming.shortest = 1000*str2double(logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT.pictureTimeFrom.Text); 
stimTiming.longest = 1000*str2double(logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT.pictureTimeTo.Text);
stimTiming.ISI = 1000*str2double(logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT.pauseTime.Text);
stimParams = logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT;
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

for i = 1:length(logStruct.VISIKOLOG.Object)
  stimulusStruct = logStruct.VISIKOLOG.Object{i};
  if isfield(stimulusStruct,'PictureCompleted')
    continue
  end
  tmp = regexp(stimulusStruct.Pictures.Picture.Filename.Text,'\','split');
  taskEventIDs = vertcat(taskEventIDs, strtrim(tmp{end})); %note: for some reason, filenames have trailing whitespace; trim it off
  stimFramesLost(i) = str2double(stimulusStruct.Frameslost.Text);
  stimJumps(i,:) = [str2double(stimulusStruct.Pictures.Picture.Jump.x.Text), str2double(stimulusStruct.Pictures.Picture.Jump.y.Text)];
  taskEventStartTimesLog(i) = str2double(stimulusStruct.Start.Text);
  taskEventEndTimesLog(i) = str2double(stimulusStruct.End.Text);
end
Output.VERBOSE(sprintf('number of stimulus trials: %s',length(taskEventIDs)));
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
  assert(length(taskEventStartTimesBlk) - length(taskEventStartTimesLog) <= 1, 'Error: Start triggers missing from log file');
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
taskData.stimParams = stimParams;
taskData.RFmap = str2double(logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT.activateJumping.Text);
% if RF mapping task, convert the picture positions to x-y degrees from fovea
if taskData.RFmap
  % in case these become useful sometime...
%   gridSizeX = str2double(logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT.jumpExtentHorizontal.Text); % in degrees
%   gridSizeY = str2double(logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT.jumpExtentVertical.Text); % in degrees
%   gridSpacing = str2double(logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT.jumpGridSize.Text); % in degrees
  gridPointsPixX = unique(taskData.stimJumps(:,1)); % in pixels; note unique() returns in low->high sorted order
  gridPointsPixY = unique(taskData.stimJumps(:,2)); % in pixels; note unique() returns in low->high sorted order
%calculate screen diagonal length in pixels for pixel/degree conversion
  screenWidthPix = str2double(logStruct.VISIKOLOG.SESSIONDATA.DisplayMode.Width.Text);
  screenHeightPix = str2double(logStruct.VISIKOLOG.SESSIONDATA.DisplayMode.Height.Text);
  screenDpix = sqrt(screenWidthPix^2 + screenHeightPix^2);
  screenDcm = str2double(logStruct.VISIKOLOG.SESSIONDATA.Monitor.DiagonalSize.Text);
  screenEyeDistance = str2double(logStruct.VISIKOLOG.SESSIONDATA.Monitor.Distance.Text);
  pixelConversionCoef = screenDcm/(screenDpix*screenEyeDistance);
  taskData.gridPointsDegX = atan(gridPointsPixX*pixelConversionCoef)*180/pi; %pix to cm, then rad to deg
  taskData.gridPointsDegY = atan(gridPointsPixY*pixelConversionCoef)*180/pi; %pix to cm, then rad to deg
  % now, convert each picture jump entry
  taskData.stimJumps = atan(taskData.stimJumps*pixelConversionCoef)*180/pi;
end
end
%
