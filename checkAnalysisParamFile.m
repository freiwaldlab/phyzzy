function [ ] = checkAnalysisParamFile( analysisParamFilename )
%checkAnalysisParamFile defines, checks, and applies the parameter file specification
%   - For variables with a default value, applies the default if no or invalid value given
%   - For variables with no default, throws informative error if no or invalid value given
% 
load(analysisParamFilename);
logString = 'Function checkAnalysisParamFile assigned default values to: \n';

% runNum
assert(logical(exist('runNum','var')),'Invalid analysis parameter file: must specify a run number as runNum.');
assert(ischar(runNum),'Invalid analysis parameter file: runNum must be a string'); 

% dateSubject
assert(logical(exist('dateSubject','var')),'Invalid analysis parameter file: must specify a date and subject as dateSubject.');
assert(ischar(dateSubject),'Invalid analysis parameter file: dateSubject must be a string'); 

% lfpFilename
assert(logical(exist('lfpFilename','var')),'Invalid analysis parameter file: must specify the lfp data file as lfpFilename.');
assert(ischar(lfpFilename),'Invalid analysis parameter file: lfpFilename must be a string'); 

% spikeFilename
assert(logical(exist('spikeFilename','var')),'Invalid analysis parameter file: must specify the spike data file as spikeFilename.');
assert(ischar(spikeFilename),'Invalid analysis parameter file: spikeFilename must be a string'); 

% taskLogFilename
assert(logical(exist('taskFilename','var')),'Invalid analysis parameter file: must specify the task log file as taskFilename.');
assert(ischar(taskFilename),'Invalid analysis parameter file: taskFilename must be a string');

% outputVolume
assert(logical(exist('outputVolume','var')),'Invalid analysis parameter file: must specify an output volume as outputVolume.');
assert(ischar(outputVolume),'Invalid analysis parameter file: outputVolume must be a string'); 

% stimParamsFilename
assert(logical(exist('stimParamsFilename','var')),'Invalid analysis parameter file: must specify a stimulus/task log filename as stimParamsFilename.');
assert(ischar(stimParamsFilename),'Invalid analysis parameter file: stimParamsFilename must be a string');

%todo: add check for preprocessed data filename? only relevant if saving...

%%% SAVE SWITCHES %%%
% saveFig
if ~logical(exist('saveFig','var')) || ~ismember(saveFig,[0 1])
  saveFig = 1; %#ok
  logString = strcat(logString,'saveFig\n');
end

% closeFig
if ~logical(exist('closeFig','var')) || ~ismember(closeFig,[0 1])
  closeFig = 0; %#ok
  logString = strcat(logString,'closeFig\n');
end

% exportFig
if ~logical(exist('exportFig','var')) || ~ismember(exportFig,[0 1])
  exportFig = 0; %#ok
  logString = strcat(logString,'exportFig\n');
end

% saveFigData
if ~logical(exist('saveFigData','var')) || ~ismember(saveFigData,[0 1])
  saveFigData = 0; %#ok
  logString = strcat(logString,'saveFigData\n');
end

% savePreprocessed
if ~logical(exist('savePreprocessed','var')) || ~ismember(savePreprocessed,[0 1])
  savePreprocessed = 0; %#ok
  logString = strcat(logString,'savePreprocessed\n');
end

%%% VERBOSITY %%%
if ~logical(exist('verbosity','var')) || ~ischar(verbosity) || ~ismember(verbosity,{'INFO','DEBUG','VERBOSE'})
  verbosity = 'INFO'; %#ok
  logString = strcat(logString,'verbosity\n');
end

%%% EPHYS PARAMS %%%
if ~logical(exist('ephysParams','var')) || ~isstruct(ephysParams)
  ephysParams.needLFP = 0;
  ephysParams.needSpikes = 0;
  ephysParams.spikeChannels = [];
  ephysParams.lfpChannels = [];
  ephysParams.channelNames = {};
  logString = strcat(logString,'ephysParams\n');
else
  if ~isfield(ephysParams,'needLFP') || ~ismember(ephysParams.needLFP,[0 1])
    ephysParams.needLFP = 0;
    logString = strcat(logString,'ephysParams.needLFP\n');
  end
  if ~isfield(ephysParams,'needSpikes') || ~ismember(ephysParams.needSpikes,[0 1])
    ephysParams.needSpikes = 0;
    logString = strcat(logString,'ephysParams.needSpikes\n');
  end
  if ~isfield(ephysParams,'spikeChannels') || ~isnumeric(ephysParams.spikeChannels)
    ephysParams.spikeChannels = [];
    logString = strcat(logString,'ephysParams.spikeChannels\n');
  end
  if ~isfield(ephysParams,'lfpChannels') || ~isnumeric(ephysParams.lfpChannels)
    ephysParams.lfpChannels = [];
    logString = strcat(logString,'ephysParams.lfpChannels\n');
  end
  if ~isfield(ephysParams,'commonRef') || ~isnumeric(ephysParams.commonRef)
    ephysParams.commonRef = zeros(size(ephysParams.lfpChannels));
    logString = strcat(logString,'ephysParams.commonRef\n');
  end
  if ephysParams.needLFP && ephysParams.needSpikes
    assert(length(ephysParams.spikeChannels) == length(ephysParams.lfpChannels),...
      'Invalid analysis parameter file: spikeChannels and lfpChannels must be the same length if analyzing both');
    assert(all(ephysParams.spikeChannels == ephysParams.lfpChannels),...
      'Invalid analysis parameter file: spikeChannels and lfpChannels must be in the same order, if analyzing both');
  end
  % channelNames
  if (ephysParams.needLFP || ephysParams.needSpikes)
    numChannels = max(ephysParams.needLFP*length(ephysParams.lfpChannels),ephysParams.needSpikes*length(ephysParams.spikeChannels));
    if isfield(ephysParams,'channelNames')
      assert(length(ephysParams.channelNames) == numChannels,'Invalid analysis parameter file: channelNames length must match length of channels to analyze, or be zero.')
    else
      ephysParams.channelNames = cell(numChannels,1);
      for channel_i = 1:numChannels
        ephysParams.channelNames{channel_i} = sprintf('ch%d',channel_i);
      end
      logString = strcat(logString,'ephysParams.channelNames\n');
    end
  end
  % lfpChannelScaleBy
  if ~isfield(ephysParams,'lfpChannelScaleBy') || ~isnumeric(ephysParams.lfpChannelScaleBy)
    ephysParams.lfpChannelScaleBy = ones(size(ephysParams.lfpChannels));
    logString = strcat(logString,'ephysParams.lfpChannelScaleBy (set to one)\n');
  end
  if length(ephysParams.lfpChannelScaleBy) == 1 && length(ephysParams.lfpChannels) > 1
    ephysParams.lfpChannelScaleBy = ephysParams.lfpChannelScaleBy*ones(size(ephysParams.lfpChannels));
    logString = strcat(logString,'ephysParams.lfpChannelScaleBy (applied single given value to all channels)\n');
  end
  assert(isfield(ephysParams,'cPtCal') && isnumeric(ephysParams.cPtCal),'Invalid analysis parameter file: must specify conversion from spike time units to decimated LFP indices');
  if ~isfield(ephysParams,'decimateFactorPass1') || ~isnumeric(ephysParams.decimateFactorPass1)
    ephysParams.decimateFactorPass1 = 6;
    logString = strcat(logString,'ephysParams.decimateFactorPass1\n');
  end
  if ~isfield(ephysParams,'decimateFactorPass2') || ~isnumeric(ephysParams.decimateFactorPass2)
    ephysParams.decimateFactorPass2 = 5;
    logString = strcat(logString,'ephysParams.decimateFactorPass2\n');
  end
  if ~isfield(ephysParams,'samPerMS') || ~isnumeric(ephysParams.samPerMS)
    ephysParams.samPerMS = 1;
    logString = strcat(logString,'ephysParams.samPerMS\n');
  elseif ephysParams.samPerMS ~= 1
    disp('Warning: lfp samples per ms ~= 1. In current implementation, this will lead to axes in samples incorrectly labeled as in ms');
  end
end

if ~isfield(ephysParams,'unitsToUnsort') || ~iscell(ephysParams.unitsToUnsort) || ~length(ephysParams.unitsToUnsort) == length(ephysParams.spikeChannels)
  ephysParams.unitsToUnsort = cell(length(ephysParams.spikeChannels));
  logString = strcat(logString,'ephysParams.unitsToUnsort\n');
end
if ~isfield(ephysParams,'unitsToDiscard') || ~iscell(ephysParams.unitsToDiscard) || ~length(ephysParams.unitsToDiscard) == length(ephysParams.spikeChannels)
  ephysParams.unitsToDiscard = cell(length(ephysParams.spikeChannels));
  logString = strcat(logString,'ephysParams.unitsToDiscard\n');
end
if ~isfield(ephysParams,'spikeWaveformPca') || ~ismember(ephysParams.spikeWaveformPca,[0 1])
  ephysParams.spikeWaveformPca = 0;
  logString = strcat(logString,'ephysParams.spikeWaveformPca\n');
end
if ~isfield(ephysParams,'plotSpikeWaveforms') || ~ismember(ephysParams.plotSpikeWaveforms,[0 1])
  ephysParams.plotSpikeWaveforms = 0;
  logString = strcat(logString,'ephysParams.plotSpikeWaveforms\n');
end
if ~isfield(ephysParams,'shiftSpikeWaveforms') || ~ismember(ephysParams.shiftSpikeWaveforms,[0 1])
  ephysParams.shiftSpikeWaveforms = 0;
  logString = strcat(logString,'ephysParams.shiftSpikeWaveforms\n');
end
if ~isfield(ephysParams,'filter') %todo: add check for digitalFilter object type?
  ephysParams.filter = '';
  logString = strcat(logString,'ephysParams.filter\n');
end
if ~isfield(ephysParams,'plotFilterResult') || ~ismember(ephysParams.plotFilterResult,[0 1])
  ephysParams.plotFilterResult = 0;
  logString = strcat(logString,'ephysParams.plotFilterResult\n');
end


%%% AnalogInParams %%%
if ~logical(exist('analogInParams','var')) || ~isstruct(analogInParams)
  analogInParams.needAnalogIn = 0;
  analogInParams.analogInChannels = [];
  analogInParams.channelNames = {};
  logString = strcat(logString,'analogInParams\n');
else
  if ~isfield(analogInParams,'needAnalogIn') || ~ismember(analogInParams.needAnalogIn,[0 1])
    analogInParams.needAnalogIn = 0;
    logString = strcat(logString,'analogInParams.needAnalogIn\n');
  end
  if ~isfield(analogInParams,'analogInChannels') || ~isnumeric(analogInParams.analogInChannels)
    analogInParams.analogInChannels = [];
    logString = strcat(logString,'analogInParams.analogInChannels\n');
  end
  assert(isfield(analogInParams,'channelNames') && length(analogInParams.channelNames) == length(analogInParams.analogInChannels),...
    'Invalid analysis parameter file: analogInParams.channelNames length must match length of analogIn channels to analyze.');
  if ~isfield(analogInParams,'analogInChannelScaleBy') || ~isnumeric(analogInParams.analogInChannelScaleBy)
    analogInParams.analogInChannelScaleBy = ones(size(analogInParams.analogInChannels));
    logString = strcat(logString,'analogInParams.analogInChannelScaleBy (set to one)\n');
  end
  if ~isfield(analogInParams,'decimateFactorPass1') || ~isnumeric(analogInParams.decimateFactorPass1)
    analogInParams.decimateFactorPass1 = 1;
    logString = strcat(logString,'analogInParams.decimateFactorPass1\n');
  end
  if ~isfield(analogInParams,'decimateFactorPass2') || ~isnumeric(analogInParams.decimateFactorPass2)
    analogInParams.decimateFactorPass2 = 1;
    logString = strcat(logString,'analogInParams.decimateFactorPass2\n');
  end
  if ~isfield(analogInParams,'samPerMS') || ~isnumeric(analogInParams.samPerMS)
    analogInParams.samPerMS = 1;
    logString = strcat(logString,'analogInParams.samPerMS\n'); 
  else
    if analogInParams.samPerMS ~= 1
      disp('Warning: analog input samples per ms ~= 1. In current implementation, this will lead to axes in samples incorrectly labeled as in ms');
    end
  end
  if ~isfield(analogInParams,'filters') || ~(length(analogInParams.filters) == length(analogInParams.analogInChannels))  %todo: add check for digitalFilter object type?
    analogInParams.filters = cell(length(analogInParams.analogInChannels),1);
    logString = strcat(logString,'analogInParams.filters\n');
  end
  if ~isfield(analogInParams,'plotFilterResult') || ~ismember(analogInParams.plotFilterResult,[0 1])
    analogInParams.plotFilterResult = 0;
    logString = strcat(logString,'analogInParams.plotFilterResult\n');
  end
end

%%% photodiodeParams %%%
if ~logical(exist('photodiodeParams','var')) || ~isstruct(photodiodeParams)
  photodiodeParams.needPhotodiode = 0;
  photodiodeFilename = '';
  logString = strcat(logString,'photodiodeParams\n');
else
  if ~isfield(photodiodeParams,'needPhotodiode') || ~ismember(photodiodeParams.needPhotodiode,[0 1])
    photodiodeParams.needPhotodiode = 0;
    logString = strcat(logString,'photodiodeParams.needPhotodiode\n');
  end
  if photodiodeParams.needPhotodiode
    assert(logical(exist('photodiodeFilename','var')),'Invalid analysis parameter file: must specify an photodiode file as photodiodeFilename if photodiodeParams.needPhotodiode.');
    assert(ischar(photodiodeFilename),'Invalid analysis parameter file: photodiodeFilename must be a string');
  end
  if photodiodeParams.needPhotodiode
    assert(isnumeric(photodiodeParams.frameTriggerChannel));
    assert(isnumeric(photodiodeParams.stimulusTriggerChannel));
  end
end

%%%
if ~logical(exist('stimSyncParams','var')) || ~isstruct(stimSyncParams)
  stimSyncParams.usePhotodiode = 0;
  logString = strcat(logString,'stimSyncParams\n');
else
  if ~isfield(stimSyncParams,'usePhotodiode') || ~ismember(stimSyncParams.usePhotodiode,[0 1])
    stimSyncParams.usePhotodiode = 0;
    logString = strcat(logString,'stimSyncParams.usePhotodiode\n');
  end
end

%%%
if ~logical(exist('eyeCalParams','var')) || ~isstruct(eyeCalParams)
  eyeCalParams.needEyeCal = 0;
  logString = strcat(logString,'eyeCalParams\n');
else
  if ~isfield(eyeCalParams,'needEyeCal') || ~ismember(eyeCalParams.needEyeCal,[0 1])
    eyeCalParams.needEyeCal = 0;
    logString = strcat(logString,'eyeCalParams.needEyeCal\n');
  end
  if ~isfield(eyeCalParams,'method') || ~ischar(eyeCalParams.method) || ~ismember(eyeCalParams.method,{'autoZeroSingle','zeroEachFixation','hardcodeZero','fromFile','monkeyLogic'})
    eyeCalParams.method = 'autoZeroSingle';
    logString = strcat(logString,'eyeCalParams.method\n');
  end
  if ~isfield(eyeCalParams,'makePlots') || ~ismember(eyeCalParams.needEyeCal,[0 1])
    eyeCalParams.needEyeCal = 0;
    logString = strcat(logString,'eyeCalParams.makePlots\n');
  end
  assert(isfield(eyeCalParams,'eyeXChannelInd'),'Invalid analysis parameter file: if needEyeCal, must supply eye X input channel index as eyeCalParams.eyeXChannelInd.')
  assert(isnumeric(eyeCalParams.eyeXChannelInd),'Invalid analysis parameter file: eyeCalParams.eyeXChannelInd must be numeric')
  assert(isfield(eyeCalParams,'eyeYChannelInd'),'Invalid analysis parameter file: if needEyeCal, must supply eye Y input channel index as eyeCalParams.eyeYChannelInd.')
  assert(isnumeric(eyeCalParams.eyeYChannelInd),'Invalid analysis parameter file: eyeCalParams.eyeYChannelInd must be numeric')
  % todo: make eyeD optional? Might not be supplied by all eye trackers...
  assert(isfield(eyeCalParams,'eyeDChannelInd'),'Invalid analysis parameter file: if needEyeCal, must supply eye diameter input channel index as eyeCalParams.eyeDChannelInd.')
  assert(isnumeric(eyeCalParams.eyeDChannelInd),'Invalid analysis parameter file: eyeCalParams.eyeDChannelInd must be numeric') %todo: implement flip with negative gain?
  assert(isfield(eyeCalParams,'gainX'),'Invalid analysis parameter file: if needEyeCal, must supply X gain as eyeCalParams.gainX.')
  assert(isnumeric(eyeCalParams.gainX),'Invalid analysis parameter file: eyeCalParams.gainX must be numeric')
  assert(isfield(eyeCalParams,'gainY'),'Invalid analysis parameter file: if needEyeCal, must supply Y gain as eyeCalParams.gainY.')
  assert(isnumeric(eyeCalParams.gainY),'Invalid analysis parameter file: eyeCalParams.gainY must be numeric')
  if ~isfield(eyeCalParams,'flipX') || ~ismember(eyeCalParams.flipX,[0 1])
    eyeCalParams.flipX = 0;
    logString = strcat(logString,'eyeCalParams.flipX\n');
  end
  if ~isfield(eyeCalParams,'flipY') || ~ismember(eyeCalParams.flipY,[0 1])
    eyeCalParams.flipY = 0;
    logString = strcat(logString,'eyeCalParams.flipY\n');
  end
  if strcmp(eyeCalParams.method,'hardcodeZero')
    assert(isfield(eyeCalParams,'offsetX'),'Invalid analysis parameter file: if needEyeCal and eye cal method is hardcodeZero, must supply X offset as eyeCalParams.offsetX.');
    assert(isnumeric(eyeCalParams.offsetX),'Invalid analysis parameter file: eyeCalParams.offsetX must be numeric');
    assert(isfield(eyeCalParams,'offsetY'),'Invalid analysis parameter file: if needEyeCal and eye cal method is hardcodeZero, must supply Y offset as eyeCalParams.offsetY.');
    assert(isnumeric(eyeCalParams.offsetY),'Invalid analysis parameter file: eyeCalParams.offsetY must be numeric');
  end
  if strcmp(eyeCalParams.method,'zeroEachFixation')
    assert(isfield(eyeCalParams,'minFixZeroTime'),...
      'Invalid analysis parameter file: if needEyeCal and eye cal method is zeroEachFixation, must supply min fixation time for re-zero as eyeCalParams.minFixZeroTime.');
    assert(isnumeric(eyeCalParams.minFixZeroTime),'Invalid analysis parameter file: eyeCalParams.minFixZeroTime must be numeric');
  end
  if strcmp(eyeCalParams.method,'fromFile')
    assert(isfield(eyeCalParams,'calFile') && ~isempty(isfield(eyeCalParams.calFile)),...
      'Invalid analysis parameter file: if needEyeCal and eye cal method is fromFile, must supply calibration file name as eyeCalParams.calFile.');
  end
end
%%%
if ~logical(exist('accelParams','var')) || ~isstruct(accelParams)
  accelParams.needAccelCal = 0;
  logString = strcat(logString,'accelParams\n');
else
  if ~isfield(accelParams,'needAccelCal') || ~ismember(accelParams.needAccelCal,[0 1])
    accelParams.needAccelCal = 0;
    logString = strcat(logString,'accelParams.needAccelCal\n');
  end
  if ~isfield(accelParams,'accelChannels') || ~iscell(accelParams.accelChannels)
    accelParams.accelChannels = {};
    logString = strcat(logString,'accelParams.accelChannels\n');
  end
  for accel_i = 1:length(accelParams.accelChannels)
    if (~accelParams.needAccelCal) || isempty(accelParams.accelChannels)
      break
    end
    assert(isfield(accelParams,'calMethods') && ismember(accelParams.calMethod,{'hardcode','calFile'}),...
      'Invalid analysis parameter file: if needAccelCal, must specify a method as hardcode or calFile, for each accelerometer.');
    if strcmp(accelParams.calMethod,'hardcode')
      if isempty(accelParams.channelGains{accel_i})
        accelParams.channelGains{accel_i} = ones(size(accelParams.accelChannels{accel_i}));
        logString = strcat(logString,'accelParams.channelGains\n');
      else
        assert(length(accelParams.accelChannels{accel_i}) == length(accelParams.channelGains{accel_i}),...
          'Invalid analysis parameter file: if supplying accelerometer gains, must provide one for each channel on a given accelerometer');
      end
    else
      assert(isfield(accelParams,'calFiles') && length(accelParams.calFiles) == length(accelParams.accelChannels) && ischar(accelParams.calFiles{accel_i}),...
        'Invalid analysis parameter file: if using a cal file for any accelerometer, must provide a cell array wtih entries for each accelerometer, and the cal filename must be a string.');
    end
  end
end
%%%
if ~logical(exist('excludeStimParams','var')) || ~isstruct(excludeStimParams)
  excludeStimParams.fixPre = 0; 
  excludeStimParams.fixPost = 0; 
  excludeStimParams.flashPre = 0;  
  excludeStimParams.flashPost = 0; 
  excludeStimParams.juicePre = 0; 
  excludeStimParams.juicePost = 0; 
  excludeStimParams.DEBUG = 0;
  logString = strcat(logString,'excludeStimParams\n');
else
  if ~isfield(excludeStimParams,'fixPre') || ~isnumeric(excludeStimParams.fixPre)
    excludeStimParams.fixPre = 0;
    logString = strcat(logString,'excludeStimParams.fixPre\n');
  end
  if ~isfield(excludeStimParams,'fixPost') || ~isnumeric(excludeStimParams.fixPost)
    excludeStimParams.fixPost = 0;
    logString = strcat(logString,'excludeStimParams.fixPost\n');
  end
  if ~isfield(excludeStimParams,'flashPre') || ~isnumeric(excludeStimParams.flashPre)
    excludeStimParams.flashPre = 0;
    logString = strcat(logString,'excludeStimParams.flashPre\n');
  end
  if ~isfield(excludeStimParams,'flashPost') || ~isnumeric(excludeStimParams.flashPost)
    excludeStimParams.flashPost = 0;
    logString = strcat(logString,'excludeStimParams.flashPost\n');
  end
  if ~isfield(excludeStimParams,'DEBUG') || ~ismember(excludeStimParams.DEBUG,[0,1])
    excludeStimParams.DEBUG = 0;
    logString = strcat(logString,'excludeStimParams.DEBUG\n');
  end
end
%%%
if ~logical(exist('psthParams','var')) || ~isstruct(psthParams)
  psthParams.psthPre = 100; % if e.g. +200, then start psth 200ms before trial onset; 
  psthParams.psthImDur = 0;  % only need to set this for variable length stim runs; else, comes from log file
  psthParams.psthPost = 300;
  psthParams.smoothingWidth = 10;  %psth smoothing width, in ms
  psthParams.errorType = 1; %chronux convention: 1 is poisson, 2 is trialwise bootstrap, 3 is across trial std for binned spikes, bootstrap for spike times 
  psthParams.errorRangeZ = 1; %how many standard errors to show
  psthParams.psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'
  psthParams.bootstrapSamples = 100;
  logString = strcat(logString,'psthParams\n');
else
  if ~isfield(psthParams,'psthPre') || ~isnumeric(psthParams.psthPre)
    psthParams.psthPre = 100;
    logString = strcat(logString,'psthParams.psthPre\n');
  end
  if ~isfield(psthParams,'psthImDur') || ~isnumeric(psthParams.psthImDur)
    psthParams.psthImDur = 0; % only need to set this for variable length stim runs; else, comes from log file
    logString = strcat(logString,'psthParams.psthImDur\n');
  end
  if ~isfield(psthParams,'psthPost') || ~isnumeric(psthParams.psthPost)
    psthParams.psthPost = 100;
    logString = strcat(logString,'psthParams.psthPost\n');
  end
  if ~isfield(psthParams,'smoothingWidth') || ~isnumeric(psthParams.smoothingWidth)
    psthParams.smoothingWidth = 10;
    logString = strcat(logString,'psthParams.smoothingWidth\n');
  end
  if ~isfield(psthParams,'errorType') || ~ismember(psthParams.errorType,[1,2,3])
    psthParams.errorType = 1;
    logString = strcat(logString,'psthParams.psthErrorType\n');
  end
  if ~isfield(psthParams,'errorRangeZ') || ~isnumeric(psthParams.errorRangeZ)
    psthParams.errorRangeZ = 1;
    logString = strcat(logString,'psthParams.errorRangeZ\n');
  end
  if (~isfield(psthParams,'bootstrapSamples') || ~isnumeric(psthParams.bootstrapSamples)) && psthParams.errorType > 1
    psthParams.bootstrapSamples = 100;
    logString = strcat(logString,'psthParams.bootstrapSamples\n');
  end
  if ~isfield(psthParams,'psthColormapFilename') || ~ischar(psthParams.psthColormapFilename)
    psthParams.psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'
    logString = strcat(logString,'psthParams.psthParams.psthColormapFilename\n');
  end
end
%
if ~exist('psthColormap','var')  % todo: integrate with psthParams
  load(psthParams.psthColormapFilename);
  psthColormap = map;
end
%
%%%
if ~logical(exist('chronuxParams','var')) || ~isstruct(chronuxParams)
  % TW=3 with T=.2, then W = 15 Hz (5 tapers)
  % TW=1.5 with T=.1, then W = 15 Hz (2 tapers)
  % TW = 1.5 with T=.2, then W = 7.5 Hz (2 tapers)
  chronuxParams.tapers = [3 5]; %[3 5] is chronux default; 
  chronuxParams.pad = 1;
  chronuxParams.fs = 1;
  chronuxParams.trialave = 1;
  chronuxParams.err = [1 .05];  %note: first entry will be automatically switched to 2 if calcSwitch.useJacknife == 1
  chronuxParams.fpass = [0 .1]; 
  tfParams.movingWin = [200 5]; 
  tfParams.specgramRowAve = 0;
  logString = strcat(logString,'chronuxParams\n');
else
  if ~isfield(chronuxParams,'tapers') || ~isnumeric(chronuxParams.tapers) || ~length(chronuxParams.tapers) == 2 
    chronuxParams.tapers = [3 5];
    logString = strcat(logString,'chronuxParams.tapers\n');
  end
  if ~isfield(chronuxParams,'pad') || ~isnumeric(chronuxParams.pad) || chronuxParams.pad < -1
    chronuxParams.pad = 1;
    logString = strcat(logString,'chronuxParams.pad\n');
  end
  if ~isfield(chronuxParams,'fs') || ~isnumeric(chronuxParams.fs)
    chronuxParams.fs = 1;
    logString = strcat(logString,'chronuxParams.fs\n');
  end
  if ~isfield(chronuxParams,'trialave') || ~isnumeric(chronuxParams.trialave)  %note lower case 'a', to match chronux convention
    chronuxParams.trialave = 1;
    logString = strcat(logString,'chronuxParams.trialave\n');
  end
  if ~isfield(chronuxParams,'err') || ~isnumeric(chronuxParams.err) || ~length(chronuxParams.err) == 2 
    chronuxParams.err = [1 .05];
    logString = strcat(logString,'chronuxParams.err\n');
  end
  if ~isfield(chronuxParams,'fpass') || ~isnumeric(chronuxParams.fpass) || ~length(chronuxParams.fpass) == 2 
    chronuxParams.fpass = [0 .1];
    logString = strcat(logString,'chronuxParams.fpass\n');
  end
end
%%%
if~logical(exist('tfParams','var')) || ~isstruct(tfParams)
  tfParams.movingWin = [200 5]; 
  tfParams.specgramRowAve = 0;
  logString = strcat(logString,'tfParams\n');
else
  if ~isfield(tfParams,'movingWin') || ~isnumeric(tfParams.movingWin) || ~length(tfParams.movingWin) == 2 
    tfParams.movingWin = [200 5];
    logString = strcat(logString,'tfParams.movingWin\n');
  end
  if ~isfield(tfParams,'movingWin') || ~ismember(tfParams.specgramRowAve,[0,1])
    tfParams.specgramRowAve = 0;
    logString = strcat(logString,'tfParams.specgramRowAve\n');
  end
end
%%%
if~logical(exist('correlParams','var')) || ~isstruct(correlParams)
  correlParams.maxShift = []; % a number, or empty
  correlParams.matchTimeRanges = 1;
  correlParams.timeDifferenceBound = [0,200];
  correlParams.normalize = 1;
  correlParams.useJacknife = 0;
  correlParams.jacknifeDraws = 100;
  correlParams.jacknifeParallelWorkers = 0;   
  logString = strcat(logString,'correlParams\n');
else
  if ~isfield(correlParams,'maxShift') || ~(isnumeric(correlParams.maxShift) || isempty(correlParams.maxShift)) 
    correlParams.maxShift = [];
    logString = strcat(logString,'correlParams.maxShift\n');
  end
  if ~isfield(correlParams,'matchTimeRanges') || ~ismember(correlParams.matchTimeRanges,[0,1])
    correlParams.matchTimeRanges = 1;
    logString = strcat(logString,'correlParams.matchTimeRanges\n');
  end
  if ~isfield(correlParams,'timeDifferenceBound') || ~isnumeric(correlParams.timeDifferenceBound) || ~length(correlParams.timeDifferenceBound) == 2
    correlParams.timeDifferenceBound = [0,200]; %todo choose better defaults
    logString = strcat(logString,'correlParams.timeDifferenceBound\n');
  end
  if ~isfield(correlParams,'useJacknife') || ~ismember(correlParams.useJacknife,[0,1])
    correlParams.useJacknife = 0;
    logString = strcat(logString,'correlParams.useJacknife\n');
  end
  if ~isfield(correlParams,'jacknifeDraws') || ~isnumeric(correlParams.jacknifeDraws) 
    correlParams.jacknifeDraws = 100;
    logString = strcat(logString,'correlParams.jacknifeDraws\n');
  end
  if ~isfield(correlParams,'jacknifeParallelWorkers') || ~isnumeric(correlParams.jacknifeParallelWorkers) 
    correlParams.jacknifeParallelWorkers = 0;
    logString = strcat(logString,'correlParams.jacknifeParallelWorkers\n');
  end
  if ~isfield(correlParams,'smoothingFilter') || ~isnumeric(correlParams.smoothingFilter) 
    spikeCorrelSmoothingWidth = 5; %ms
    filterPoints = -20*spikeCorrelSmoothingWidth:20*spikeCorrelSmoothingWidth;
    smoothingFilter = exp(-1*filterPoints.^2/(2*spikeCorrelSmoothingWidth^2));
    correlParams.smoothingFilter = smoothingFilter/sum(smoothingFilter);
    logString = strcat(logString,'correlParams.smoothingFilter\n');
  end
end
%%%
if~logical(exist('lfpAlignParams','var')) || ~isstruct(lfpAlignParams)
  lfpAlignParams.samPerMS = 1; % because this is after decimation
  lfpAlignParams.msPreAlign = psthParams.psthPre+tfParams.movingWin(1)/2; 
  lfpAlignParams.msPostAlign = psthParams.psthImDur+psthParams.psthPost+tfParams.movingWin(1)/2;
  lfpAlignParams.DCSUB_SAM = 0;
  logString = strcat(logString,'lfpAlignParams\n');
else
  if ~isfield(lfpAlignParams,'samPerMS') || ~isnumeric(lfpAlignParams.samPerMS) 
    lfpAlignParams.samPerMS = 1;
    logString = strcat(logString,'lfpAlignParams.samPerMS\n');
  end
  if ~isfield(lfpAlignParams,'msPreAlign') || ~isnumeric(lfpAlignParams.msPreAlign) 
    lfpAlignParams.msPreAlign = psthParams.psthPre+tfParams.movingWin(1)/2; ;
    logString = strcat(logString,'lfpAlignParams.msPreAlign\n');
  end
  if ~isfield(lfpAlignParams,'msPostAlign') || ~isnumeric(lfpAlignParams.msPostAlign) 
    lfpAlignParams.msPostAlign = psthParams.psthImDur+psthParams.psthPost+tfParams.movingWin(1)/2;;
    logString = strcat(logString,'lfpAlignParams.msPostAlign\n');
  end
  if ~isfield(lfpAlignParams,'DCSUB_SAM') || ~isnumeric(lfpAlignParams.DCSUB_SAM) 
    lfpAlignParams.DCSUB_SAM = 0;
    logString = strcat(logString,'lfpAlignParams.DCSUB_SAM\n');
  end
  
  
end
%%%
if~logical(exist('spikeAlignParams','var')) || ~isstruct(spikeAlignParams)
  spikeAlignParams.preAlign = psthParams.psthPre+3*psthParams.smoothingWidth;
  spikeAlignParams.postAlign = psthParams.psthImDur+psthParams.psthPost+3*psthParams.smoothingWidth;
  logString = strcat(logString,'spikeAlignParams\n');
else
  if ~isfield(spikeAlignParams,'preAlign') || ~isnumeric(spikeAlignParams.preAlign) 
    spikeAlignParams.preAlign = psthParams.psthPre+3*psthParams.smoothingWidth;
    logString = strcat(logString,'spikeAlignParams.preAlign\n');
  end
  if ~isfield(spikeAlignParams,'postAlign') || ~isnumeric(spikeAlignParams.postAlign) 
    spikeAlignParams.postAlign = psthParams.psthImDur+psthParams.psthPost+3*psthParams.smoothingWidth;
    logString = strcat(logString,'spikeAlignParams.postAlign\n');
  end
end
%%%
if~logical(exist('frEpochsCell','var')) || ~iscell(frEpochsCell)  %todo: improve variable name
  % firing rate calculation epochs. Can provide either time (ms from stim onset),
  % or function handle, which will receive the minimum stimulus duration in the run as an input
  frEpochsCell = {{60, @(stimDur) stimDur+60}};
  for epoch_i = 1:length(frEpochsCell)
    assert(length(frEpochsCell{epoch_i}) == 2, 'Invalid analysis parameter file: frEpochsCell entries must have length 2.');
  end
end
%%%
% note: runAnalyses will likely use fields of an analysisGroups structure,
% and will catch errors and define defaults internally
if ~exist('analysisGroups','var') || ~isstruct(analysisGroups) 
  analysisGroups = struct();
  logString = strcat(logString,'analysisGroups\n');
end
%%%
% note: runAnalyses will likely have the option to use other fields of
% plotSwitch, and will catch errors and define defaults internally
if ~exist('plotSwitch','var') || ~isstruct(plotSwitch) 
  plotSwitch = struct();
  logString = strcat(logString,'plotSwitch\n');
end
%%%
% note: runAnalyses will likely have the option to use other fields of
% calcSwitch, and will catch errors and define defaults internally
if ~exist('calcSwitch','var') || ~isstruct(calcSwitch) 
  calcSwitch.spikeTimes = 0;
  logString = strcat(logString,'calcSwitch\n');
else
  if ~isfield(calcSwitch,'spikeTimes') || ~ismember(calcSwitch.spikeTimes,[0,1]) 
    calcSwitch.spikeTimes = 0;
    logString = strcat(logString,'calcSwitch.spikeTimes\n');
  end
end
if ~strcmp(logString,'checkAnalysisParamFile assigned default values to: \n')
  analysisParamFilenameParts = strsplit(analysisParamFilename, '.');
  analysisParamFilenameStem = analysisParamFilenameParts{1};
  movefile(analysisParamFilename, strcat(analysisParamFilenameStem,'Uncorrected.mat'));
  save(analysisParamFilename);
  fprintf(logString);
end
end

