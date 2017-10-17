function [ output_args ] = checkAnalysisParamFile( analysisParamFilename )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
load(analysisParamFilename);
logString = 'Assigned default values to: \n';

% runNum
assert(logical(exist('runNum','var')),'Invalid analysis parameter file: must specify a run number as runNum.');
assert(ischar(runNum),'Invalid analysis parameter file: runNum must be a string'); 

% dateSubject
assert(logical(exist('dateSubject','var')),'Invalid analysis parameter file: must specify a date and subject as dateSubject.');
assert(ischar(dateSubject),'Invalid analysis parameter file: dateSubject must be a string'); 

% ephysVolume
assert(logical(exist('ephysVolume','var')),'Invalid analysis parameter file: must specify the ephys data folder as ephysVolume.');
assert(ischar(ephysVolume),'Invalid analysis parameter file: ephysVolume must be a string'); 

% stimulusLogVolume
assert(logical(exist('stimulusLogVolume','var')),'Invalid analysis parameter file: must specify the stimulus log folder as stimulusLogVolume.');
assert(ischar(stimulusLogVolume),'Invalid analysis parameter file: stimulusLogVolume must be a string');

% outputVolume
assert(logical(exist('runNum','var')),'Invalid analysis parameter file: must specify a run number as runNum.');
assert(ischar(runNum),'Invalid analysis parameter file: runNum must be a string'); 

% stimParamsFilename
assert(logical(exist('outputVolume','var')),'Invalid analysis parameter file: must specify an output folder as outputVolume.');
assert(ischar(outputVolume),'Invalid analysis parameter file: outputVolume must be a string');

% stimParamsFilename
assert(logical(exist('outputVolume','var')),'Invalid analysis parameter file: must specify an output folder as outputVolume.');
assert(ischar(outputVolume),'Invalid analysis parameter file: outputVolume must be a string');

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
if ~logical(exist('verbosity','var')) || ~ischar(verbostity) || ~ismember(verbosity,{'INFO','DEBUG','VERBOSE'})
  verbosity = 'INFO'; %#ok
  logString = strcat(logString,'verbosity\n');
end

%%% EPHYS PARAMS %%%
if ~logical(exist('ephysParams','var')) || ~isstruct(ephysParams)
  ephysParams.needLFP = 0;
  ephysParams.needSpikes = 0;
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
      assert(length(channelNames) == numChannels,'Invalid analysis parameter file: channelNames length must match length of channels to analyze, or be zero.')
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
  assert(isfield(ephysParams,'cPtCal') && isnumeric(ephysParams.),'Invalid analysis parameter file: must specify conversion from spike time units to decimated LFP indices');
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
  else
    if ephysParams.samPerMS ~= 1
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
    ephysParams.needSpikes = 0;
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
    ephysParams.shiftSpikeWaveforms = 0;
    logString = strcat(logString,'ephysParams.plotFilterResult\n');
  end
end

%%% AnalogInParams %%%
analogInParams.needAnalogIn = 0;
analogInParams.analogInChannels = [138,139,140]; 
analogInParams.channelNames = {'eyeX','eyeY','eyeD'};
analogInParams.analogInChannelScaleBy = [5/32764 5/32764 5/32764]; %converts raw values to volts
analogInParams.decimateFactorPass1 = 1; 
analogInParams.decimateFactorPass2 = 1;
analogInParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to analogIn (should be raw rate/productOfDecimateFactors)
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
butter200Hz_v1 = designfilt('lowpassiir', 'PassbandFrequency', 120, 'StopbandFrequency', 480, 'PassbandRipple', 1,...
  'StopbandAttenuation', 60, 'SampleRate', 1000, 'MatchExactly', 'passband');  %this returns a 3rd order iir filter
analogInParams.filters = {0,0,0};%{butter200Hz_v1;butter200Hz_v1;butter200Hz_v1}; %filter channel i if filters{i} is digital filter or 1x2 numeric array
analogInParams.plotFilteredSignal = 1; %#ok

photodiodeParams.needPhotodiode = 0;
photodiodeParams.channels = [1;2]; %#ok

% parameters preprocessLogFile, see function for details
stimSyncParams.usePhotodiode = 0;        %#ok
%
eyeCalParams.needEyeCal = 0;
eyeCalParams.method = 'zeroEachFixation';
eyeCalParams.makePlots = 0;
eyeCalParams.eyeXChannelInd = 1;
eyeCalParams.eyeYChannelInd = 2;
eyeCalParams.eyeDChannelInd = 3;
eyeCalParams.gainX = 112;
eyeCalParams.gainY = 107;
eyeCalParams.flipX = 1;
eyeCalParams.flipY = 1; 
eyeCalParams.offsetX = -6.4;
eyeCalParams.offsetY = -5.6; 
eyeCalParams.minFixZeroTime = 1000; %#ok

accelParams.needAccelCal = 0;
accelParams.accelChannels = {[4;5;6]};
accelParams.channelGains = {[1/.666 1/.666 1/.666]};
accelParams.calMethods = {'hardcode'}; %other option is 'calFile'; calibration method
accelParams.calFiles = {''}; %if method is 'calFile', an ns2 filename

% parameters for excludeStimuli, see function for details
excludeStimParams.fixPre = 100; %ms
excludeStimParams.fixPost = 100; %ms
excludeStimParams.flashPre = 0;  %ms
excludeStimParams.flashPost = 0; %ms
excludeStimParams.juicePre = 0; % optional, ms
excludeStimParams.juicePost = 0; % optional, ms
excludeStimParams.DEBUG = 0; % makes exclusion criterion plots if true
% additional optional excludeStimParams: accel1, accel2, minStimDur (ms)

psthParams.psthPre = 100; % if e.g. +200, then start psth 200ms before trial onset; 
psthParams.psthImDur = 0;  % only need to set this for variable length stim runs; else, comes from log file
psthParams.psthPost = 300;
psthParams.smoothingWidth = 10;  %psth smoothing width, in ms
psthParams.errorType = 1; %chronux convention: 1 is poisson, 2 is trialwise bootstrap, 3 is across trial std for binned spikes, bootstrap for spike times 
psthParams.errorRangeZ = 1; %how many standard errors to show
psthParams.bootstrapSamples = 100;
psthParams.psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'


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

correlParams.maxShift = []; % a number, or empty
correlParams.matchTimeRanges = 1;
correlParams.timeDifferenceBound = [0,200];
correlParams.normalize = 1;
correlParams.useJacknife = 0;
correlParams.jacknifeDraws = 100;
switch machine
  case 'laptop'
    correlParams.jacknifeParallelWorkers = 0;    
  case 'hopper'
    correlParams.jacknifeParallelWorkers = 0;   
  case 'turing'
    correlParams.jacknifeParallelWorkers = 20;    
end
spikeCorrelSmoothingWidth = 5; %ms
filterPoints = -20*spikeCorrelSmoothingWidth:20*spikeCorrelSmoothingWidth;
smoothingFilter = exp(-1*filterPoints.^2/(2*spikeCorrelSmoothingWidth^2));
correlParams.smoothingFilter = smoothingFilter/sum(smoothingFilter); %#ok
%
lfpAlignParams.samPerMS = 1; % because this is after decimation
lfpAlignParams.msPreAlign = psthParams.psthPre+tfParams.movingWin(1)/2; 
lfpAlignParams.msPostAlign = psthParams.psthImDur+psthParams.psthPost+tfParams.movingWin(1)/2;
%
spikeAlignParams.preAlign = psthParams.psthPre+3*psthParams.smoothingWidth;
spikeAlignParams.postAlign = psthParams.psthImDur+psthParams.psthPost+3*psthParams.smoothingWidth;   %#ok
% for lfps, constrain first and (optional) last [n m] samples to 0 mean
useDCSUB = 0;
if useDCSUB
  %lfpAlignParams.DCSUB_SAM = [lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10; 0, 0 ]; % 0-th order 
  lfpAlignParams.DCSUB_SAM = [lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10;lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10 ]; % 1st order 
else
  lfpAlignParams.DCSUB_SAM = 0;   %#ok
end
% firing rate calculation epochs. Can provide either time (ms from stim onset),
% or function handle, which will receive the minimum stimulus duration in the run as an input
frEpochsCell = {{60, @(stimDur) stimDur+60};...
                {60, 260}; ...
                {260, @(stimDur) stimDur+60}}; %#ok


% Boolean variables to specify which computations to perform; TODO: read
% from config file, eventually with conditional on log file info
calcCoherenceRFcpt = 0;  %#ok
calcCoherenceRFcc = 0;   %#ok
calcCoherenceRFptpt = 0; %#ok
calcGrangerRF = 0;       %#ok 



%%%% note: all analysisGroups cell arrays are nx1, NOT 1xn
analysisGroups.selectivityIndex.groups = {{'face';'nonface'},{'face';'object'},{'face';'body'}};
%
analysisGroups.stimPrefBarPlot.groups = {{{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'};{'face';'object';'body'}}};
analysisGroups.stimPrefBarPlot.colors  = {{{'b';'c';'y';'g';'m';'r';'k'};{'b';'g';'r'}}};
analysisGroups.stimPrefBarPlot.names = {'fobPlus'};
%
analysisGroups.stimulusLabelGroups.groups = {{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'}};
analysisGroups.stimulusLabelGroups.names = {'fobPlus'};
analysisGroups.stimulusLabelGroups.colors = {{'b';'c';'y';'g';'m';'r';'k'}};
%
analysisGroups.evokedPotentials.groups = {{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'}};
analysisGroups.evokedPotentials.names = {'fobPlus'};
analysisGroups.evokedPotentials.colors = {{'b';'c';'y';'g';'m';'r';'k'}};
%
analysisGroups.analogInPotentials.groups = {{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'}};
analysisGroups.analogInPotentials.channels = {[1; 2]};
analysisGroups.analogInPotentials.names = {'eyePositions,fobPlus'};
analysisGroups.analogInPotentials.units = {'degrees visual angle'};
analysisGroups.analogInPotentials.colors = {{'b';'c';'y';'g';'m';'r';'k'}};
%
analysisGroups.analogInDerivatives.groups = {{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'}};
analysisGroups.analogInDerivatives.channels = {[1; 2]};
analysisGroups.analogInDerivatives.names = {'eyeVelocity,fobPlus'};
analysisGroups.analogInDerivatives.units = {'degrees visual angle/sec'};
analysisGroups.analogInDerivatives.colors = {{'b';'c';'y';'g';'m';'r';'k'}};
%
analysisGroups.colorPsthEvoked.groups = {{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'}};
analysisGroups.colorPsthEvoked.names = {'fobPlus'};
analysisGroups.colorPsthEvoked.colors = {{'b';'c';'y';'g';'m';'r';'k'}};
%
analysisGroups.linePsthEvoked.groups = {{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'}};
analysisGroups.linePsthEvoked.names = {'fobPlus'};
analysisGroups.linePsthEvoked.colors = {{'b';'c';'y';'g';'m';'r';'k'}};
%
analysisGroups.evokedPsthOnePane.groups = {{'face';'nonface'}};
analysisGroups.evokedPsthOnePane.names = {'faceVnon'};
%
analysisGroups.tuningCurves.groups = {{'humanFaceL90','humanFaceL45','humanFaceFront','humanFaceR45','humanFaceR90'},...
  {'monkeyFaceL90','monkeyFaceL45','monkeyFaceFront','monkeyFaceR45','monkeyFaceR90'}}; %can be images or categories
analysisGroups.tuningCurves.paramValues = {[-90 -45 0 45 90], [-90 -45 0 45 90]};
analysisGroups.tuningCurves.paramLabels = {'viewing angle (degrees)','viewing angle (degrees)'};
analysisGroups.tuningCurves.names = {'Human face view','Monkey face view'};
%
analysisGroups.spectraByCategory.groups = {{'face';'nonface'}};  %todo: add spectra diff?
analysisGroups.spectraByCategory.names = {'faceVnon'};
analysisGroups.spectraByCategory.colors = {{'r';'b'}};
%
analysisGroups.tfSpectraByCategory.groups = {{'face'};{'nonface'}};%{'object'};{'body'}      %todo: add tf spectra diff?
analysisGroups.tfSpectraByCategory.names = {'face','nonface'};%'nonface';'object';'body'
%
analysisGroups.lfpSingleTrialsByCategory.groups = {{'face';'nonface'}};
analysisGroups.lfpSingleTrialsByCategory.names = {'faceVnon'};
%
analysisGroups.coherenceByCategory.groups = {{'face';'nonface'}}; %{'face';'object';'body'};{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'hand';'techno'}
analysisGroups.coherenceByCategory.colors = {{'r';'b'}}; %{'r';'g';'b'};{'b';'c';'y';'g';'m';'r';'k';'k'}
analysisGroups.coherenceByCategory.names = {'faceVnon'}; %'fob';'slimCats'
%
analysisGroups.tfCouplingByCategory.groups = {{'face'};{'nonface'};{'object'};{'body'}};

analysisGroups.byImage = {};      %#ok
analysisGroupColors.byImage = {}; %#ok
%%%%%

calcSwitch.categoryPSTH = 1;
calcSwitch.imagePSTH = 1;
calcSwitch.faceSelectIndex = 1;
calcSwitch.faceSelectIndexEarly = 1;
calcSwitch.faceSelectIndexLate = 1;
calcSwitch.inducedTrialMagnitudeCorrection = 0;
calcSwitch.evokedSpectra = 0;
calcSwitch.inducedSpectra = 0;
calcSwitch.evokedImageTF = 0;
calcSwitch.inducedImageTF = 0;
calcSwitch.evokedCatTF = 0;
calcSwitch.inducedCatTF = 0;
calcSwitch.meanEvokedTF = 0;
calcSwitch.trialMeanSpectra = 0;
calcSwitch.coherenceByCategory = 0;
calcSwitch.spikeTimes = 0;
calcSwitch.useJacknife = 0;      

if calcSwitch.useJacknife
  chronuxParams.err(1) = 2; %#ok
end

%%% set paths and directories, EDIT RARELY %%%
analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);   %#ok
lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);        %#ok
spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
taskFilename = sprintf('%s/%s/%s0%s.log',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
analysisParamsFilename = strcat(outDir,analysisParamsFilenameStem);
preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok
%
load('cocode2.mat');
psthColormap = map;  %#ok
%
if ~exist(outDir,'dir')
  mkdir(outDir);
end
save(analysisParamsFilename);

end

