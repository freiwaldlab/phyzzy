function [ analysisParamFilename ] = buildAnalysisParamFileSocialMovies( )    
%buildAnalysisParamFile saves a mat file of parameters, which control the
%behavior of analyzeSession
%   todo: option to load 'fixed' params from file, for ease accross days


%%%%%%%  USER PARAMETERS, EDIT ROUTINELY %%%%%%%%
runNum = '010';
dateSubject = '180221ALAN'; 
machine = 'laptop';

switch machine
  case 'rig'
    ephysVolume = '/Volumes/Users-1/User/Desktop';
    stimulusLogVolume = '/Volumes/Users/FreiwaldLab/Desktop';
    outputVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Analyzed';
    stimParamsFilename = '/Users/stephenserene/Desktop/Freiwald/AnalysisSerene/StimParamsFullSocialInt.mat';   %#ok
  case 'laptop'
    ephysVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Blackrock'; 
    stimulusLogVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Visiko';
    outputVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Analyzed';
    stimParamsFilename = '/Users/stephenserene/Desktop/Freiwald/AnalysisSerene/StimParamsFullSocialInt.mat';   %#ok
  case 'hopper'
    ephysVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Blackrock'; 
    stimulusLogVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Visiko';
    outputVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Analyzed';
    stimParamsFilename = '/Freiwald/sserene/ephys/AnalysisSerene/StimParamsFullSocialInt.mat';   %#ok
  case 'turing'
    ephysVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Blackrock'; 
    stimulusLogVolume = '/Freiwald/ephys/sserene/ALAN_DATA/Visiko';
    outputVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Analyzed';
    stimParamsFilename = '/Freiwald/sserene/ephys/AnalysisSerene/StimParamsFullSocialInt.mat';   %#ok
end
analysisLabel = 'Basic';
analysisParamFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'
preprocessedDataFilenameStem = 'preprocessedData.mat';
saveFig = 1;           %#ok
closeFig = 0;          %#ok
exportFig = 1;         %#ok
saveFigData = 0;       %#ok
savePreprocessed = 0;  %#ok
verbosity = 'INFO'; %other options, 'DEBUG', 'VERBOSE';


% parameters preprocessSpikes and preprocessLFP, see functions for details
ephysParams.needLFP = 1;
ephysParams.needSpikes = 1;
ephysParams.spikeChannels = [1]; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
ephysParams.lfpChannels = [1]; 
ephysParams.channelNames = {'ML'};
ephysParams.lfpChannelScaleBy = [8191/32764, 8191/32764, 8191/32764]; %converts raw values to microvolts
ephysParams.commonRef = [0 0 0]; %not yet implemented; will allow software re-refrence across headstages
ephysParams.stimulationChannels = []; %not yet implemented; will read stimulation currents recorded at headstage
ephysParams.cPtCal = 1/30; % conversion from spike sample indices to timestep of decimated LFP
ephysParams.decimateFactorPass1 = 6; %note: product of the two decimate factors should be 30, if 1 khz samples desired
ephysParams.decimateFactorPass2 = 5;
ephysParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)
%note: use Blackrock indexing for unitsToUnsort and unitsToDiscard, so unsorted is 0, first defined unit is 1, etc.
ephysParams.unitsToUnsort = {[],[],[]}; %these units will be re-grouped with u0
ephysParams.unitsToDiscard = {[],[],[]}; %these units will be considered noise and discarded
ephysParams.spikeWaveformPca = 0;
ephysParams.plotSpikeWaveforms = 0; %0, 1 to build then close, 2 to build and leave open
ephysParams.shiftSpikeWaveforms = 0;
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
hp1Hz = designfilt('highpassiir', 'FilterOrder',8,'PassbandFrequency',1, ...
  'StopbandAttenuation',100,'PassbandRipple',0.5,'SampleRate',1000);     %#ok
% note: with these specifications, returns a 48th order butterworth filter
butter1Hz200Hz_v1 = designfilt('bandpassiir','DesignMethod','butter','PassbandFrequency1',1,'PassbandFrequency2',200,...
  'SampleRate',1000,'MatchExactly','passband','StopbandFrequency1',0.67,'StopbandFrequency2',250);
[tmp1,tmp2] = butter(4,[1/500,200/500]);
butter1Hz200Hz_v2 = [tmp1,tmp2];        %#ok
ephysParams.filter = butter1Hz200Hz_v1; % if filtering desired, ephysFilter is a digitalFilter
ephysParams.plotFilterResult = 0; %#ok

% parameters preprocessAnalogIn, see function for details
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
analogInParams.plotFilterResult = 1; %#ok

photodiodeParams.needPhotodiode = 0;
photodiodeParams.frameRate = 85;
photodiodeParams.levelCalibType = 'autoAndCheck';
photodiodeParams.numLevels = 1;
photodiodeParams.saveCalibFile = 0;
photodiodeParams.centerCornerOffset = 5.3;
photodiodeParams.frameTriggerChannel = 129; 
photodiodeParams.stimulusTriggerChannel = []; %#ok

% parameters preprocessLogFile, see function for details
stimSyncParams.logProcessor = @preprocessLogFileVisikoMovie;
stimSyncParams.tryPreparsedLogFile = 1;
stimSyncParams.keepTriggersAndSubTriggers = 0;
stimSyncParams.subTriggerArrayFilenames = {'socialSceneConcatSubTriggers.mat'};
stimSyncParams.usePhotodiode = 0;        %#ok
%
eyeCalParams.needEyeCal = 0;
eyeCalParams.method = 'hardcodeZero'; %'zeroEachFixation'
eyeCalParams.makePlots = 1;
eyeCalParams.eyeXChannelInd = 1;
eyeCalParams.eyeYChannelInd = 2;
eyeCalParams.eyeDChannelInd = 3;
eyeCalParams.gainX = 112;
eyeCalParams.gainY = 107;
eyeCalParams.flipX = 1;
eyeCalParams.flipY = 1; 
eyeCalParams.offsetX = -6.4;
eyeCalParams.offsetY = -5.6; 
eyeCalParams.calFile = ''; %note: needed only when method = fromFile
eyeCalParams.fixOutLag = 10; 
eyeCalParams.minFixZeroTime = 1000; %#ok

accelParams.needAccelCal = 0;
accelParams.accelChannels = {[4;5;6]};
accelParams.channelGains = {[1/.666 1/.666 1/.666]};
accelParams.calMethods = {'hardcode'}; %other option is 'calFile'; calibration method
accelParams.calFiles = {''}; %if method is 'calFile', an ns2 filename

% parameters for excludeStimuli, see function for details
excludeStimParams.needExcludeTrials = 0;
excludeStimParams.fixPre = 100; %ms
excludeStimParams.fixPost = 100; %ms
excludeStimParams.flashPre = 0;  %ms
excludeStimParams.flashPost = 0; %ms
excludeStimParams.juicePre = 0; % optional, ms
excludeStimParams.juicePost = 0; % optional, ms
excludeStimParams.DEBUG = 0; % makes exclusion criterion plots if true
% additional optional excludeStimParams: accel1, accel2, minStimDur (ms)

psthParams.psthPre = 300; % if e.g. +200, then start psth 200ms before trial onset; 
psthParams.psthImDur = 2480;  % only need to set this for variable length stim runs; else, comes from log file
psthParams.psthPost = 300;
psthParams.smoothingWidth = 100;  %psth smoothing width, in ms
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
tfParams.movingWin = [650 5]; %was [200 5]; 
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

plotSwitch.imagePsth = 0;
plotSwitch.categoryPsth = 0;
plotSwitch.prefImRaster = 0;
plotSwitch.prefImRasterEvokedOverlay = 0;
plotSwitch.prefImMultiChRasterEvokedOverlay = 0;
plotSwitch.imageTuningSorted = 0;
plotSwitch.stimPrefBarPlot = 0;
plotSwitch.stimPrefBarPlotEarly = 0;
plotSwitch.stimPrefBarPlotLate = 0;
plotSwitch.tuningCurves = 0;
plotSwitch.tuningCurvesEarly = 0;
plotSwitch.tuningCurvesLate = 0;
plotSwitch.RF = 1;
plotSwitch.rfEarly = 0;
plotSwitch.rfLate = 0;
plotSwitch.latencyRF = 0;
plotSwitch.evokedPowerRF = 0;
plotSwitch.evokedPsthMuaMultiCh = 0;
plotSwitch.evokedByCategory = 0;
plotSwitch.analogInByItem = 0;
plotSwitch.analogInDerivativesByItem = 0;
plotSwitch.colorPsthEvoked = 1;
plotSwitch.linePsthEvoked = 0;
plotSwitch.colorPsthItemMarginals = 1;
plotSwitch.runSummary = 0;
plotSwitch.runSummaryImMeanSub = 0;
plotSwitch.runSummaryImMeanSubDiv = 0;
plotSwitch.lfpPowerMuaScatter = 0; 
plotSwitch.lfpPeakToPeakMuaScatter = 0;
plotSwitch.lfpLatencyMuaLatency = 0;
plotSwitch.lfpPowerAcrossChannels = 0;
plotSwitch.lfpPeakToPeakAcrossChannels = 0;
plotSwitch.lfpLatencyShiftAcrossChannels = 0;
plotSwitch.singleTrialLfpByCategory = 0;
plotSwitch.lfpSpectraByCategory = 0;
plotSwitch.spikeSpectraByCategory = 0;
plotSwitch.SpikeSpectraTfByImage = 0;
plotSwitch.lfpSpectraTfByImage = 0;
plotSwitch.couplingPhasesUnwrapped = 0;
plotSwitch.couplingPhasesAsOffsets = 0;
plotSwitch.couplingPhasesPolar = 0;
plotSwitch.tfSpectraByCategory = 0;
plotSwitch.tfErrs = 0;           %#ok

%%%% note: all analysisGroups cell arrays are nx1, NOT 1xn
analysisGroups.selectivityIndex.groups = {{'face';'nonface'},{'face';'object'},{'face';'body'}};
%
analysisGroups.stimPrefBarPlot.groups = {{{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'};{'face';'object';'body'}}};
analysisGroups.stimPrefBarPlot.colors  = {{{'b';'c';'y';'g';'m';'r';'k'};{'b';'g';'r'}}};
analysisGroups.stimPrefBarPlot.names = {'fobPlus'};
analysisGroups.stimPrefBarPlot.groupDepth = 2;
%
analysisGroups.stimulusLabelGroups.groups = {{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'}};
analysisGroups.stimulusLabelGroups.names = {'fobPlus'};
analysisGroups.stimulusLabelGroups.colors = {{'b';'c';'y';'g';'m';'r';'k'}};
%
analysisGroups.evokedPotentials.groups = {{'socialInteraction';'objects';'fighting'}};
analysisGroups.evokedPotentials.names = {'socVobj'};
analysisGroups.evokedPotentials.colors = {{'b';'g';'c'}};
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
analysisGroups.colorPsthEvoked.groups = {{'socialInteraction';'objects';'fighting'}};
analysisGroups.colorPsthEvoked.names = {'socVobj'};
analysisGroups.colorPsthEvoked.colors = {{'b';'g';'c'}};
%
analysisGroups.linePsthEvoked.groups = {{'socialInteraction';'objects'}};
analysisGroups.linePsthEvoked.names = {'socVobj'};
analysisGroups.linePsthEvoked.colors = {{'b';'g'}};
%
analysisGroups.colorPsthItemMarginals.groups = {{'fighting';'mounting';'grooming';'chasing';'objects';'scramble'}};
analysisGroups.colorPsthItemMarginals.names = {'socVids'};
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
analysisGroups.tfCouplingByCategory.groups = {{'face'};{'nonface'};{'object'};{'body'}}; %#ok
%%%%%

calcSwitch.categoryPSTH = 1;
calcSwitch.imagePSTH = 1;
calcSwitch.faceSelectIndex = 1;
calcSwitch.faceSelectIndexEarly = 0;
calcSwitch.faceSelectIndexLate = 0;
calcSwitch.inducedTrialMagnitudeCorrection = 0;
calcSwitch.evokedSpectra = 1;
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
lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);        
spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
taskFilename = sprintf('%s/%s/%s0%s.log',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
photodiodeFilename = lfpFilename;                                                           %#ok
outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
analysisParamFilename = strcat(outDir,analysisParamFilenameStem);
preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok
%
load('cocode2.mat');
psthColormap = map;  %#ok
%
if ~exist(outDir,'dir')
  mkdir(outDir);
end
save(analysisParamFilename);
end

