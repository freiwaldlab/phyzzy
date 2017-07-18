function [ analysisParamsFilename ] = buildAnalysisParamFile( )    
%buildAnalysisParamFile saves a mat file of parameters, which control the
%behavior of analyzeSession
%   todo: option to load 'fixed' params from file, for ease accross days


%%%%%%%  USER PARAMETERS, EDIT ROUTINELY %%%%%%%%
runNum = '002';
dateSubject = '170405ALAN'; 
machine = 'laptop';

switch machine
  case 'rig'
    ephysVolume = '/Volumes/Users-1/User/Desktop';
    stimulusLogVolume = '/Volumes/Users/FreiwaldLab/Desktop';
    outputVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Analyzed';
    stimParamsFilename = '/Users/stephenserene/Desktop/Freiwald/AnalysisSerene/StimParamsFullFOB3.mat';   %#ok
  case 'laptop'
    ephysVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Blackrock'; 
    stimulusLogVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Visiko';
    outputVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Analyzed';
    stimParamsFilename = '/Users/stephenserene/Desktop/Freiwald/AnalysisSerene/StimParamsFullFOB3.mat';   %#ok
  case 'hopper'
    ephysVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Blackrock'; 
    stimulusLogVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Visiko';
    outputVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Analyzed';
    stimParamsFilename = '/Freiwald/sserene/ephys/AnalysisSerene/StimParamsFullFOB3.mat';   %#ok
  case 'turing'
    ephysVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Blackrock'; 
    stimulusLogVolume = '/Freiwald/ephys/sserene/ALAN_DATA/Visiko';
    outputVolume = '/Freiwald/sserene/ephys/ALAN_DATA/Analyzed';
    stimParamsFilename = '/Freiwald/sserene/ephys/AnalysisSerene/StimParamsFullFOB3.mat';   %#ok
end
analysisLabel = 'Basic';
analysisParamsFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'
preprocessedDataFilenameStem = 'preprocessedData.mat';
%categoryListSlim = {'humanFace','monkeyFace','place','fruit','humanBody','monkeyBody','techno'}; %minimal cat list for clean plots
saveFig = 1;           %#ok
exportFig = 0;         %#ok
saveFigData = 0;       %#ok
savePreprocessed = 0;  %#ok
verbosity = 'INFO'; %other options, 'DEBUG', 'VERBOSE';


% parameters preprocessSpikes and preprocessLFP, see functions for details
ephysParams.needLFP = 1;
ephysParams.needSpikes = 1;
ephysParams.spikeChannels = [1,35,33]; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
ephysParams.lfpChannels = [1,35,33]; 
ephysParams.channelNames = {'ML','AL','AM'};
ephysParams.unitsToIgnore = {{},{},{}}; %note: use Blackrock indexing, so unsorted is 0, first defined unit is 1, etc.
ephysParams.lfpChannelScaleBy = [8191/32764 8191/32764 8191/32764]; %converts raw values to microvolts
ephysParams.common_ref = [0, 35, 35]; %not yet implemented; will allow software re-refrence across headstages
ephysParams.stimulationChannels = []; %not yet implemented; will read stimulation currents recorded at headstage
ephysParams.cPtCal = 1/30; % conversion from spike sample indices to timestep of decimated LFP
ephysParams.decimateFactorPass1 = 6; %note: product of the two decimate factors should be 30, if 1 khz samples desired
ephysParams.decimateFactorPass2 = 5;
ephysParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
hp1Hz = designfilt('highpassiir', 'FilterOrder',8,'PassbandFrequency',1, ...
  'StopbandAttenuation',100,'PassbandRipple',0.5,'SampleRate',1000);     %#ok
% note: with these specifications, returns a 48th order butterworth filter
butter1Hz200Hz_v1 = designfilt('bandpassiir','DesignMethod','butter','PassbandFrequency1',1,'PassbandFrequency2',200,...
  'SampleRate',1000,'MatchExactly','passband','StopbandFrequency1',0.67,'StopbandFrequency2',250);
[tmp1,tmp2] = butter(4,[1/500,200/500]);
butter1Hz200Hz_v2 = [tmp1,tmp2];        %#ok
ephysParams.filter = butter1Hz200Hz_v1; % if filtering desired, ephysFilter is a digitalFilter 

% parameters preprocessAnalogIn, see function for details
analogInParams.needAccel = 0;
analogInParams.needEyes = 0;
analogInParams.accelChannels = []; %2d array, numAccelerometers x numChannelsPerAccelerometer
analogInParams.accelChannelNames = {}; % row vector, whose length matches num rows in accelChannels
analogInParams.eyeChannels = [];

% parameters preprocessLogFile, see function for details
stimSyncParams.usePhotodiode = 0;        %#ok

% parameters for excludeStimuli, see function for details
excludeStimParams.fixPre = 0; %ms
excludeStimParams.fixPost = 0; %ms
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

correlParams.maxShift = 50;
correlParams.matchTimeRanges = 1;
correlParams.timeDifferenceBound = [0,200];
correlParams.normalize = 1;
correlParams.useJacknife = 1;
correlParams.jacknifeDraws = 100;
switch machine
  case 'laptop'
    correlParams.jacknifeParallelWorkers = 0;    %#ok
  case 'hopper'
    correlParams.jacknifeParallelWorkers = 0;   %#ok
  case 'turing'
    corelParams.jacknifeParallelWorkers = 20;    %#ok
end

lfpAlignParams.samPerMS = 1; % because this is after decimation
lfpAlignParams.msPreAlign = psthParams.psthPre+tfParams.movingWin(1)/2; 
lfpAlignParams.msPostAlign = psthParams.psthImDur+psthParams.psthPost+tfParams.movingWin(1)/2;
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
plotSwitch.calcLatencyRF = 0;
plotSwitch.calcEvokedPowerRF = 0;
plotSwitch.evokedPsthMuaMultiCh = 0;
plotSwitch.evokedByCategory = 0;
plotSwitch.colorPsthEvoked = 0;
plotSwitch.linePsthEvoked = 0;
plotSwitch.runSummary = 0;
plotSwitch.runSummaryImMeanSub = 0;
plotSwitch.runSummaryImMeanSubDiv = 0;
plotSwitch.lfpPowerMuaScatter = 1; 
plotSwitch.lfpPeakToPeakMuaScatter = 1;
plotSwitch.lfpLatencyMuaLatency = 1;
plotSwitch.lfpPowerAcrossChannels = 1;
plotSwitch.lfpPeakToPeakAcrossChannels = 1;
plotSwitch.lfpLatencyShiftAcrossChannels = 1;
plotSwitch.singleTrialLfpByCategory = 1;
plotSwitch.lfpSpectraByCategory = 1;
plotSwitch.spikeSpectraByCategory = 1;
plotSwitch.SpikeSpectraTfByImage = 1;
plotSwitch.lfpSpectraTfByImage = 1;
plotSwitch.tfSpectraByCategory = 0;
plotSwitch.tfErrs = 1;           %#ok

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
calcSwitch.faceSelectIndex = 0;
calcSwitch.faceSelectIndexEarly = 0;
calcSwitch.faceSelectIndexLate = 0;
calcSwitch.inducedTrialMagnitudeCorrection = 0;
calcSwitch.evokedSpectra = 1;
calcSwitch.inducedSpectra = 0;
calcSwitch.evokedImageTF = 0;
calcSwitch.inducedImageTF = 0;
calcSwitch.evokedCatTF = 1;
calcSwitch.inducedCatTF = 0;
calcSwitch.meanEvokedTF = 1;
calcSwitch.trialMeanSpectra = 0;
calcSwitch.coherenceByCategory = 1;
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

