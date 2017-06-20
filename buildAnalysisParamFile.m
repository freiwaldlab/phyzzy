function [ analysisParamsFilename ] = buildAnalysisParamFile( )
%buildAnalysisParamFile saves a mat file of parameters, which control the
%behavior of analyzeSession
%   todo: option to load 'fixed' params from file, for ease accross days


%%%%%%%  USER PARAMETERS, EDIT ROUTINELY %%%%%%%%
runNum = '002';
dateSubject = '170405ALAN'; 
inRig = 0;
if inRig
  ephysVolume = '/Volumes/Users-1/User/Desktop';
  stimulusLogVolume = '/Volumes/Users/FreiwaldLab/Desktop';
else
  ephysVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Blackrock'; 
  stimulusLogVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Visiko';
end
outputVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Analyzed';
analysisLabel = 'Basic';
analysisParamsFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'
picParamsFilename = '/Users/stephenserene/Desktop/Freiwald/AnalysisSerene/StimParamsFullFOB3.mat';
categoryListSlim = {'humanFace','monkeyFace','place','fruit','humanBody','monkeyBody','techno'}; %minimal cat list for clean plots
saveFig = 1;
exportFig = 0;
saveFigData = 0;
verbosity = 'INFO'; %other options, 'DEBUG', 'VERBOSE';


% parameters preprocessSpikes and preprocessLFP, see functions for details
ephysParams.needLFP = 1;
ephysParams.needSpikes = 1;
ephysParams.spikeChannels = [1,35,33]; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
ephysParams.lfpChannels = [1,35,33]; 
ephysParams.channelNames = {'ML','AL','AM'};
ephysParams.lfpChannelScaleBy = [8191/32764 8191/32764 8191/32764]; %converts raw values to microvolts
ephysParams.common_ref = [0, 35, 35]; %not yet implemented; will allow software re-refrence across headstages
ephysParams.stimulationChannels = []; %not yet implemented; will read stimulation currents recorded at headstage
ephysParams.cPtCal = 1/30; % conversion from spike sample indices to timestep of decimated LFP
ephysParams.decimateFactorPass1 = 6; %note: product of the two decimate factors should be 30, if 1 khz samples desired
ephysParams.decimateFactorPass2 = 5;
ephysParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
hp1Hz = designfilt('highpassiir', 'FilterOrder',8,'PassbandFrequency',1, ...
  'StopbandAttenuation',100,'PassbandRipple',0.5,'SampleRate',1000); 
% note: with these specifications, returns a 48th order butterworth filter
butter1Hz200Hz_v1 = designfilt('bandpassiir','DesignMethod','butter','PassbandFrequency1',1,'PassbandFrequency2',200,...
  'SampleRate',1000,'MatchExactly','passband','StopbandFrequency1',0.67,'StopbandFrequency2',250);
[tmp1,tmp2] = butter(4,[1/500,200/500]);
butter1Hz200Hz_v2 = [tmp1,tmp2];
ephysParams.filter = butter1Hz200Hz_v1; % if filtering desired, ephysFilter is a digitalFilter 

% parameters preprocessAnalogIn, see function for details
analogInParams.needAccel = 0;
analogInParams.needEyes = 0;
analogInParams.accelChannels = []; %2d array, numAccelerometers x numChannelsPerAccelerometer
analogInParams.accelChannelNames = {}; % row vector, whose length matches num rows in accelChannels
analogInParams.eyeChannels = [];

% parameters preprocessLogFile, see function for details
stimSyncParams.usePhotodiode = 0;

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
psthParams.smoothingWidth = 10; %psth smoothing width, in ms
% TW=3 with T=.2, then W = 15 Hz (5 tapers)
% TW=1.5 with T=.1, then W = 15 Hz (2 tapers)
% TW = 1.5 with T=.2, then W = 7.5 Hz (2 tapers)
chr_params.tapers = [3 5]; %[3 5] is chronux default; 
chr_params.pad = 1;
chr_params.fs = 1;
chr_params.trialave = 1;
chr_params.err = [1 .05];
chr_params.fpass = [0 .1];
tfParams.movingWin = [200 5];
tfParams.specgramRowAve = 0;

lfpAlignParams.samPerMS = 1; % because this is after decimation
lfpAlignParams.msPreAlign = psthParams.psthPre+tfParams.movingWin(1)/2; 
lfpAlignParams.msPostAlign = psthParams.psthImDur+psthParams.psthPost+tfParams.movingWin(1)/2;
spikeAlignParams.preAlign = psthParams.psthPre+3*psthParams.smoothingWidth;
spikeAlignParams.postAlign = psthParams.psthImDur+psthParams.psthPost+3*psthParams.smoothingWidth;
% for lfps, constrain first and (optional) last [n m] samples to 0 mean
useDCSUB = 0;
if useDCSUB
  %lfpAlignParams.DCSUB_SAM = [lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10; 0, 0 ]; % 0-th order 
  lfpAlignParams.DCSUB_SAM = [lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10;lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10 ]; % 1st order 
else
  lfpAlignParams.DCSUB_SAM = 0;
end
%
frCalcOn = 60;
frCalcOff = 0; %note: if frCalcOff < frCalcOn, will update when psthImDur is set after reading log file 
frCalcOnEarly = 60;
frCalcOffEarly = 260;
frCalcOnLate = 260;
frCalcOffLate = 460;
%
psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'

% Boolean variables to specify which computations to perform; TODO: read
% from config file, eventually with conditional on log file info
makeImPSTH = 1;
makeCatPSTH = 1;

calcCoherenceRFcpt = 0;
calcCoherenceRFcc = 0;
calcCoherenceRFptpt = 0;
calcGrangerRF = 0;

plotSwitch.prefImRaster = 0;
plotSwitch.prefImRasterEvokedOverlay = 0;
plotSwitch.prefImMultiChRasterEvokedOverlay = 0;
plotSwitch.imageTuningSorted = 0;
plotSwitch.categoryPrefBarPlot = 0;
plotSwitch.categoryPrefBarPlotEarly = 0;
plotSwitch.categoryPrefBarPlotLate = 0;
plotSwitch.tuningCurves = 0;
plotSwitch.tuningCurvesEarly = 0;
plotSwitch.tuningCurvesLate = 0;
plotSwitch.calcLatencyRF = 0;
plotSwitch.calcEvokedPowerRF = 0;
plotSwitch.faceVnonEvokedPotential = 0;
plotSwitch.faceVnonEvokedMuaMultiCh = 0;
plotSwitch.evokedByCategory = 0;
plotSwitch.psthEvokedByCategory = 0;
plotSwitch.runSummary = 0;
plotSwitch.runSummaryImMeanSub = 0;
plotSwitch.runSummaryImMeanSubDiv = 0;
plotSwitch.lfpPowerMuaScatterAll = 0; 
plotSwitch.lfpPeakToPeakMuaScatterAll = 0;
plotSwitch.lfpPowerMuaScatterAllEarly = 0;
plotSwitch.lfpPeakToPeakMuaScatterAllEarly = 0;
plotSwitch.lfpPowerMuaScatterAllLate = 0;
plotSwitch.lfpPeakToPeakMuaScatterAllLate = 0;
plotSwitch.lfpLatencyMuaLatency = 0;
plotSwitch.lfpLatencyMuaLatencyEarly = 0;
plotSwitch.lfpLatencyMuaLatencyLate = 0;
plotSwitch.lfpPowerAcrossChannels = 0;
plotSwitch.lfpPeakToPeakAcrossChannels = 0;
plotSwitch.lfpLatencyShiftAcrossChannels = 0;
plotSwitch.singleTrialLfpByCategory = 1;
plotSwitch.lfpSpectraByCategory = 1;
plotSwitch.spikeSpectraByCategory = 1;
plotSwitch.SpikeSpectraTfByImage = 1;
plotSwitch.lfpSpectraTfByImage = 1;
plotSwitch.tfErrs = 1;

%%%%
analysisGroups.spectraByCategory.groups = {{'face','nonface'}};  %todo: add spectra diff?
analysisGroups.spectraByCategory.names = {{'faceVnon'}};
analysisGroups.spectraByCategory.colors = {{'faceVnon'}};
%
analysisGroups.tfSpectraByCategory.groups = {{'face'},{'nonface'},{'object'},{'body'}}; %todo: add tf spectra diff?
analysisGroups.tfSpectraByCategory.names = {{'face'},{'nonface'},{'object'},{'body'}};
%
analysisGroups.lfpSingleTrialsByCategory.groups = {{'face','nonface'}};
analysisGroups.lfpSingleTrialsByCategory.names = {{'faceVnon'}};
%
analysisGroups.coherenceByCategory.groups = {{'face','nonface'},{'face','object','body'},{'humanFace','monkeyFace','place','fruit','humanBody','monkeyBody','hand','techno'}};
analysisGroups.coherenceByCategory.colors = {{'r','b'},{'r','g','b'},{'b','c','y','g','m','r','k','k'}};
analysisGroups.coherenceByCategory.names = {'faceVnon','fob','slimCats'};
%
analysisGroups.tfCouplingByCategory.groups = {{'face'},{'nonface'},{'object'},{'body'}};

analysisGroups.byImage = {};
analysisGroupColors.byImage = {};
%%%%%

calcSwitch.faceSelectIndex = 0;
calcSwitch.faceSelectIndexEarly = 0;
calcSwitch.faceSelectIndexLate = 0;
calcSwitch.inducedTrialMagnitudeCorrection = 0;
calcSwitch.evokedSpectra = 1;
calcSwitch.inducedSpectra = 1;
calcSwitch.evokedImageTF = 0;
calcSwitch.inducedImageTF = 0;
calcSwitch.evokedCatTF = 0;
calcSwitch.inducedCatTF = 0;
calcSwitch.meanEvokedTF = 0;
calcSwitch.coherenceByCategory = 0;
calcSwitch.spikeTimes = 0;
calcSwitch.useJacknife = 0;

%%% set paths and directories, EDIT RARELY %%%
analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);
lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);
spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
taskFilename = sprintf('%s/%s/%s0%s.log',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
analysisParamsFilename = strcat(outDir,analysisParamsFilenameStem);
%
load('cocode2.mat');
psthColormap = map;
%
if ~exist(outDir,'dir')
  mkdir(outDir);
end
save(analysisParamsFilename);
end

