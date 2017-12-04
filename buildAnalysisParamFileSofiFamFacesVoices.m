function [ analysisParamsFilename ] = buildAnalysisParamFileSofiFamFacesVoices( )    
%buildAnalysisParamFile saves a mat file of parameters, which control the
%behavior of analyzeSession
%   todo: option to load 'fixed' params from file, for ease accross days


%%%%%%%  USER PARAMETERS, EDIT ROUTINELY %%%%%%%%
runNum = '004';
dateSubject = '171115Buster'; 
machine = 'turing';

switch machine
  case 'rig'
    ephysVolume = '';
    stimulusLogVolume = '/';
    outputVolume = '';
    stimParamsFilename = '/Volumes/FreiwaldShares/slandi/ephys/tests_2017/AnalysisSofi/StimParamsFamAudioVisual_v1.mat';   %#ok
  case 'laptop'
    ephysVolume = ''; 
    stimulusLogVolume = '';
    outputVolume = '';
    stimParamsFilename = '/Volumes/FreiwaldShares/slandi/ephys/tests_2017/AnalysisSofi/StimParamsFamAudioVisual_v1.mat';   %#ok
  case 'turing'
    ephysVolume = '/Volumes/FreiwaldShares/slandi/ephys/tests_2017/Blackrock'; 
    stimulusLogVolume = '/Volumes/FreiwaldShares/slandi/ephys/tests_2017/Behavior';
    outputVolume = '/Volumes/FreiwaldShares/slandi/ephys/tests_2017/Analyzed';
    stimParamsFilename = '/Volumes/FreiwaldShares/slandi/ephys/tests_2017/AnalysisSofi/StimParamsFamAudioVisual_v1.mat';   %#ok
end
analysisLabel = 'Basic';
analysisParamsFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'
preprocessedDataFilenameStem = 'preprocessedData.mat';
%categoryListSlim = {'humanFace','monkeyFace','place','fruit','humanBody','monkeyBody','techno'}; %minimal cat list for clean plots
saveFig = 1;           %#ok
exportFig = 0;         %#ok
saveFigData = 0;       %#ok
savePreprocessed = 1;  %#ok
verbosity = 'INFO'; %other options, 'DEBUG', 'VERBOSE';


% parameters preprocessSpikes and preprocessLFP, see functions for details
ephysParams.needLFP = 1;
ephysParams.needSpikes = 1;
ephysParams.spikeChannels = [1];%,35,33]; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
ephysParams.lfpChannels = [1];%;,35,33]; 
ephysParams.channelNames = {'TP'};%,'AL','AM'};
ephysParams.lfpChannelScaleBy = [8191/32764];% 8191/32764 8191/32764]; %converts raw values to microvolts
ephysParams.common_ref = [0];%, 35, 35]; %not yet implemented; will allow software re-refrence across headstages
ephysParams.stimulationChannels = []; %not yet implemented; will read stimulation currents recorded at headstage
ephysParams.cPtCal = 1/30; % conversion from spike sample indices to timestep of decimated LFP
ephysParams.decimateFactorPass1 = 6; %note: product of the two decimate factors should be 30, if 1 khz samples desired
ephysParams.decimateFactorPass2 = 5;
ephysParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)
%note: use Blackrock indexing for unitsToUnsort and unitsToDiscard, so unsorted is 0, first defined unit is 1, etc.
ephysParams.unitsToUnsort = {[]};%,[],[]}; %these units will be re-grouped with u0
ephysParams.unitsToDiscard = {[]};%,[],[]}; %these units will be considered noise and discarded
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
analogInParams.needAnalogIn = 1;
analogInParams.analogInChannels = [129,130,131];  % caspar's cff 
analogInParams.channelNames = {'eyeX','eyeY','eyeD'};
analogInParams.analogInChannelScaleBy = [5/32764 5/32764 5/32764]; %converts raw values to volts
analogInParams.decimateFactorPass1 = 1; %note: product of the two decimate factors should be 30, if 1 khz samples desired
analogInParams.decimateFactorPass2 = 1;
analogInParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to analogIn (should be raw rate/productOfDecimateFactors)
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
butter200Hz_v1 = designfilt('lowpassiir', 'PassbandFrequency', 120, 'StopbandFrequency', 480, 'PassbandRipple', 1,...
  'StopbandAttenuation', 60, 'SampleRate', 1000, 'MatchExactly', 'passband');  %this returns a 3rd order iir filter
analogInParams.filters = {0,0,0};%{butter200Hz_v1;butter200Hz_v1;butter200Hz_v1}; %filter channel i if filters{i} is digital filter or 1x2 numeric array
analogInParams.plotFilteredSignal = 1; %#ok

photodiodeParams.needPhotodiode = 0;
photodiodeParams.channels = [137]; %#ok
photodiodeParams.frameTriggerChannel = 1;
% parameters preprocessLogFile, see function for details
stimSyncParams.usePhotodiode = 0;        %#ok

%
audioParams.needAudio = 1;
audioParams.channels = 132;
stimSyncParams.useAudio = 1;        %#ok

eyeCalParams.needEyeCal = 0;
eyeCalParams.method = 'zeroEachFixation'; %hardcodeZero %autoZeroSingle
eyeCalParams.makePlots = 1;
eyeCalParams.eyeXChannelInd = 1;
eyeCalParams.eyeYChannelInd = 2;
eyeCalParams.eyeDChannelInd = 3;
eyeCalParams.gainX = 312;
eyeCalParams.gainY = 164;
eyeCalParams.flipX = 1;
eyeCalParams.flipY = 1; 
eyeCalParams.offsetX = -6.9;
eyeCalParams.offsetY = -6.8; 
eyeCalParams.minFixZeroTime = 1000; %ok
eyeCalParams.degrees = 1;

accelParams.needAccelCal = 0;
accelParams.accelChannels = {[4;5;6]};
accelParams.channelGains = {[1/.666 1/.666 1/.666]};
accelParams.calMethods = {'hardcode'}; %other option is 'calFile'; calibration method
accelParams.calFiles = {''}; %if method is 'calFile', an ns2 filename

% parameters for excludeStimuli, see function for details
excludeStimParams.fixPre = 0; %ms
excludeStimParams.fixPost = 0; %ms
excludeStimParams.flashPre = 0;  %ms
excludeStimParams.flashPost = 0; %ms
excludeStimParams.juicePre = 0; % optional, ms
excludeStimParams.juicePost = 0; % optional, ms
excludeStimParams.DEBUG = 0; % makes exclusion criterion plots if true
excludeStimParams.excludeFixation = 0; % do not exclude any trials due to bad fixation
% additional optional excludeStimParams: accel1, accel2, minStimDur (ms)

psthParams.psthPre = 100; % if e.g. +200, then start psth 200ms before trial onset; 
psthParams.psthImDur = 0;  % only need to set this for variable length stim runs; else, comes from log file
psthParams.psthPost = 500;
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
frEpochsCell = {{30, @(stimDur) stimDur+60};...
                {30, 300}; ...
                {300, 700}}; %#ok


% Boolean variables to specify which computations to perform; TODO: read
% from config file, eventually with conditional on log file info
calcCoherenceRFcpt = 0;  %#ok
calcCoherenceRFcc = 0;   %#ok
calcCoherenceRFptpt = 0; %#ok
calcGrangerRF = 0;       %#ok 

plotSwitch.imagePsth = 1;
plotSwitch.categoryPsth = 1;
plotSwitch.prefImRaster = 1;
plotSwitch.prefImRasterEvokedOverlay = 1;
plotSwitch.prefImMultiChRasterEvokedOverlay = 0;
plotSwitch.imageTuningSorted = 0;
plotSwitch.stimPrefBarPlot = 1;
plotSwitch.stimPrefBarPlotEarly = 1;
plotSwitch.stimPrefBarPlotLate = 1;
plotSwitch.tuningCurves = 0;
plotSwitch.tuningCurvesEarly = 0;
plotSwitch.tuningCurvesLate = 0;
plotSwitch.calcLatencyRF = 0;
plotSwitch.calcEvokedPowerRF = 0;
plotSwitch.evokedPsthMuaMultiCh = 0;
plotSwitch.evokedByCategory = 1;
plotSwitch.analogInByItem = 0;
plotSwitch.analogInDerivativesByItem = 0;
plotSwitch.colorPsthEvoked = 1;
plotSwitch.linePsthEvoked = 1;
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
analysisGroups.selectivityIndex.groups = {{'coo';'grunt';'faceVoice';'voice';'sound';'face';'grayBG'}};
%
analysisGroups.stimPrefBarPlot.groups = {{{'coo';'grunt';'faceVoice';'voice';'sound';'face';'grayBG'}}};                        
analysisGroups.stimPrefBarPlot.colors  = {{{'r';'c';'m';'g';'y';'b';'k'}}};
analysisGroups.stimPrefBarPlot.names = {'fobPlus'};
analysisGroups.stimPrefBarPlot.groupDepth = 2;

analysisGroups.stimulusLabelGroups.groups = {{'coo';'grunt';'faceVoice';'voice';'sound';'face';'grayBG'}};
analysisGroups.stimulusLabelGroups.names = {'fobPlus'};
analysisGroups.stimulusLabelGroups.colors = {{'r';'c';'m';'g';'y';'b';'k'}};
%
analysisGroups.evokedPotentials.groups = analysisGroups.stimulusLabelGroups.groups;
analysisGroups.evokedPotentials.names = {'fobPlus'};
analysisGroups.evokedPotentials.colors = analysisGroups.stimulusLabelGroups.colors;
%
analysisGroups.analogInPotentials.groups = analysisGroups.stimulusLabelGroups.groups;
analysisGroups.analogInPotentials.channels = {[1; 2]};
analysisGroups.analogInPotentials.names = {'eyePositions,fobPlus'};
analysisGroups.analogInPotentials.units = {'degrees visual angle'};
analysisGroups.analogInPotentials.colors = analysisGroups.stimulusLabelGroups.colors;
%
analysisGroups.analogInDerivatives.groups = analysisGroups.stimulusLabelGroups.groups;
analysisGroups.analogInDerivatives.channels = {[1; 2]};
analysisGroups.analogInDerivatives.names = {'eyeVelocity,fobPlus'};
analysisGroups.analogInDerivatives.units = {'degrees visual angle/sec'};
analysisGroups.analogInDerivatives.colors = analysisGroups.stimulusLabelGroups.colors;
%
analysisGroups.colorPsthEvoked.groups = analysisGroups.stimulusLabelGroups.groups;
analysisGroups.colorPsthEvoked.names = {'fobPlus'};
analysisGroups.colorPsthEvoked.colors = analysisGroups.stimulusLabelGroups.colors;
%
analysisGroups.linePsthEvoked.groups = analysisGroups.stimulusLabelGroups.groups;
analysisGroups.linePsthEvoked.names = {'fobPlus'};
analysisGroups.linePsthEvoked.colors = analysisGroups.stimulusLabelGroups.colors;
%
analysisGroups.evokedPsthOnePane.groups = analysisGroups.stimulusLabelGroups.groups;
analysisGroups.evokedPsthOnePane.names = {'fobPlus'};
%
analysisGroups.tuningCurves.groups = {{'Familiar','nonFamiliar'},...
  {'monkeyFace','humanFace'}}; %can be images or categories
analysisGroups.tuningCurves.paramValues = {[1 10], [1 10]};
analysisGroups.tuningCurves.paramLabels = {['fam' 'nonfam'],['monkey' 'human']};
analysisGroups.tuningCurves.names = {'Fam vs nonFam','Monkey vs Human'};
%
analysisGroups.spectraByCategory.groups = analysisGroups.stimulusLabelGroups.groups;  %todo: add spectra diff?
analysisGroups.spectraByCategory.names = {'fobPlus'};
analysisGroups.spectraByCategory.colors = analysisGroups.stimulusLabelGroups.colors;
%
analysisGroups.tfSpectraByCategory.groups = analysisGroups.stimulusLabelGroups.groups;  %{'object'};{'body'}      %todo: add tf spectra diff?
analysisGroups.tfSpectraByCategory.names = analysisGroups.stimulusLabelGroups.groups;  %'nonface';'object';'body'
%
analysisGroups.lfpSingleTrialsByCategory.groups = analysisGroups.stimulusLabelGroups.groups;  
analysisGroups.lfpSingleTrialsByCategory.names = analysisGroups.stimulusLabelGroups.groups;
%
analysisGroups.coherenceByCategory.groups = analysisGroups.stimulusLabelGroups.groups;  %{'face';'object';'body'};{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'hand';'techno'}
analysisGroups.coherenceByCategory.colors = analysisGroups.stimulusLabelGroups.colors; %{'r';'g';'b'};{'b';'c';'y';'g';'m';'r';'k';'k'}
analysisGroups.coherenceByCategory.names = {'fobPlus'}; %'fob';'slimCats'
%
analysisGroups.tfCouplingByCategory.groups = analysisGroups.stimulusLabelGroups.groups;

%analysisGroups.byImage = {};      %#ok
%analysisGroupColors.byImage = {}; %#ok
%%%%%

calcSwitch.categoryPSTH = 1;
calcSwitch.imagePSTH = 1;
calcSwitch.faceSelectIndex = 0;
calcSwitch.faceSelectIndexEarly = 0;
calcSwitch.faceSelectIndexLate = 0;
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
analogInFilename = sprintf('%s/%s/%s_%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);   %#ok
lfpFilename = sprintf('%s/%s/%s_%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);        %#ok
spikeFilename = sprintf('%s/%s/%s_%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
taskFilename = sprintf('%s/%s/%s0%s.log',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
photodiodeFilename = lfpFilename;                                                           %#ok

outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
analysisParamsFilename = strcat(outDir,analysisParamsFilenameStem);
preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok
%
load('cocode2.mat');
%load(psthParams.psthColormapFilename);
psthColormap = map;  %#ok
%
if ~exist(outDir,'dir')
  mkdir(outDir);
end
save(analysisParamsFilename);
end
