function [ analysisParamFilename ] = buildAnalysisParamFileSocialVids( )    
%buildAnalysisParamFile saves a mat file of parameters, which control the
%behavior of analyzeSession
%   todo: option to load 'fixed' params from file, for ease accross days

%%%%%%%  USER PARAMETERS, EDIT ROUTINELY %%%%%%%%
runNum = '003';
dateSubject = '20180828Mo';
[~, machine] = system('hostname');
machine = machine(~isspace(machine));

switch machine
  case 'Alienware_FA'
    ephysVolume = 'D:/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Data'; 
    stimulusLogVolume = 'D:/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Data';
    outputVolume = 'D:/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Analyzed';
    stimParamsFilename = 'D:/Onedrive/Lab/ESIN_Ephys_Files/Analysis/phyzzy/StimParamFileSocialVids.mat';   %#ok    
  case 'DESKTOP-MAT9KQ7'
    ephysVolume = 'C:/Users/aboha/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Data';
    stimulusLogVolume = 'C:/Users/aboha/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Data';
    outputVolume = 'C:/Users/aboha/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Analyzed';
    stimParamsFilename = 'C:/Users/aboha/Onedrive/Lab/ESIN_Ephys_Files/Analysis/phyzzy/StimParamFileSocialVids.mat';   %#ok
  case 'FA_Desktop_Home'
    ephysVolume = 'E:/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Data';
    stimulusLogVolume = 'E:/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Data';
    outputVolume = 'E:/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Analyzed';
    stimParamsFilename = 'E:/Onedrive/Lab/ESIN_Ephys_Files/Analysis/phyzzy/StimParamFileSocialVids.mat';   %#ok
end
analysisLabel = 'BasicSorted';
analysisParamFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'
preprocessedDataFilenameStem = 'preprocessedData.mat';
saveFig = 1;                %#ok
closeFig = 0;               %#ok
exportFig = 0;              %#ok
saveFigData = 0;            %#ok
savePreprocessed = 1;       %#ok
verbosity = 'INFO';         %other options, 'DEBUG', 'VERBOSE';


% parameters preprocessSpikes and preprocessLFP, see functions for details
ephysParams.needLFP = 1;
ephysParams.needSpikes = 1;
ephysParams.spikeChannels = [1]; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
ephysParams.lfpChannels = [1]; 
ephysParams.channelNames = {'8Bm'};
ephysParams.lfpChannelScaleBy = [8191/32764]; %converts raw values to microvolts
ephysParams.offlineSorted = 1; %Checks for a '*s.mat' Structure in the folder, with resorted spikes.
ephysParams.commonRef = [0]; %not yet implemented; will allow software re-refrence across headstages
ephysParams.stimulationChannels = []; %not yet implemented; will read stimulation currents recorded at headstage
ephysParams.cPtCal = 1/30; % conversion from spike sample indices to timestep of decimated LFP
ephysParams.decimateFactorPass1 = 6; %note: product of the two decimate factors should be 30, if 1 khz samples desired
ephysParams.decimateFactorPass2 = 5;
ephysParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)
%note: use Blackrock indexing for unitsToUnsort and unitsToDiscard, so unsorted is 0, first defined unit is 1, etc.
ephysParams.unitsToUnsort = {[]}; %these units will be re-grouped with u0
ephysParams.unitsToDiscard = {[]}; %these units will be considered noise and discarded
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
analogInParams.analogInChannels = [129,130]; 
analogInParams.channelNames = {'eyeX','eyeY'};
analogInParams.analogInChannelScaleBy = [5/32764 5/32764 5/32764]; %converts raw values to volts
analogInParams.decimateFactorPass1 = 1; 
analogInParams.decimateFactorPass2 = 1;
analogInParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to analogIn (should be raw rate/productOfDecimateFactors)
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
butter200Hz_v1 = designfilt('lowpassiir', 'PassbandFrequency', 120, 'StopbandFrequency', 480, 'PassbandRipple', 1,...
  'StopbandAttenuation', 60, 'SampleRate', 1000, 'MatchExactly', 'passband');  %this returns a 3rd order iir filter
analogInParams.filters = {0,0,0};%{butter200Hz_v1;butter200Hz_v1;butter200Hz_v1}; %filter channel i if filters{i} is digital filter or 1x2 numeric array
analogInParams.plotFilterResult = 1; %#ok

photodiodeParams.needStrobe = 1;
photodiodeParams.inputDataType = 'blackrockFilename';
photodiodeParams.dataChannel = 131;
photodiodeParams.dataLoader = []; %Incase you're using something besides a raw array or a blackrock file.
photodiodeParams.peakTimeOffset = 0; %this is the offset, in ms, of the peak from the event it's coupled to (note: will be subtracted, so should be > 0 if peak follows event, type: numeric)
photodiodeParams.strobeTroughs = 1; %Strobe causes troughs.
photodiodeParams.peakFreq = 85; %approximate number of peaks per second
photodiodeParams.levelCalibType = 'autoAndPlot';
%'hardcode', 'hardcodeAndPlot', 'hardcodeAndCheck', 'auto', 'autoAndPlot',
%'autoAndCheck', 'manual'
photodiodeParams.numLevels = 1;
photodiodeParams.saveCalibFile = 0;
photodiodeParams.centerCornerOffset = 5.3;
photodiodeParams.stimulusTriggerChannel = [];
photodiodeParams.peaksToPlot = 50; %number of peaks to show in calibration plots
photodiodeParams.cleanPeaks = 1;
photodiodeParams.useRisingEdge = 0; % 0 = peaks, 1 = rising edge.
photodiodeParams.numLevels = 3; %Number of levels the strobe is acting at.
photodiodeParams.checkHighLowAlternation = 0;
photodiodeParams.minPeakNumInLevel = 30;
%photodiodeParams.samplingFreq = 30000; %Defined in Blackrock
photodiodeParams.hardcodeFromFile = 0; %(only needed for hardcode calib types)
photodiodeParams.saveCalibFile = 0;
photodiodeParams.saveFigures = 0;
photodiodeParams.displayStats = 1;

lineNoiseTriggerParams.needStrobe = 0;


%     - cleanPeaks: if true, keep peaks only if the signal dropped below the low-peak threshold between them;
%                   for a series of peaks without a sub-threshold trough, keep only the highest (type: logical)
%     - useRisingEdge: define threshold using peaks, then use rising edge
%       threshold crossing for trigger times; optional, default 0 (type: logical)
%     for trigger times
%     - checkHighLowAlternation: (optional, default zero, only used if numLevels == 2)
%     - numLevels: (type: int)
%     - minPeakNumInLevel: min number of peaks required to define a level (type: int)
%     - samplingFreq: not needed if inputDataType is blackrockFilename, or
%                     a custom handle that sets sampling freq (numeric)
%     - peaksToPlot: number of peaks to show in calibration plots; suggested value > 10 (int)
%     - peakFreq: approximate number of peaks per second
%     - hardcodeFromFile: (type: logical) (only needed for hardcode calib types)
%     - saveCalibFile: (type: logical)
%     - inputCalibrationFile: filename ending in .mat, contains field lowThreshold (numeric)
%                             if numLevels >= 2, also contains field highThreshold (numeric)
%                             if numLevels == 3, also contains field midThreshold (numeric)
%                             (required only if hardcodeFromFile): (type: string)    
%     - outputCalibrationFile: filename ending in .mat (required only if saveCalibFile): (type:string)
%     - saveFigures: (type: logical)
%     - calibFigFname: figure filename, without .fig stem; required only if saveFigures (type:string)
%     - triggersFigFname: figure filename, without .fig stem; required only if saveFigures (type:string)
%     - dateSubject: MMDDYYNAME or similar; required only if saveFigures (type: string)
%     - runNum: eg 002; required only if saveFigures (type: string)
%     - outDir: directory for figure and calib. output; include the final '/' in the path (type: string) 
%               (required only if saveFigures or saveCalibFile is 1)
%

% parameters preprocessLogFile, see function for details
stimSyncParams.logProcessor = @preprocessLogFileMonkeyLogic;
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
%excludeStimParams.excludeTrials = @excludeTrialsMonkeyLogic;
excludeStimParams.needExcludeTrials = 0;
excludeStimParams.fixPre = 100; %ms
excludeStimParams.fixPost = 100; %ms
excludeStimParams.flashPre = 0;  %ms
excludeStimParams.flashPost = 0; %ms
excludeStimParams.juicePre = 0; % optional, ms
excludeStimParams.juicePost = 0; % optional, ms
excludeStimParams.DEBUG = 0; % makes exclusion criterion plots if true
% additional optional excludeStimParams: accel1, accel2, minStimDur (ms)

psthParams.psthPre = 800; % if e.g. +200, then start psth 200ms before trial onset; 
psthParams.psthImDur = 2800;  % only need to set this for variable length stim runs; else, comes from log file
psthParams.psthPost = 0;
psthParams.smoothingWidth = 50;  %psth smoothing width, in ms
psthParams.errorType = 1; %chronux convention: 1 is poisfStimson, 2 is trialwise bootstrap, 3 is across trial std for binned spikes, bootstrap for spike times 
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
% or function handle, which will receive the minimum stimulus duration in
% the run as an input.
frEpochsCell = {{60, @(stimDur) stimDur+60};...
                {60, 260}; ...
                {260, @(stimDur) stimDur+60}}; %#ok

plotSwitch.imagePsth = 1;
plotSwitch.categoryPsth = 1;
plotSwitch.prefImRaster = 0;
plotSwitch.prefImRasterEvokedOverlay = 0;
plotSwitch.prefImMultiChRasterEvokedOverlay = 0;
plotSwitch.imageTuningSorted = 1;
plotSwitch.stimPrefBarPlot = 1;
plotSwitch.stimPrefBarPlotEarly = 0;
plotSwitch.stimPrefBarPlotLate = 0;
plotSwitch.tuningCurves = 0;
plotSwitch.tuningCurvesEarly = 0;
plotSwitch.tuningCurvesLate = 0;
plotSwitch.RF = 0;
plotSwitch.rfEarly = 0;
plotSwitch.rfLate = 0;
plotSwitch.latencyRF = 0;
plotSwitch.evokedPowerRF = 0;
plotSwitch.evokedPsthMuaMultiCh = 0;
plotSwitch.evokedByCategory = 0;
plotSwitch.analogInByItem = 0;
plotSwitch.analogInDerivativesByItem = 0;
plotSwitch.colorPsthEvoked = 0;
plotSwitch.linePsthEvoked = 0;
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
%Defined for Groups of 2, A-B/A+B type index.
analysisGroups.selectivityIndex.groups = {{'socialInteraction';'nonInteraction'},{'twoMonkeys';'twoHumans'},{'socialInteraction';'goalDirected'}};

%Barplots showing average activity across all members of a catagory
analysisGroups.stimPrefBarPlot.groups = {{{'socialInteraction';'nonInteraction';'goalDirected';'idle';'scene';'scramble'}}};
analysisGroups.stimPrefBarPlot.colors  = {{{'b';'c';'y';'g';'m';'g'}}};
analysisGroups.stimPrefBarPlot.names = {'Barplots per label'};
analysisGroups.stimPrefBarPlot.groupDepth = 2; %2 subplots, each figure is defined by a cell array in the first item (groups).

%Relevant for things like Scatterplots, where you want a color for each
%stimulus.
analysisGroups.stimulusLabelGroups.groups = {{'socialInteraction';'nonInteraction';'goalDirected';'idle';'scene';'scramble'}};
analysisGroups.stimulusLabelGroups.names = {'Preferred Stimulus'};
analysisGroups.stimulusLabelGroups.colors = {{'b';'c';'y';'g';'m';'g'}};

%Essentially LFP selectivity/strength/quality
analysisGroups.evokedPotentials.groups = {{'socialInteraction';'nonInteraction';'objects'}};
analysisGroups.evokedPotentials.names = {'socVobj'};
analysisGroups.evokedPotentials.colors = {{'b';'g';'c'}};

%Looks like evokedpotentials, but pulls from a different analog channels
%(like pupul or eye).
analysisGroups.analogInPotentials.groups = {{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'}};
analysisGroups.analogInPotentials.channels = {[1; 2]};
analysisGroups.analogInPotentials.names = {'eyePositions,fobPlus'};
analysisGroups.analogInPotentials.units = {'degrees visual angle'};
analysisGroups.analogInPotentials.colors = {{'b';'c';'y';'g';'m';'r';'k'}};

%Same thing, but derivatives.
analysisGroups.analogInDerivatives.groups = {{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'techno'}};
analysisGroups.analogInDerivatives.channels = {[1; 2]};
analysisGroups.analogInDerivatives.names = {'eyeVelocity,fobPlus'};
analysisGroups.analogInDerivatives.units = {'degrees visual angle/sec'};
analysisGroups.analogInDerivatives.colors = {{'b';'c';'y';'g';'m';'r';'k'}};

%Makes subplots w/ PSTH on top and evoked potential on the bottom
analysisGroups.colorPsthEvoked.groups = {{'socialInteraction';'nonInteraction';'objects'}};
analysisGroups.colorPsthEvoked.names = {'socVobj'};
analysisGroups.colorPsthEvoked.colors = {{'b';'g';'c'}};

%same as above, but shows error bars, harder to see catagory selectivity
%though
analysisGroups.linePsthEvoked.groups = {{'socialInteraction';'nonInteraction';'objects'}};
analysisGroups.linePsthEvoked.names = {'socVobj'};
analysisGroups.linePsthEvoked.colors = {{'b';'g'}};

%Same as above, but everything ontop of eachother in 1 panel.
analysisGroups.evokedPsthOnePane.groups = {{'face';'nonface'}};
analysisGroups.evokedPsthOnePane.names = {'faceVnon'};

%Creates tuning curves for units, you should have some meaningful numeric
%descriptor.
analysisGroups.tuningCurves.groups = {{'humanFaceL90','humanFaceL45','humanFaceFront','humanFaceR45','humanFaceR90'},...
  {'monkeyFaceL90','monkeyFaceL45','monkeyFaceFront','monkeyFaceR45','monkeyFaceR90'}}; %can be images or categories
analysisGroups.tuningCurves.paramValues = {[-90 -45 0 45 90], [-90 -45 0 45 90]};
analysisGroups.tuningCurves.paramLabels = {'viewing angle (degrees)','viewing angle (degrees)'};
analysisGroups.tuningCurves.names = {'Human face view','Monkey face view'};

%Power spectra of SPIKE TIMES (i.e. 1's at spike times, flat elsewhere).
%Used for both Spikes and LFPs.
analysisGroups.spectraByCategory.groups = {{'socialInteraction';'nonInteraction';}};  %todo: add spectra diff?
analysisGroups.spectraByCategory.names = {'SocialVNonSocial'};
analysisGroups.spectraByCategory.colors = {{'r';'b'}};

%Calculates the same spectra, but w/ sliding windows. sometimes the Power
%spectrum changes overtime.
analysisGroups.tfSpectraByCategory.groups = {{'socialInteraction';'nonInteraction';}};%{'object'};{'body'}      %todo: add tf spectra diff?
analysisGroups.tfSpectraByCategory.names = {'socialInteraction';'interaction'};%'nonface';'object';'body'

%Evoked potential plot, shows individual traces for a bunch of trials.
analysisGroups.lfpSingleTrialsByCategory.groups = {{'socialInteraction';'nonInteraction';}};
analysisGroups.lfpSingleTrialsByCategory.names = {'SocialVNonSocial'};

%Coherence between LFP time series and spike time series w/i single
%channel. Does this on a trial by trial basis, and then averages across all
%members of each group. 
analysisGroups.coherenceByCategory.groups = {{'socialInteraction';'nonInteraction';}}; %{'face';'object';'body'};{'humanFace';'monkeyFace';'place';'fruit';'humanBody';'monkeyBody';'hand';'techno'}
analysisGroups.coherenceByCategory.colors = {{'r';'b'}}; %{'r';'g';'b'};{'b';'c';'y';'g';'m';'r';'k';'k'}
analysisGroups.coherenceByCategory.names = {'SocialVNonSocial'}; %'fob';'slimCats'

%Calculates the same as above but in sliding windows.
analysisGroups.tfCouplingByCategory.groups = {{'socialInteraction'};{'nonInteraction';};{'objects'};{'goalDirected'}}; %#ok
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
calcSwitch.spikeTimes = 1;
calcSwitch.useJacknife = 0;      

if calcSwitch.useJacknife
  chronuxParams.err(1) = 2; %#ok
end

%%% set paths and directories, EDIT RARELY %%%
analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);   %#ok
lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);        
spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
taskFilename = sprintf('%s/%s/%s%s.mat',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
photodiodeFilename = lfpFilename;                %#ok
lineNoiseTriggerFilename = lfpFilename; %#ok
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

