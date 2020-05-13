function [analysisParamFilename] = buildAnalysisParamFileSocialVids( varargin )    
%buildAnalysisParamFile saves a mat file of parameters, which control the
%behavior of processRun, runAnalysis

% %%%%%%%  USER PARAMETERS, EDIT ROUTINELY %%%%%%%%
runNum = '003';
dateSubject = '20180628Mo';
assert(~isempty(str2double(runNum)), 'Run number had letters, likely not normal run') %Just for batch runs where unique runs follow unconventional naming scheme.

[~, machine] = system('hostname');
machine = machine(~isspace(machine));

switch machine
  case 'homeDesktopWork'
    ephysVolume = slashSwap('C:\OneDrive\Lab\ESIN_Ephys_Files\Data');
    stimulusLogVolume = ephysVolume;
    outputVolume = 'C:/Analyzed';
    stimParamsFilename = slashSwap('C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
    stimDir = slashSwap('C:\OneDrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories');
    stimParamsFilename = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy\stimParamFileLib\StimParamFileSocialVids_Full.mat';   %#ok
  case 'Alienware_FA'
    ephysVolume = slashSwap('D:\Onedrive\Lab\ESIN_Ephys_Files\Data');
    stimulusLogVolume = ephysVolume;
    outputVolume = slashSwap('D:\DataAnalysis\March2020');
    stimParamsFilename = slashSwap('D:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
    stimDir = slashSwap('D:\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code');
  case {'turing.rockefeller.edu','hopper.rockefeller.edu'}
    ephysVolume = '/Freiwald/lab_files/raw_data/EPHYS/Farid_ESINRec/Data2018';
    stimulusLogVolume = ephysVolume;
    outputVolume = '/Freiwald/faboharb/analysis/Analyzed';
  case 'DataAnalysisPC'
    ephysVolume = slashSwap('\\BlackrockPC\nsp\Data');
    stimulusLogVolume = slashSwap('\\CONTROLLERPC\Monkeylogic Experiments');
    outputVolume = slashSwap('C:\Users\Farid\OneDrive\Lab\ESIN_Ephys_Files\Analysis\Analyzed_Rig');
    stimParamsFilename = slashSwap('C:\Users\Farid\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy\stimParamFileLib\StimParamFileSocialVids_Full.mat');                  %#ok
    stimDir = slashSwap('C:\Users\Farid\OneDrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories');
  case 'DESKTOP-MAT9KQ7'
    ephysVolume = 'C:/Users/aboha/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Data';
    stimulusLogVolume = ephysVolume;
    outputVolume = 'C:/Users/aboha/Onedrive/Lab/ESIN_Ephys_Files/Analysis/Analyzed';
    stimParamsFilename = 'C:/Users/aboha/Onedrive/Lab/ESIN_Ephys_Files/Analysis/phyzzy/StimParamFileSocialVids_V2.mat';   %#ok
  case 'SurfaceBook2FA'
    ephysVolume = 'C:/Users/Farid Aboharb/OneDrive/Lab/ESIN_Ephys_Files/Analysis/Data 2018';
    stimulusLogVolume = ephysVolume;
    outputVolume = 'C:/Users/Farid Aboharb/OneDrive/Lab/ESIN_Ephys_Files/Analysis/Analyzed';
    stimParamsFilename = 'C:/Users/Farid Aboharb/OneDrive/Lab/ESIN_Ephys_Files/Analysis/phyzzy/StimParamFileSocialVids_Full.mat';   %#ok
    stimDir = 'G:/StimuliForFaridfromJulia/SocialCategories';       
end

analysisLabel = 'Basic';
preprocessedDataFilenameStem = 'preprocessedData.mat';
analysisParamFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'

figStruct = struct();
figStruct.saveFig = 1;                
figStruct.closeFig = 0;               %#ok
figStruct.exportFig = 0;              %#ok
figStruct.saveFigData = 0;            %#ok
figStruct.figPos = [0 0 .6 0.7];      % Normalized units for figure position

savePreprocessed = 1;       %#ok
verbosity = 'INFO';         %other options, 'DEBUG', 'VERBOSE';

%% Plot switches
plotSwitch.pupilDilation = 1;               % plots image which illustrates continuous values for pupil dilation. 
plotSwitch.subEventAnalysis = 0;            % plot traces comparing activity surrounding an event (defined in eventData, generated w/ eventDetectionApp), vs null.
plotSwitch.imageEyeMap = 0;                 
plotSwitch.eyeCorralogram = 0;              % Eye Gram
plotSwitch.eyeStatsAnalysis = 1;            % use ClusterFix to generate a vector characterizing eye movements.
plotSwitch.attendedObject = 1;              % Vectors to distinguish where subject is looking. Required for prefImRasterColorCoded.
plotSwitch.eyeStimOverlay = 0;              % Visualize eye traces on stimuli.
plotSwitch.spikePupilCorr = 1;              % see the correlation between single trial PSTHes and pupil size.

plotSwitch.clusterOnEyePaths = 0;           % Resort spikes based on distinct eye paths, make "New events".
plotSwitch.stimPSTHoverlay = 0;             % grabs stimuli and plots the PSTH underneath.
plotSwitch.imagePsth = 1;
plotSwitch.categoryPsth = 0;
plotSwitch.stimCatANOVA = 1;
plotSwitch.prefImRaster = 0;                % Raster, Not color coded.
plotSwitch.prefImRasterColorCoded = 4;      % Raster, uses info from attendedObj switch. 1 is colored spikes, 2 is colored background, 3 is Saccade Image, 4 is pupil img.
plotSwitch.topStimToPlot = 5;               % Determines number of stimuli for which rasters are plotted.
plotSwitch.prefImRasterEvokedOverlay = 0;   % Produces images for MUA and Unsorted even if the same. Relies on sometihng in CatPSTH.
plotSwitch.prefImRasterAverageEvokedOverlay = 1;
plotSwitch.prefImMultiChRasterEvokedOverlay = 0;
plotSwitch.imageTuningSorted = 1;           % Barplot per image, Required for stimPSTHoverlay, sigStruct
plotSwitch.stimPrefBarPlot = 0;             % Per event bar graph.
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
plotSwitch.lfpSpectraByCategory = 1; % LFP Comparison 
plotSwitch.lfpAutocorrelTfByItem = 0;
plotSwitch.lfpAutocorrelByItem = 0;
plotSwitch.spikeSpectraByCategory = 0;
plotSwitch.spikeAutocorrelTfByItem = 0;
plotSwitch.spikeAutocorrelByItem = 0;
plotSwitch.spikeSpectraTfByImage = 0;
plotSwitch.lfpSpectraTfByImage = 0;
plotSwitch.couplingPhasesUnwrapped = 0;
plotSwitch.couplingPhasesAsOffsets = 0;
plotSwitch.couplingPhasesPolar = 0;
plotSwitch.tfSpectraByCategory = 1; %Do I want this?
plotSwitch.tfErrs = 0;           %#ok


%% Calc switches

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
calcSwitch.evokedCatSpikeTF = 0; %Required for one of the above plot switches to actually produce the figure, but crashes @ "spikesByItemBinned = spikesByCategoryBinned;" in the 2k lines.
calcSwitch.inducedCatSpikeTF = 0;
calcSwitch.evokedCatLFPTF = 0; %Required for one of the above plot switches to actually produce the figure, but crashes @ "spikesByItemBinned = spikesByCategoryBinned;" in the 2k lines.
calcSwitch.inducedCatLFPTF = 0;
calcSwitch.evokedCoupling = 0;
calcSwitch.inducedCoupling = 0;
calcSwitch.meanEvokedTF = 0;
calcSwitch.trialMeanSpectra = 0;
calcSwitch.coherenceByCategory = 0;
calcSwitch.spikeTimes = 0;
calcSwitch.useJacknife = 0;      

%% Parameters
% parameters preprocessSpikes and preprocessLFP, see functions for details
ephysParams.needLFP = 1;
ephysParams.needSpikes = 1;
autoChannelDetect = 1;
ephysParams.spikeChannels = [1]; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
ephysParams.lfpChannels = [1];
ephysParams.channelNames = {'Ch1'};
% ephysParams.channelNames = {'8B'};
ephysParams.lfpChannelScaleBy = [8191/32764, 8191/32764]; %converts raw values to microvolts
ephysParams.offlineSorted = 0;        % Checks for a '*.xls' Structure in the folder, with resorted spikes.
ephysParams.waveClus = 0;             % Does automated offline sorting using wave_clus.
ephysParams.paramHandle = @set_parameters; %Function which produces param struct for wave_clus. in wave_clus folder.
ephysParams.waveClusReclass = 1;      % Reclassify clusters (as defined by mean waveform proximity to threshold) to MUA.
ephysParams.waveClusMUAThreshold = 1.25; %scaling for reclassification of clusters as MUA. 1 = 0 scaling = no reclassification of clusters.
ephysParams.waveClusProjPlot = 1;     % Plots all the clusters in higher dimensional space (defined in type and number in wave_clus parameters).
ephysParams.waveClusClear = 1;        % 1 deletes all files associated with analysis (leaves processed NSX files), 2 deletes the entire associated '_parsed' folder.
ephysParams.commonRef = [0];          % not yet implemented; will allow software re-refrence across headstages
ephysParams.stimulationChannels = []; % not yet implemented; will read stimulation currents recorded at headstage
ephysParams.cPtCal = 1/30;            % conversion from spike sample indices to timestep of decimated LFP
ephysParams.decimateFactorPass1 = 6;  % note: product of the two decimate factors should be 30, if 1 khz samples desired
ephysParams.decimateFactorPass2 = 5;
ephysParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)
%note: use Blackrock indexing for unitsToUnsort and unitsToDiscard, so unsorted is 0, first defined unit is 1, etc.
ephysParams.unitsToUnsort = {[], []}; %these units will be re-grouped with u0
ephysParams.unitsToDiscard = {[], []}; %these units will be considered noise and discarded
ephysParams.spikeWaveformPca = 0;
ephysParams.plotSpikeWaveforms = 2; %0, 1 to build then close, 2 to build and leave open
ephysParams.spikeWaveformsColors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483]];
ephysParams.shiftSpikeWaveforms = 0;
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
hp1Hz = designfilt('highpassiir', 'FilterOrder',8,'PassbandFrequency',1, ...
  'StopbandAttenuation',100,'PassbandRipple',0.5,'SampleRate',1000);     %#ok
% note: with these specifications, returns a 48th order butterworth filter
butter1Hz200Hz_v1 = designfilt('bandpassiir','DesignMethod','butter','PassbandFrequency1',1,'PassbandFrequency2',200,...
  'SampleRate',1000,'MatchExactly','passband','StopbandFrequency1',0.67,'StopbandFrequency2',250);
[tmp1,tmp2] = butter(4,[1/500,200/500]);
butter1Hz200Hz_v2 = [tmp1,tmp2];        %#ok
clear('tmp1','tmp2')
ephysParams.filter = butter1Hz200Hz_v1; % if filtering desired, ephysFilter is a digitalFilter
ephysParams.plotFilterResult = 0; 
ephysParams.outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
ephysParams.saveFig = figStruct.saveFig;

% parameters preprocessAnalogIn, see function for details
analogInParams.needAnalogIn = 1;
analogInParams.analogInChannels = [129,130,137]; 
analogInParams.channelNames = {'eyeX','eyeY','pupil'};
analogInParams.channelUnits = {'dva','dva','au'};
analogInParams.analogInChannelScaleBy = [5/32764 5/32764 5/32764]; %converts raw values to volts
analogInParams.decimateFactorPass1 = 1; 
analogInParams.decimateFactorPass2 = 1;
analogInParams.filterPad = 0;
analogInParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to analogIn (should be raw rate/productOfDecimateFactors)
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
butter200Hz_v1 = designfilt('lowpassiir', 'PassbandFrequency', 120, 'StopbandFrequency', 480, 'PassbandRipple', 1,...
  'StopbandAttenuation', 60, 'SampleRate', 1000, 'MatchExactly', 'passband');  %this returns a 3rd order iir filter
nyqfrq = 1000 ./ 2;
%lowPass30Hz = fir2(60,[0,30./nyqfrq,30./nyqfrq,1],[1,1,0,0]); %60th order, 30 Hz low pass filter
lowPass30Hz = designfilt('lowpassfir', 'FilterOrder', 3, 'PassbandFrequency', 30, 'StopbandFrequency', 31, 'SampleRate', 1000);
analogInParams.filters = {0, 0, 0}; %{lowPass30Hz, lowPass30Hz, 0};%{butter200Hz_v1;butter200Hz_v1;butter200Hz_v1}; %filter channel i if filters{i} is digital filter or 1x2 numeric array
analogInParams.plotFilterResult = 1;

photodiodeParams.needStrobe = 1;
photodiodeParams.inputDataType = 'blackrockFilename';
photodiodeParams.dataChannel = 131;
photodiodeParams.dataLoader = [];         % Incase you're using something besides a raw array or a blackrock file.
photodiodeParams.peakTimeOffset = 0;      % this is the offset, in ms, of the peak from the event it's coupled to (note: will be subtracted, so should be > 0 if peak follows event, type: numeric)
photodiodeParams.strobeTroughs = 1;       % Strobe causes troughs.
photodiodeParams.peakFreq = 85;           % approximate number of peaks per second
photodiodeParams.levelCalibType = 'auto'; % 'hardcode', 'hardcodeAndPlot', 'hardcodeAndCheck', 'auto', 'autoAndPlot','autoAndCheck', 'manual'
photodiodeParams.numLevels = 1;
photodiodeParams.saveCalibFile = 0;
photodiodeParams.centerCornerOffset = 5.3;
photodiodeParams.stimulusTriggerChannel = [];
photodiodeParams.peaksToPlot = 50;        % number of peaks to show in calibration plots
photodiodeParams.cleanPeaks = 1;
photodiodeParams.useRisingEdge = 0;       % 0 = peaks, 1 = rising edge.
photodiodeParams.numLevels = 3;           % Number of levels the strobe is acting at.
photodiodeParams.levelHigh = 5350;
photodiodeParams.levelMid = 4000;
photodiodeParams.levelLow = 1250;
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
stimSyncParams.usePhotodiode = 1;
stimSyncParams.stimDir = stimDir;
stimSyncParams.outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum); %#ok
%
eyeCalParams.needEyeCal = 1;
eyeCalParams.method = 'monkeyLogic'; %'zeroEachFixation', 'monkeyLogic'
eyeCalParams.monkeyLogicShift = 1;
eyeCalParams.makePlots = 1;
eyeCalParams.eyeXChannelInd = 1;
eyeCalParams.eyeYChannelInd = 2;
eyeCalParams.eyeDChannelInd = 3;
eyeCalParams.samplingRate = 120; %Sampling rate in Hz of equipment
eyeCalParams.gainX = 0;
eyeCalParams.gainY = 0;
eyeCalParams.flipX = 0;
eyeCalParams.flipY = 0; 
eyeCalParams.offsetX = 0;
eyeCalParams.offsetY = 0; 
eyeCalParams.calFile = ''; %note: needed only when method = fromFile
eyeCalParams.fixOutLag = 10; 
eyeCalParams.minFixZeroTime = 1000; %#ok
eyeCalParams.analogInParams = analogInParams;

accelParams.needAccelCal = 0;
accelParams.accelChannels = {[4;5;6]};
accelParams.channelGains = {[1/.666 1/.666 1/.666]};
accelParams.calMethods = {'hardcode'}; %other option is 'calFile'; calibration method
accelParams.calFiles = {''}; %if method is 'calFile', an ns2 filename

% parameters for excludeStimuli, see function for details
%excludeStimParams.excludeTrials = @excludeTrialsMonkeyLogic;
excludeStimParams.excludeProcessor = @excludeTrialsMonkeyLogic;
excludeStimParams.needExcludeTrials = 1;
excludeStimParams.excludeFailed = 1;
excludeStimParams.excludeAfterFailed = 0;
excludeStimParams.frameDropThreshold = 5;
excludeStimParams.fixPre = 100; %ms
excludeStimParams.fixPost = 100; %ms
excludeStimParams.flashPre = 0;  %ms
excludeStimParams.flashPost = 0; %ms
excludeStimParams.juicePre = 0; % optional, ms
excludeStimParams.juicePost = 0; % optional, ms
excludeStimParams.DEBUG = 0; % makes exclusion criterion plots if true
% additional optional excludeStimParams: accel1, accel2, minStimDur (ms)

% TW=3 with T=.2, then W = 15 Hz (5 tapers)
% TW=1.5 with T=.1, then W = 15 Hz (2 tapers)
% TW = 1.5 with T=.2, then W = 7.5 Hz (2 tapers)
chronuxParams.tapers = [3 5]; %[3 5] is chronux default; 
chronuxParams.pad = 1;
chronuxParams.fs = 1;
chronuxParams.trialave = 1;
chronuxParams.err = [1 .05];  %note: first entry will be automatically switched to 2 if calcSwitch.useJacknife == 1
chronuxParams.fpass = [0 .1]; 
tfParams.movingWin = [300 5]; 
tfParams.specgramRowAve = 0;

psthParams.type = 'normal'; %options are 'normal', 'baselineSub', 'meanWhite'
psthParams.psthPre = 800; % if e.g. +200, then start psth 200ms before trial onset; 
psthParams.psthImDur = 2800;  % only need to set this for variable length stim runs; else, comes from log file
psthParams.psthPost = 500;
psthParams.latency = 0;
psthParams.movingWin = tfParams.movingWin;
psthParams.smoothingWidth = 40;  %psth smoothing width, in ms
psthParams.Zscore = 0;  % 0: raw PSTH, 1: pre-trial baseline subtraction Z Scored PSTH, 2: whole trial baseline subtracted Z Scored PSTH
psthParams.errorType = 2; %chronux convention: 1 is poisfStimson, 2 is trialwise bootstrap, 3 is across trial std for binned spikes, bootstrap for spike times 
psthParams.errorRangeZ = 1; %how many standard errors to show
psthParams.bootstrapSamples = 100;
psthParams.sortStim = 1;
psthParams.sortOrder = {'socialInteraction';'goalDirected';'idle';'objects';'scene';'scramble';'headTurn';'subEvents';'allStim'};
psthParams.psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'
load(psthParams.psthColormapFilename);
psthParams.colormap = map;

eyeStatsParams.psthPre = psthParams.psthPre;
eyeStatsParams.psthImDur = psthParams.psthImDur;
eyeStatsParams.stimDir = stimDir;
eyeStatsParams.lfpPaddedBy = tfParams.movingWin(1)/2;
eyeStatsParams.outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
eyeStatsParams.clusterFixLPFilterIn = 25;                 % Filter used for saccade detection internally in ClusterFix.

genStatsParams.ANOVAParams.target = 'socialInteraction';    % When performing a one way ANOVA, the label from groups which is used. the rest are 'non-' label.

subEventAnalysisParams.preAlign = 300;
subEventAnalysisParams.postAlign = 300;
subEventAnalysisParams.nullAllStim = 1;
subEventAnalysisParams.nullSampleMult = 10;       % For every n stimuli with the event, sample all other stimuli this number of times for the null distribution.
subEventAnalysisParams.psthParams = psthParams;
subEventAnalysisParams.psthParams.psthPre = subEventAnalysisParams.preAlign - 100;
subEventAnalysisParams.psthParams.psthImDur = subEventAnalysisParams.postAlign - 100;
subEventAnalysisParams.psthParams.psthPost = 0;
subEventAnalysisParams.psthParams.smoothingWidth = 10;
subEventAnalysisParams.testPeriod = [0 200];                            % The period during which spikes are counted and a t-test is performed.
subEventAnalysisParams.psthParams.movingWin = psthParams.movingWin;
subEventAnalysisParams.stimPlotParams.psthPre = psthParams.psthPre;
subEventAnalysisParams.stimPlotParams.psthPre = psthParams.psthPre;
subEventAnalysisParams.stimPlotParams.psthImDur = psthParams.psthImDur;
subEventAnalysisParams.stimPlotParams.psthPost = psthParams.psthPost;
subEventAnalysisParams.saveFig = 1;                
subEventAnalysisParams.closeFig = 0;               
subEventAnalysisParams.exportFig = 0;              
%subEventAnalysisParams.stimPlotParams.lineprops = [];
subEventAnalysisParams.stimDir = stimDir;

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
% frEpochsCell = {{60, @(stimDur) stimDur+60};...
%                 {-800, 60}; ...
%                 {100, 1100}}; %#ok
%               
frEpochsCell = {{60, @(stimDur) stimDur+60}};...
%                 {-800, 60}; ...
%                 {@(stimDur) stimDur+60, @(stimDur) stimDur+460}}; %#ok
              
epochLabels = {'Presentation'};%,'Fixation','Reward'};
            
assert(length(frEpochsCell) == length(epochLabels), 'Epoch time bins and epochLabel lengths must match')
%%%% note: all analysisGroups cell arrays are nx1, NOT 1xn
%Defined for Groups of 2, A-B/A+B type index.
analysisGroups.selectivityIndex.groups = {{'socialInteraction';'nonInteraction'},{'socialInteraction';'agents'}};

%Barplots showing average activity across all members of a catagory
analysisGroups.stimPrefBarPlot.groups = {{{'socialInteraction';'goalDirected';'idle';'objects';'scene';'scramble';'headTurn';'bodyTurn'}}};
analysisGroups.stimPrefBarPlot.colors  = {{{[0.55 0.13 0.16];[0.93 .2 0.15];[.98 0.65 0.13];[0 0.55 0.25];[0.15, 0.20, 0.5];[0.15, 0.20, 0.5]; [.5 0 .5]; [0 0.5 0.5]}}};
analysisGroups.stimPrefBarPlot.names = {'Barplots per label'};
analysisGroups.stimPrefBarPlot.groupDepth = 2; %2 subplots, each figure is defined by a cell array in the first item (groups).

%
analysisGroups.stimulusLabelGroups.groups = {{'socialInteraction';'goalDirected';'idle';'objects';'scene';'scramble'; 'headTurn'}};
analysisGroups.stimulusLabelGroups.names = {'Preferred Stimulus', 'Preferred Stimulus'};
analysisGroups.stimulusLabelGroups.colors = {{[0.55 0.13 0.16];[0.93 .2 0.15];[.98 0.65 0.13];[0 0.55 0.25];[0.15, 0.20, 0.5];[0.15, 0.20, 0.5]; [.5 0 .5]}, {[0.55 0.13 0.16];[0.93 .2 0.15]}};

%Essentially LFP selectivity/strength/quality
analysisGroups.evokedPotentials.groups = {{'socialInteraction';'nonInteraction';'objects'};{'socialInteraction';'nonInteraction'}};
analysisGroups.evokedPotentials.names = {'socVobj';'Social v Non-Soc'};
analysisGroups.evokedPotentials.colors = {{[0.1 0.1 1];[0.1 .7 0.1];[1 0.1 0.1]};{[0.1 0.1 1];[0.1 .7 0.1];[1 0.1 0.1]}};

%Looks like evokedpotentials, but pulls from a different analog channels
%(like pupul or eye).
analysisGroups.analogInPotentials.groups = {};
analysisGroups.analogInPotentials.channels = {[1; 2]};
analysisGroups.analogInPotentials.names = {'eyePositions,fobPlus'};
analysisGroups.analogInPotentials.units = {'degrees visual angle'};
analysisGroups.analogInPotentials.colors = {{[0.1 0.1 1];[1 0.1 0.1];'k';[0.1 .7 0.1];'m';[1 0.1 0.1];'k'}};

%Same thing, but derivatives.
analysisGroups.analogInDerivatives.groups = {};
analysisGroups.analogInDerivatives.channels = {[1; 2]};
analysisGroups.analogInDerivatives.names = {'eyeVelocity,fobPlus'};
analysisGroups.analogInDerivatives.units = {'degrees visual angle/sec'};
analysisGroups.analogInDerivatives.colors = {{[0.1 0.1 1];[0.1 .7 0.1];[1 0.1 0.1];[0.1 .7 0.1];'m';[1 0.1 0.1];'k'}};

%Makes subplots w/ PSTH on top and evoked potential on the bottom
analysisGroups.colorPsthEvoked.groups = {{'socialInteraction';'nonInteraction';'objects';'landscapes'}};
analysisGroups.colorPsthEvoked.names = {'socVobj'};
analysisGroups.colorPsthEvoked.colors = {{[0.1 0.1 1];[0 .6 0];[1 0.1 0.1]}};

%same as above, but shows error bars, harder to see catagory selectivity
%though
analysisGroups.linePsthEvoked.groups = {{'socialInteraction';'nonInteraction';'objects'}};
analysisGroups.linePsthEvoked.names = {'socVobj'};
analysisGroups.linePsthEvoked.colors = {{[0.1 0.1 1];[0 .6 0];[1 0.1 0.1]}};

%Same as above, but everything ontop of eachother in 1 panel.
analysisGroups.evokedPsthOnePane.groups = {};
analysisGroups.evokedPsthOnePane.names = {'faceVnon'};

%Creates tuning curves for units, you should have some meaningful numeric
%descriptor.
analysisGroups.tuningCurves.groups = {}; %can be images or categories
analysisGroups.tuningCurves.paramValues = {[-90 -45 0 45 90], [-90 -45 0 45 90]};
analysisGroups.tuningCurves.paramLabels = {'viewing angle (degrees)','viewing angle (degrees)'};
analysisGroups.tuningCurves.names = {'Human face view','Monkey face view'};

%Power spectra of SPIKE TIMES (i.e. 1's at spike times, flat elsewhere).
%Used for both Spikes and LFPs.
analysisGroups.spectraByCategory.groups = {{'interaction';'scrambles';}};  %todo: add spectra diff?
analysisGroups.spectraByCategory.names = {'Interaction V Scrambles'};
analysisGroups.spectraByCategory.colors = {{[1 0.1 0.1];[0.1 0.1 1]}};

%Calculates the same spectra, but w/ sliding windows. sometimes the Power
%spectrum changes overtime.
analysisGroups.tfSpectraByCategory.groups = {{'socialInteraction';'nonInteraction';}};%{'object'};{'body'}      %todo: add tf spectra diff?
analysisGroups.tfSpectraByCategory.names = {'socialInt';'Int'};%'nonface';'object';'body'

%Evoked potential plot, shows individual traces for a bunch of trials.
analysisGroups.lfpSingleTrialsByCategory.groups = {{'socialInteraction';'nonInteraction';}};
analysisGroups.lfpSingleTrialsByCategory.names = {'SocialVNonSocial'};

%Coherence between LFP time series and spike time series w/i single
%channel. Does this on a trial by trial basis, and then averages across all
%members of each group. 
analysisGroups.coherenceByCategory.groups = {{'interaction';'scrambles';}};
analysisGroups.coherenceByCategory.colors = {{[1 0.1 0.1];[0.1 0.1 1]}}; 
analysisGroups.coherenceByCategory.names = {'Interaction V Scrambles'}; %'fob';'slimCats'

%Calculates the same as above but in sliding windows.
analysisGroups.tfCouplingByCategory.groups = {{'socialInteraction'};{'nonInteraction';};{'objects'};{'goalDirected'}}; %#ok
%%%%%

if calcSwitch.useJacknife
  chronuxParams.err(1) = 2; %#ok
end

%% set paths and directories, EDIT RARELY %%
analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);
lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);        
spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
taskFilename = sprintf('%s/%s/%s%s.mat',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance 
[photodiodeFilename, lineNoiseTriggerFilename] = deal(lfpFilename);                %#ok
outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
analysisParamFilename = strcat(outDir,analysisParamFilenameStem);
preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok

% In case difference logfile is being used.
if ~logical(exist(taskFilename,'file'))
  [A, B, C] = fileparts(taskFilename);
  switch C
    case '.mat'
      taskFilename = [A '/' B '.bhv2'];
    case '.bhv2'
      taskFilename = [A '/' B '.mat'];
  end
end

% Check if the key file types exist.
assert(logical(exist(analogInFilename,'file')),'The analog input file you requested does not exist.');
assert(logical(exist(taskFilename,'file')),'The log file you requested does not exist.');

% Auto channel detect - if set, find the parsed file (or parse it), and
% reassign spike channels on this.
if autoChannelDetect
  %If the file is parsed, retrieve channels present
  parsedFolderName = sprintf('%s/%s/%s%s_parsed',ephysVolume,dateSubject,dateSubject,runNum);
  [ephysParams.spikeChannels, ephysParams.lfpChannels, ephysParams.channelNames] = autoDetectChannels(parsedFolderName);
end

if ~exist(outDir,'dir')
  mkdir(outDir);
end
if isempty(varargin) 
  save(analysisParamFilename);
elseif strcmp(varargin,'saveNoPreprocParams')
  save(analysisParamFilename,'calcSwitch','analysisGroups','plotSwitch','-append');
end
end

function swappedString = slashSwap(pathString)
%Swaps direction of slashes to match Unix/Phyzzy, from Windows Path.
  stringParts = split(pathString, '\');
  swappedString = char(join(stringParts, '/'));
end

function [spike, LFP, names] = autoDetectChannels(parsedFolderName)
  if exist(parsedFolderName,'dir') == 7
    channelFiles = dir([parsedFolderName '/*.NC5']);
    channelNames = {channelFiles.name};
    channelsTmp = [];
    channelNameTmp = {};
    for chan_ind = 1:length(channelNames)
      startInd = regexp(channelNames{chan_ind}, 'Ch') + 2;
      stopInd = regexp(channelNames{chan_ind},'.NC5') - 1;
      chanNum = channelNames{chan_ind}(startInd:stopInd);
      channelNameTmp{chan_ind} = ['Ch' chanNum];
      channelsTmp(chan_ind) = str2double(channelNames{chan_ind}(startInd:stopInd));
    end
    [spike, LFP] = deal(channelsTmp); %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
    names = channelNameTmp;
  else
  end
end