function [ analysisParamsFilename ] = buildAnalysisParamFileRISE( )
%buildAnalysisParamFile saves a mat file of parameters, which control the
%behavior of analyzeSession
%   todo: option to load 'fixed' params from file, for ease accross days


%%%%%%%  USER PARAMETERS, EDIT ROUTINELY %%%%%%%%
runNum = '007';
dateSubject = '161011ALAN';
ephysVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Blackrock'; %'/Volumes/Users-1/User/Desktop%' %
stimulusLogVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Visiko';
outputVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Analyzed';
analysisLabel = 'Basic';
analysisParamsFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'
picParamsFilename = '/Users/stephenserene/Desktop/Freiwald/AnalysisSerene/StimParamsRISE.mat';
categoryListSlim = {'humanFace','monkeyFace','place','fruit','humanBody','monkeyBody','techno'}; %minimal cat list for clean plots
saveFig = 0;
exportFig = 0;
saveFigData = 0;
verbosity = 'INFO'; %other options, 'DEBUG', 'VERBOSE';

% parameters preprocessSpikes and preprocessLFP, see functions for details
ephysParams.needLFP = 1;
ephysParams.needSpikes = 1;
ephysParams.spikeChannels = [1,33,34]; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
ephysParams.lfpChannels = [1,33,34]; 
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
ephysParams.filter = ''; % if filtering desired, ephysFilter is a digitalFilter 

% parameters preprocessAnalogIn, see function for details
analogInParams.needAccel = 0;
analogInParams.needEyes = 0;
analogInParams.accelChannels = []; %2d array, numAccelerometers x numChannelsPerAccelerometer
analogInParams.accelChannelNames = {}; % row vector, whose length matches num rows in accelChannels
analogInParams.eyeChannels = [];

% parameters preprocessLogFile, see function for details
stimSyncParams.usePhotodiode = 0;

% parameters for excludeStimuli, see function for details
excludeStimParams.fixPre = 200; %ms
excludeStimParams.fixPost = 200; %ms
excludeStimParams.flashPre = 200;  %ms
excludeStimParams.flashPost = 200; %ms
excludeStimParams.juicePre = 500; % optional, ms
excludeStimParams.juicePost = 500; % optional, ms
excludeStimParams.DEBUG = 0; % makes exclusion criterion plots if true
% additional optional excludeStimParams: accel1, accel2, minStimDur (ms)


psthParams.psthPre = 100; % if e.g. +200, then start psth 200ms before trial onset; 
psthParams.psthImDur = 0;  % only need to set this for variable length stim runs; else, comes from log file
psthParams.psthPost = 200;
psthParams.smoothingWidth = 10; %psth smoothing width, in ms
% TW=3 with T=.2, then W = 15 Hz (5 tapers)
% TW=1.5 with T=.1, then W = 15 Hz (2 tapers)
% TW = 1.5 with T=.2, then W = 7.5 Hz (2 tapers)
chr_params.tapers = [1 1]; %[3 5] is chronux default; 
chr_params.pad = 1;
chr_params.fs = 1;
chr_params.trialave = 1;
chr_params.err = [1 .05];
chr_params.fpass = [0 .1];
tfParams.movingWin = [50 5];
tfParams.specgramRowAve = 0;

lfpAlignParams.samPerMS = 1; % because this is after decimation
lfpAlignParams.msPreAlign = psthParams.psthPre+tfParams.movingWin(1)/2; 
lfpAlignParams.msPostAlign = psthParams.psthImDur+psthParams.psthPost+tfParams.movingWin(1)/2;
spikeAlignParams.preAlign = psthParams.psthPre+3*psthParams.smoothingWidth;
spikeAlignParams.postAlign = psthParams.psthImDur+psthParams.psthPost+3*psthParams.smoothingWidth;
% for lfps, constrain first and (optional) last [n m] samples to 0 mean
lfpAlignParams.DCSUB_SAM = [lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10;lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10 ]; % for no DC sub, [0 ...], for both ends, [sam sam], for start only [sam 0] 
%
frCalcOn = 80;
frCalcOff = 0; %note: if frCalcOff < frCalcOn, will update when psthImDur is set after reading log file 
frCalcOnEarly = 80;
frCalcOffEarly = 250;
frCalcOnLate = 250;
frCalcOffLate = 425;
%
psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'

% Boolean variables to specify which computations to perform; TODO: read
% from config file, eventually with conditional on log file info
makeImPSTH = 1;
makeCatPSTH = 0;
imageTF = 1;
catTF = 0;
calcLatencyRF = 0;
crossTF = 1;
calcEvokedPowerRF = 0;
calcCoherenceRFcpt = 0;
calcCoherenceRFcc = 0;
calcCoherenceRFptpt = 0;
calcGrangerRF = 0;

%%% set paths and directories, EDIT RARELY %%%
analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);
lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);
spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
taskFilename = sprintf('%s/%s/%s0%s.log',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
outDir = strcat(sprintf(strcat(outputVolume,'/%s/'),dateSubject),analysisLabel,'/');
analysisParamsFilename = strcat(outDir,analysisParamsFilenameStem);
%
load('cocode2.mat');
psthColormap = map;
%
if ~exist(outDir,'dir')
  mkdir(outDir);
end
buildStimParamFileRISE(taskFilename, picParamsFilename)
save(analysisParamsFilename);
end

