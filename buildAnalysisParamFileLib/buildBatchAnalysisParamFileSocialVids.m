function [batchAnalysisParamFilename] = buildBatchAnalysisParamFileSocialVids( varargin )
% Generate a file which specifies the parameters for batch analysis. 

[~, machine] = system('hostname');
machine = machine(~isspace(machine));

switch machine
  case 'Alienware_FA'
    analysisDirectory = slashSwap('D:\DataAnalysis\FullTime');
    outputDir = [analysisDirectory '/batchAnalysis'];
    stimParamsFilename = slashSwap('D:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
  case 'HomeDesktop'
    analysisDirectory = slashSwap('E:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\AnalysisSample');
    outputDir = [analysisDirectory '/batchAnalysis'];
    stimParamsFilename = slashSwap('E:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
end

analysisLabel = 'Basic';
preprocessedDataFilenameStem = 'preprocessedData.mat';
analysisParamFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'

saveFig = 1;                
closeFig = 0;               %#ok
exportFig = 0;              %#ok
saveFigData = 0;            %#ok
verbosity = 'INFO';         %other options, 'DEBUG', 'VERBOSE';

%% Switches
calcSwitch.excludeRepeats = 1;
calcSwitch.cellCount = 1;

plotSwitch.stimPresCount = 1;
plotSwitch.slidingWindowANOVA = 1;

%% Parameters
preprocessParams.spikeDataFileName = 'spikeDataBank'; %File ending in .mat, not included to allow for slicing (e.g. 'spikeDataBank_1.mat'...)
%3 variables below are not used, as defining the variables in a load
%command w/ a cell array (or other structures) doesn't work.
preprocessParams.preprocessedVars = {'spikesByEvent','eventIDs','eventCategories','preAlign','postAlign'}; %Variables extracted from preprocessedData.mat
preprocessParams.analyzedVars = {'analysisParamFilename','dateSubject', 'runNum', 'groupLabelsByImage','psthByImage','attendedObjData'}; %Variables extracted from analyzedData.mat
preprocessParams.analysisParamVars = {'psthParams'}; %Variables extracted from analysisParam.mat

cellCountParams.excludePhase2 = 1; % a switch which can be used to remove data from the same neuron collected in subsequent runs. Good for getting accurate counts.
cellCountParams.batchRunxls = fullfile(analysisDirectory,'BatchRunResults.xlsx');                         %Batch analysis xlsx produced by processRunBatch.
cellCountParams.recordingLogxls = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Data\RecordingsMoUpdated.xlsx';  %Used to exclude phase 2 to give accurate unit counts.

slidingANOVAParams.binSize = 100;
slidingANOVAParams.binStep = 25;
slidingANOVAParams.Omega = 1; %switch to convert ANOVA curves to Omega curves.
slidingANOVAParams.target = {'socialInteraction','agents','interaction'};  %Labels which must exist in the stimParamFile associated with the runs. 
slidingANOVAParams.stimParamFile = stimParamsFilename;

psthParams.type = 'normal'; %options are 'normal', 'baselineSub', 'meanWhite'
psthParams.psthPre = 800; % if e.g. +200, then start psth 200ms before trial onset; 
psthParams.psthImDur = 2800;  % only need to set this for variable length stim runs; else, comes from log file
psthParams.psthPost = 500;
psthParams.latency = 0;
psthParams.movingWin = [200 5];
psthParams.smoothingWidth = 25;  %psth smoothing width, in ms
psthParams.Zscore = 0;  % 0: raw PSTH, 1: pre-trial baseline subtraction Z Scored PSTH, 2: whole trial baseline subtracted Z Scored PSTH
psthParams.errorType = 1; %chronux convention: 1 is poisfStimson, 2 is trialwise bootstrap, 3 is across trial std for binned spikes, bootstrap for spike times 
psthParams.errorRangeZ = 1; %how many standard errors to show
psthParams.bootstrapSamples = 100;
psthParams.sortStim = 1;
psthParams.sortOrder = {'socialInteraction';'goalDirected';'idle';'objects';'scene';'scramble';'allStim'};
psthParams.psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'
load(psthParams.psthColormapFilename);
psthParams.colormap = map;
psthParams.stimParamFile = stimParamsFilename;

analysisParamFilenameStem = 'batchAnalysisParams.mat';
if ~exist([outputDir '/'])
  mkdir([outputDir '/']);
end
batchAnalysisParamFilename = [outputDir '/' analysisParamFilenameStem];
save(batchAnalysisParamFilename)

end

function swappedString = slashSwap(pathString)
%Swaps direction of slashes to match Unix/Phyzzy, from Windows Path.
  stringParts = split(pathString, '\');
  swappedString = char(join(stringParts, '/'));
end