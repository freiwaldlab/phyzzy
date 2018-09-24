function [ analysisParamFilename ] = buildAnalysisParamFileMinimal( )    
%buildAnalysisParamFile saves a mat file of parameters, which control the
%behavior of analyzeSession
%   todo: option to load 'fixed' params from file, for ease accross days


%%%%%%%  USER PARAMETERS, EDIT ROUTINELY %%%%%%%%
runNum = '001';
dateSubject = '171026ALAN'; 
ephysVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Blackrock'; 
stimulusLogVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Visiko';
outputVolume = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Analyzed';
stimParamsFilename = '/Users/stephenserene/Desktop/Freiwald/AnalysisSerene/StimParamsFullFOB3.mat';   %#ok
analysisLabel = 'Minimal';
analysisParamFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'
preprocessedDataFilenameStem = 'preprocessedData.mat';


% parameters preprocessSpikes and preprocessLFP, see functions for details
% note: ephysParams is also optional if you make a runAnalyses file that does not use spikes of LFPs, however that functionality isn't fully tested as of 10/26/2017
ephysParams.needLFP = 1;
ephysParams.needSpikes = 1;
ephysParams.spikeChannels = [1]; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
ephysParams.lfpChannels = [1]; 
ephysParams.cPtCal = 1/30; % conversion from spike sample indices to timestep of decimated LFP


%%% set paths and directories, EDIT RARELY %%%
analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);   %#ok
lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);        
spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
taskFilename = sprintf('%s/%s/%s0%s.log',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
analysisParamFilename = strcat(outDir,analysisParamFilenameStem);
preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok
%
if ~exist(outDir,'dir')
  mkdir(outDir);
end
save(analysisParamFilename);
end

