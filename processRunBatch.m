function [] = processRunBatch(varargin)
% applies the processRun function to a set of analysisParamFiles created
% with differente dateSubj and run variables
%   Inputs:
%   - varargin can have the following forms:
%     - two arguments, the first being a "runList" cell array of the form 
%       {{'dateSubj';{'001','002'}} i.e.{{'180314MO'; {'002';'003'}}}, the
%       second of the form paramBuilderFunctionName as a string i.e.
%       'buildAnalysisParamFileSocVid'
%     - one argument of the form paramBuilderFunctionName as a string i.e.
%       'buildAnalysisParamFileSocVid'. User will be prompted to select
%       directory with files.

%% Load Appropriate variables
if nargin == 1
    analysisParamFile = varargin{1};
    disp('Select Folder containing Data')
    runList = buildRunList(uigetdir);
elseif nargin == 2
    runList = varargin{1};
    analysisParamFile = varargin{2};
elseif nargin >= 3
    disp('Too many input arguments.')
    return
end

load(feval(analysisParamFile));
for dateSubj_i = 1:length(runList)
    dateSubject = runList{dateSubj_i}{1};
    for run_i = 1:length(runList{dateSubj_i}{2})
        runNum = runList{dateSubj_i}{2}{run_i};
        analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);   %#ok
        lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);        
        spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
        taskFilename = sprintf('%s/%s/%s%s.bhv2',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
        outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
        stimSyncParams.outDir = outDir;
        ephysParams.outDir = outDir;
        photodiodeFilename = lfpFilename;                %#ok
        lineNoiseTriggerFilename = lfpFilename; %#ok        
        analysisParamFilename = strcat(outDir,analysisParamFilenameStem);
        preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok
        if ~exist(outDir,'dir')
          mkdir(outDir);
        end
        save(analysisParamFilename);
        fprintf('Processing run for %s%s... \n', dateSubject,runNum);
        try
            processRun('paramFile', analysisParamFilename)
        catch
            disp('Error caught, going on to next loop')
            clc; close all;
            fprintf('Ran into Error, Moving to next File. \n');
            continue
        end
        clc; close all;
        fprintf('Done! \n');
    end
end
end

