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

%% Load Appropriate variables and paths
addpath('buildAnalysisParamFileLib');
addpath(genpath('dependencies'));

if nargin == 1
    analysisParamFile = varargin{1};
    disp('Select Folder containing Data')
    runList = buildRunList(uigetdir);
elseif nargin == 2
    runList = varargin{1};
    analysisParamFile = varargin{2};
elseif nargin >= 3 || nargin == 0
    disp('Must have either 1 or 2 inputs.')
    return
end

%% Create all of the appropriate AnalysisParamFiles as a cell array.
load(feval(analysisParamFile));
rmdir(outDir, 's'); %The analysis file is overwritten, so if it isn't a run you're about to do, this file is no longer correct. 
analysisParamFileList = cell(0);
meta_ind = 1;

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
        analysisParamFileList{meta_ind} = analysisParamFilename;
        meta_ind = meta_ind + 1;
    end
end

%% Process the runs
[errorsMsg, startTimes, endTimes] = deal(cell(length(analysisParamFileList),1));

if license('test','Distrib_Computing_Toolbox')
  parfor run_ind = 1:length(analysisParamFileList)
    fprintf('Processing %s... \n', analysisParamFileList{run_ind});
    startTimes{run_ind} = datetime('now');
    try
      processRun('paramFile', analysisParamFileList{run_ind})
      errorsMsg{run_ind} = 'None';
      endTimes{run_ind} = datetime('now');
    catch MyErr
      errorsMsg{run_ind} = MyErr.message;
      endTimes{run_ind} = datetime('now');
      fprintf('Ran into Error, Moving to next File. \n');
      continue
    end
    clc; close all;
    fprintf('Done! \n');
  end
else
  for run_ind = 1:length(analysisParamFileList)
    fprintf('Processing %s... \n', analysisParamFileList{run_ind});
    startTimes{run_ind} = datetime('now');
    try
      processRun('paramFile', analysisParamFileList{run_ind})
      errorsMsg{run_ind} = 'None';
      endTimes{run_ind} = datetime('now');
    catch MyErr
      errorsMsg{run_ind} = MyErr.message;
      endTimes{run_ind} = datetime('now');
      fprintf('Ran into Error, Moving to next File. \n');
      continue
    end
    clc; close all;
    fprintf('Done! \n');
  end
end

%Save Batch Run Results
errorReport = [analysisParamFileList' errorsMsg startTimes endTimes];
T = cell2table(errorReport);
T.Properties.VariableNames = {'File_Analyzed', 'Error', 'Start_Time', 'End_time'};
writetable(T,sprintf('%s/BatchRunResults.csv',outputVolume))

end
