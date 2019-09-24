function [] = processRunBatch(varargin)
% applies the processRun function to a set of analysisParamFiles created
% with differente dateSubj and run variables
%   Inputs:
%   - varargin can have the following forms:
%     - two arguments, the first being a "runList" cell array of the form 
%       {{'dateSubj';{'001','002'}} i.e.{{'180314Mo'; {'002';'003'}}}, the
%       second of the form paramBuilderFunctionName as a string i.e.
%       'buildAnalysisParamFileSocVid'
%     - one argument of the form paramBuilderFunctionName as a string i.e.
%       'buildAnalysisParamFileSocVid'. User will be prompted to select
%       directory with files.

%% Load Appropriate variables and paths
addpath(genpath('D:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy'))
% addpath('buildAnalysisParamFileLib');
% addpath(genpath('dependencies'));
    
if nargin == 1
    %If there is only 1 file, it loads the analysisParamFile and composes a
    %list from all the data files in the ephysVolume.
    analysisParamFile = varargin{1};
    load(feval(analysisParamFile));
    runListFolder = uigetdir;
    if runListFolder == 0
      disp('using ephysVolume as directory')
      runListFolder = ephysVolume;
    end
    runList = buildRunList(runListFolder, 'nev');
elseif nargin == 2
    runList = varargin{1};
    analysisParamFile = varargin{2};
elseif nargin >= 3 || nargin == 0
    disp('Must have either 1 or 2 inputs.')
    return
end

%% Create all of the appropriate AnalysisParamFiles as a cell array.
[analysisParamFileList, analysisParamFileName] = deal(cell(0));
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
        analysisParamFileName{meta_ind} = [dateSubject runNum];
        meta_ind = meta_ind + 1;
    end
end

%% Process the runs
[errorsMsg, startTimes, endTimes, analysisOutFilename] = deal(cell(length(analysisParamFileList),1));

diary batchRunLog.txt
diary on

analysisParamFileList'

if license('test','Distrib_Computing_Toolbox')
  parfor run_ind = 1:length(analysisParamFileList)
    fprintf('run_ind reads %d... \n', run_ind);
    fprintf('Processing %s... \n', analysisParamFileList{run_ind});
    startTimes{run_ind} = datetime('now');
    try
      [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
      errorsMsg{run_ind} = 'None';
      endTimes{run_ind} = datetime('now');
    catch MyErr
      errorsMsg{run_ind} = MyErr.message;
      endTimes{run_ind} = datetime('now');
      fprintf('Ran into Error, Moving to next File. \n');
      continue
    end
    close all;
    fprintf('Done! \n');
  end
else
  for run_ind = 1:length(analysisParamFileList)
    fprintf('Processing %s... \n', analysisParamFileList{run_ind});
    startTimes{run_ind} = datetime('now');
    try
      [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
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

diary off

%analysisOutFilename is now a cell array of the filepaths to the outputs of
%runAnalyses. Cycle through them and extract desired information (# of
%units, significance), which you can add to the output file.
[UnitCount, sigUnits, sigStim, sigStimLen] = deal(cell(size(analysisOutFilename)));

for ii = 1:length(analysisOutFilename)
  if ~isempty((analysisOutFilename{ii}))
    tmp = load(analysisOutFilename{ii}, 'sigStruct'); %Relies on psth Overlay function in runAnalyses.
    UnitCount{ii} = tmp.sigStruct.totalUnits{1};
    sigUnits{ii} = length(tmp.sigStruct.sigUnits);
    sigStim{ii} = unique(cat(1, tmp.sigStruct.sigStim{:}));
    sigStimLen{ii} = length(sigStim{ii});
    if ~isempty(sigStim{ii})
      sigStim{ii} = strjoin(sigStim{ii}, ' ');
    end
  end
end

% %Save Batch Run Results
errorReport = [analysisParamFileName' startTimes endTimes errorsMsg UnitCount sigUnits sigStimLen sigStim ];
T = cell2table(errorReport);
T.Properties.VariableNames = {'File_Analyzed', 'Start_Time', 'End_time', 'Error', 'Unit_Count', 'Signifiant_Unit_count', 'Stimuli_count', 'Stimuli'};
writetable(T,sprintf('%s/BatchRunResults.csv',outputVolume))

%Save Batch Run Results
% errorReport = [analysisParamFileName' startTimes endTimes errorsMsg];
% T = cell2table(errorReport);
% T.Properties.VariableNames = {'File_Analyzed', 'Start_Time', 'End_time', 'Error'};
% writetable(T,sprintf('%s/BatchRunResults.csv',outputVolume))

%PDF Summary
%Remove empty cells from analysisOutFilename
nonEmptyCellInd = ~(cellfun('isempty',analysisOutFilename));
analysisOutFilename = analysisOutFilename(nonEmptyCellInd);

%Find the directory
figDirPath = @(cell) (fileparts(cell));
analysisOutFigDirs = cellfun(figDirPath, analysisOutFilename, 'UniformOutput', 0);

%run through the createSummaryDoc Function.
for figDir_ind = 1:length(analysisOutFigDirs)
  createSummaryDoc('buildLayoutParamFile', analysisOutFigDirs{figDir_ind})
end


end
