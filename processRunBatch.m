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
usePreprocessed = 1;

if nargin == 1
    %If there is only 1 file, it loads the analysisParamFile and composes a
    %list from all the data files in the ephysVolume.
    analysisParamFile = varargin{1};
    load(feval(analysisParamFile));
    runListFolder = ephysVolume;
    runList = buildRunList(runListFolder, 'nev');
elseif nargin >= 2 || nargin == 0
    disp('Must have 1 input.')
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
      %Autodetecting spike channels
      %If the file is parsed, retrieve channels present
      parsedFolderName = sprintf('%s/%s/%s%s_parsed',ephysVolume,dateSubject,dateSubject,runNum);
      if exist(parsedFolderName,'dir') == 7
        [ephysParams.spikeChannels, ephysParams.lfpChannels, ephysParams.channelNames] = autoDetectChannels(parsedFolderName);
      end
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
if usePreprocessed
  for path_ind = 1:length(analysisParamFileList)
    analysisParamFileList{path_ind} = [fileparts(analysisParamFileList{path_ind}) filesep 'preprocessedData.mat'];
  end
end

if license('test','Distrib_Computing_Toolbox')
  parfor run_ind = 1:length(analysisParamFileList)
    fprintf('run_ind reads %d... \n', run_ind);
    fprintf('Processing %s... \n', analysisParamFileList{run_ind});
    startTimes{run_ind} = datetime('now');
    try
      if usePreprocessed
        [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids','preprocessed',analysisParamFileList{run_ind});
      else
        [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
      end
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
      if usePreprocessed
        [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids','preprocessed',analysisParamFileList{run_ind});
      else
        [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
      end
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
%%
%analysisOutFilename is now a cell array of the filepaths to the outputs of
%runAnalyses. Cycle through them and extract desired information (# of
%units, significance), which you can add to the output file.
%[UnitCount, sigUnits, sigStim, sigStimLen] = deal(cell(size(analysisOutFilename)));
replaceAnalysisOut = 1;
if replaceAnalysisOut
  addEnd = @(x) strjoin([x, {'analyzedData.mat'}], filesep);
  breakString = @(x) strsplit(x, filesep);
  joinStrings = @(x) strjoin([x(length(x)-2), x(length(x))],'');
  
  outputVolume = 'D:\DataAnalysis\FullTime';
  analysisOutFilename = dir([outputVolume '\**\analyzedData.mat']);
  analysisOutFilename = {analysisOutFilename.folder}';
  [errorsMsg, startTimes, endTimes] = deal(cell(length(analysisOutFilename),1));
  
  tmpAnalysisParamFileName = cellfun(breakString, analysisOutFilename, 'UniformOutput',0);
  analysisParamFileName = cellfun(joinStrings, tmpAnalysisParamFileName, 'UniformOutput',0);
  analysisOutFilename = cellfun(addEnd, analysisOutFilename,'UniformOutput',0);
  [startTimes{:}, endTimes{:}] = deal(datetime(now,'ConvertFrom','datenum'));
  [errorsMsg{:}] = deal('None');
end
titles = {'File_Analyzed', 'Start_Time', 'End_time', 'Error', 'Channel', 'Unit_Count', 'Sig Unit count', 'Sig Unsorted','Sig MUA', 'Stimuli_count', 'Stimuli','ANOVA','Other Info'};
sigMUAInd = strcmp(titles,'Sig MUA');
sigUnsortedInd = strcmp(titles,'Sig Unsorted');

epochCount = 3; %Hardcoding 3 Epochs here.
table = cell(epochCount,1);
[table{:}] = deal(cell(length(analysisOutFilename), length(titles)));
true_ind = 1;

for ii = 1:length(analysisOutFilename)
  if ~isempty((analysisOutFilename{ii}))
    tmp = load(analysisOutFilename{ii}, 'stimCatANOVATable','sigStruct','frEpochs'); %Relies on psth Overlay function in runAnalyses.
    true_ind_page = true_ind;
    for epoch_i = 1:length(tmp.sigStruct.sigUnits)
      for channel_ind = 1:length(tmp.sigStruct.channels) %Add channel count here.
        % Null model based significance testing
        channel = tmp.sigStruct.channels(channel_ind);
        UnitCount = tmp.sigStruct.totalUnits{channel_ind};
        sigUnits = length(tmp.sigStruct.sigUnits{epoch_i}{channel_ind});
        sigUnsorted = tmp.sigStruct.sigUnsorted{epoch_i}{channel_ind};
        sigMUA = tmp.sigStruct.sigMUA{epoch_i}{channel_ind};
        sigStim = unique(cat(1, tmp.sigStruct.sigStim{epoch_i}{channel_ind}));
        sigStimLen = length(sigStim);
        if ~isempty(sigStim)
          sigStimNames = strjoin(sigStim, ' ');
        else
          sigStimNames = ' ';
        end
        % ANOVA based significance
        unitStructs = tmp.stimCatANOVATable{channel_ind};
        ANOVASigString = ' ';
        for unit_ind = 1:length(unitStructs)
          anovaSig = 0;
          if isfield(unitStructs{unit_ind}, 'ANOVA')
            if epoch_i == 1
              anovaSig = unitStructs{unit_ind}.ANOVA.pp < 0.05;
            elseif epoch_i == 2
              anovaSig = unitStructs{unit_ind}.ANOVA.fp < 0.05;
            elseif epoch_i == 3
              anovaSig = unitStructs{unit_ind}.ANOVA.rp < 0.05;              
            end
          end
          ANOVASigString = [ANOVASigString ['[' num2str(unitStructs{unit_ind}.taskModulatedP < 0.05) ';' num2str(anovaSig(1)) ']']];
        end
        %Package
        table{epoch_i}(true_ind, :) = [analysisParamFileName(ii), startTimes(ii), endTimes(ii), errorsMsg(ii), channel, UnitCount, sigUnits, sigUnsorted, sigMUA, sigStimLen, sigStimNames, ANOVASigString, ' '];
        if ii == 1 && epoch_i == 1
%           table{epoch_i}{true_ind, 11} = sprintf('Comparison: %s Vs %s', strjoin(unitStructs{1}.ANOVA.stats.grpnames{1}), strjoin(unitStructs{1}.ANOVA.stats.grpnames{2}));
        end
        true_ind = true_ind + 1;
      end
      if epoch_i ~= length(tmp.sigStruct.sigUnits) % More pages to do = let count reset.
        true_ind = true_ind_page;
      end
    end
  end
end

%Save Batch Run Results
for table_ind = 1:length(table)
  table{table_ind}{3,end} = sprintf('%d - %d ms', tmp.frEpochs(table_ind,1), tmp.frEpochs(table_ind,2));
  T = cell2table(table{table_ind});
  T.Properties.VariableNames = titles;
  writetable(T,sprintf('%s/BatchRunResults.xlsx',outputVolume),'Sheet', sprintf('%s Epoch', tmp.sigStruct.epochLabels{table_ind}))
end
% %PDF Summary
% %Remove empty cells from analysisOutFilename
% nonEmptyCellInd = ~(cellfun('isempty',analysisOutFilename));
% analysisOutFilename = analysisOutFilename(nonEmptyCellInd);
% 
% %Find the directory
% figDirPath = @(cell) (fileparts(cell));
% analysisOutFigDirs = cellfun(figDirPath, analysisOutFilename, 'UniformOutput', 0);
% 
% %run through the createSummaryDoc Function.
% for figDir_ind = 1:length(analysisOutFigDirs)
%   createSummaryDoc('buildLayoutParamFile', analysisOutFigDirs{figDir_ind})
% end


end

function [spike, LFP, names] = autoDetectChannels(parsedFolderName)
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
end
