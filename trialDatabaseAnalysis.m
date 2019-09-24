%% trialDatabase Analysis

analysisDir = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\Analyzed';

% Find the available trials in a Analysis directory.
trialDatabaseFiles = dir([analysisDir '\**\trialDatabase.mat']);

% Load the information from these files
trialDatabaseData = cell(length(trialDatabaseFiles),1);
for file_ind = 1:length(trialDatabaseFiles)
  trialDatabaseData{file_ind} = load([trialDatabaseFiles(file_ind).folder filesep trialDatabaseFiles(file_ind).name]);
end

% Different stimuli used across days will cause issues here, so we need to
% focus on stimuli which all neurons saw. Find the stimuli shared by the
% largest number of units.
%allStim = cell(50, length(trialDatabaseData));
for ind = 1:length(trialDatabaseData)
    stimInRun = trialDatabaseData{ind}.trialDatabaseStruct.images;
    allStim(1:length(stimInRun),ind) = stimInRun;
end

tf = cellfun('isempty',allStim);      % true for empty cells
allStim(tf) = {'None'};               % replace by a cell with a 'none'
uniqueStim = unique(allStim);

bigDataStimInd = cellfun(@(c)strcmp(c,allStim),uniqueStim,'UniformOutput',false);
combinedDataMUAMat = zeros(length(uniqueStim), length(trialDatabaseData));

%create a matrix where each item is a spike count for a particular stimuli,
%by a unit.
epoch_ind = 1;
for run_ind = 1:length(trialDatabaseData)
  runData = trialDatabaseData{run_ind}.trialDatabaseStruct;
  %Find the index to use to sort trial data into larger matrix.
  sortingInd = cellfun(@(c)find(strcmp(c, uniqueStim)), runData.images,'UniformOutput',false);
  sortingInd = [sortingInd{:}]'; %Reshape and remove from cells
  %extract data corresponding to epoch of interest (mean MUA for now).
  MUAData = runData.spikeData{1}{1}.MUA;
  meanMUACounts = cellfun(@(c)mean(c.counts), MUAData,'UniformOutput',false); %Take the mean MUA across trials.
  meanMUACounts = [meanMUACounts{:}]'; %Reshape and remove from cells
  combinedDataMUAMat(sortingInd,run_ind) = meanMUACounts;
end

% See how many counts exist per stim, to create a more compact space.
stimCount = zeros(size(combinedDataMUAMat, 1),1);
for row_ind = 1:size(combinedDataMUAMat, 1)
  stimCount(row_ind) = nnz(combinedDataMUAMat(row_ind,:));
end

