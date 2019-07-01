%% Wave clus bulk preprocess
% A script which searches a directory for files of a type denoted below,
% and then applies parse_data_NSx to them.

%% Variables
openNSxPath = 'D:/Onedrive/Lab/ESIN_Ephys_Files/Analysis/phyzzy/dependencies/NPMK';
dataDir = 'C:/DataSlice';
RAM = 2; %GB of RAM allocated to process - check parse_data_NSX for more.
fileExt = '.ns5';

%Setting relevant to deleting unpackaged channels which aren't of interest.
deleteCh = 1;
relevantCh = 128; %The Threshold of Relevant channels. Blackrock Machines used in our lab have 1 - 128 channels.


%% Script
addpath(openNSxPath);
runList = buildRunList(dataDir, fileExt);
parsedData = cell(length(runList));

for file_i = 1:length(runList)
  for run_i = 1:size(runList{file_i}{2}, 2)
    filePath = [dataDir filesep runList{file_i}{1} filesep runList{file_i}{1} runList{file_i}{2}{run_i} fileExt];
    parsedData = parse_data_NSx(filePath, RAM);
    if deleteCh
      for data_i = 1:length(parsedData)
        [~, B, ~] = fileparts(parsedData{data_i});
        B = split(B,'_');
        chNum = str2double(B{3}(3:end));
        if chNum < relevantCh
          continue
        else
          delete(parsedData{data_i})
          parsedData{data_i} = [];
        end
      end
    end
  end
end

%% Individual Functions
function [runList] = buildRunList(dataDir, fileExt)
% Constructs the runList variable for the bulk_parse_data function based on
% the selection of a directory.
runList_struct = dir(dataDir);
runList_cell = struct2cell(runList_struct);

runList = runList_cell(1, cell2mat(runList_cell(5,1:end)));

%Assumes the first 2 entries are '.' and '..', and discards them. This is
%not an ideal method and may not work on non-Windows.
runList = runList(3:end);

runList_tmp = cell(1, length(runList));

for ii = 1:length(runList)
    filePattern = fullfile([dataDir filesep runList{1,ii}], ['*' fileExt]);
    theFiles = struct2cell(dir(filePattern));
    theFiles = theFiles(1,:);
    cells = cellfun(@(x) extractBetween(x, runList{1,ii}, '.'), theFiles);
    runListBlock = {runList{ii}; cells};
    runList_tmp{ii} = runListBlock;    
end

runList = runList_tmp;
end