%% Wave clus bulk preprocess
dataDir = 'C:/Data 2018';

runList = buildRunList(dataDir);

for file_i = 1:length(runList)
  for run_i = 1:size(runList{file_i}{2}, 2)
    ns5Path = [dataDir '\' runList{file_i}{1} '\' runList{file_i}{1} runList{file_i}{2}{run_i} '.ns5'];
    parse_data_NSx(ns5Path, 2)
  end
end


function [runList] = buildRunList(dataDir)
% Constructs the runList variable for the processRunBatch function based on
% the selection of a directory.
runList_struct = dir(dataDir);
runList_cell = struct2cell(runList_struct);

runList = runList_cell(1, cell2mat(runList_cell(5,1:end)));

%Assumes the first 2 entries are '.' and '..', and discards them. This is
%not an ideal method and may not work on non-Windows.
runList = runList(3:end);

runList_tmp = cell(1, length(runList));

for ii = 1:length(runList)
    filePattern = fullfile([dataDir,filesep,runList{1,ii}], '*.nev');
    theFiles = struct2cell(dir(filePattern));
    theFiles = theFiles(1,:);
    cells = cellfun(@(x) extractBetween(x, runList{1,ii}, '.'), theFiles);
    runListBlock = {runList{ii}; cells};
    runList_tmp{ii} = runListBlock;    
end

runList = runList_tmp;
end