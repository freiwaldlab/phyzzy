function [dirList] = buildRunList(dataDir, tag)
% Constructs the runList variable for the processRunBatch function based on
% the selection of a directory and a tag, which both denotes the directory
% contents of interest but also changes output formating.

filesOfInterest = dir([dataDir filesep '**' filesep '*.' tag]);

%Allows for customization of list
filesOfInterest = uiDirChoose(filesOfInterest);

fileDirs = unique({filesOfInterest.folder}');

switch tag
  case 'nev'
    %Find everywhere with data files - denoted by folders containing a .nev.
    %Change if not using Blackrock NSP.
    for ii = 1:length(fileDirs)
      [~, fileDirs{ii}] = fileparts(fileDirs{ii});
    end
    
    %Cycle through those folders
    dirList_tmp = cell(1, length(fileDirs));
    
    name = regexp(fileDirs{1},'\D*','Match'); %Find the name
    extractRunNum = @(x) extractBetween(x, name, '.');
    
    for ii = 1:length(fileDirs)
      theFiles = dir([dataDir filesep fileDirs{ii} filesep '*.nev']);
      theFiles = {theFiles(:).name};
      cells = cellfun(extractRunNum, theFiles);
      runListBlock = {fileDirs{ii}; cells};
      dirList_tmp{ii} = runListBlock;
    end
  case 'fig'
    %Returns directories where .fig files are located, used in conjunction
    %with the 'createSummaryDoc' function.
end
dirList = dirList_tmp;
end

function filesOfInterestUpdated = uiDirChoose(filesOfInterest)
theList = {filesOfInterest.name}';
listCheckFig = figure();
listCheckFig.Position(3) = 200;

listCheck = uicontrol(listCheckFig,'style', 'listbox','String',theList,'Value',1);
listCheck.Units = 'Normalized';
listCheck.Position = [.055 .05 0.9 0.8];
listCheck.Max = length(theList);

selectButton = uicontrol(listCheckFig,'style','pushbutton','String','Select','Units','Normalized','Position',[.055 .9 0.3 0.1]);
selectButton.Callback = @resumefig;
%Wait for the select button to be hit
uiwait()
filesOfInterestUpdated = filesOfInterest(listCheck.Value);
close(listCheckFig)
end

function resumefig(~, ~)
  uiresume()
end