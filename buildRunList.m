function [dirList] = buildRunList(dataDir, tag)
% Constructs the runList variable for the processRunBatch function based on
% the selection of a directory and a tag, which both denotes the directory
% contents of interest but also changes output formating.

filesOfInterest = dir([dataDir filesep '**' filesep '*.' tag]);

%Allows for customization of list
filesOfInterest = uiDirChoose(filesOfInterest);

 switch tag
  case 'nev'
    %Find the name in the filename.
    [~,B,~] = fileparts(filesOfInterest(1).folder);
    nameInd = regexp(B,'\D');
    %Extract folder names and run names.
    fileNames = {filesOfInterest.name}';
    dirNames = extractBetween(fileNames,1,max(nameInd));
    runNames = extractBetween(fileNames,max(nameInd)+1,'.');
    tmpDirName = unique(dirNames);
    %Cycle through those folders
    for ii = 1:length(tmpDirName)
      runInd = strcmp(dirNames, tmpDirName{ii});
      runListBlock = {tmpDirName{ii}, {runNames{runInd}}};
      dirList_tmp{ii} = runListBlock;
    end
  case 'fig'
    dirList = unique({filesOfInterest.folder}');
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
selectButton.Callback = @(src,event)uiresume();
%Wait for the select button to be hit
uiwait(listCheckFig, 30)
try
  filesOfInterestUpdated = filesOfInterest(listCheck.Value);
  close(listCheckFig)
catch
  error('Nothing selected');
end
end