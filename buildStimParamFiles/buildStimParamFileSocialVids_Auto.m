function [ ] = buildStimParamFileSocialVids_Auto()
%Creates a .m file containing cell arrays detailing the parameters of a
%particular stimulus set.

%File produced 2 Cell Arrays. The first is a ParamArray, with the name of
%the video. {{VideoID1; {label1}; {label2}...}{VideoID2; {label1};
%{label2}...}. The second is a TriggerLabels, containing all possible labels. 

% For experiment run 2018.10.10

%Find the folder with what you want
fprintf('Select folder containing files...\n');
%StimFolder = uigetdir();
StimFolder = 'E:\StimuliForFaridfromJulia\SocialDecomposed_renamed';

%Find every file in this folder with the extension desired.
tmpStimFiles = dir(StimFolder);
tmpStimFiles = tmpStimFiles(3:end); %Remove the '.' and '..'

%Anything that is not a directory, add to the list...

stimInd = find([tmpStimFiles.isdir] == 0);
if ~isempty(stimInd)
    stimList = cell(length(stimInd),1);
    for ii = 1:length(stimInd)
        stimList{ii} = tmpStimFiles(stimInd(ii)).name;
    end
else
    stimList = cell(0);
end

%Look through the remaining lists, going as deep as needed.
dirInd = find([tmpStimFiles.isdir] == 1);
for ii = 1:length(dirInd)
  tmpStimFolder = dir([tmpStimFiles(dirInd(ii)).folder filesep tmpStimFiles(dirInd(ii)).name]);
  subStimFolder = struct2cell(tmpStimFolder(3:end));
  subStimList = subStimFolder(1,:)';
  stimList = vertcat(stimList, subStimList);
end

%% Turn that list into the appropriate file for phyzzy.
pictureLabels = stimList; %Nice names won't exist for now;

for ii = 1:length(stimList)
  stim = stimList{ii};
  stimParts = split(stim,["_","."]);
  code = stimParts{2};
  stimLabels = cell(1);
  switch code(1)
    case '1';            stimLabels{1} = 'agents';
    case '2';            stimLabels{1} = 'objects';
    case '3';            stimLabels{1} = 'scrambles';
    case '4';            stimLabels{1} = 'landscapes';
  end
  switch code(2)
    case '0';            stimLabels = horzcat(stimLabels, 'N/A');
    case '1';            stimLabels = horzcat(stimLabels, 'interaction');
    case '2';            stimLabels = horzcat(stimLabels, 'GoalDirected/Moving');
    case '3';            stimLabels = horzcat(stimLabels, 'Still/Idle');
  end
  switch code(3)
    case '0';            stimLabels = horzcat(stimLabels, 'N/A');   
    case '1';            stimLabels = horzcat(stimLabels, 'chasing');
    case '2';            stimLabels = horzcat(stimLabels, 'fighting');
    case '3';            stimLabels = horzcat(stimLabels, 'mounting');
    case '4';            stimLabels = horzcat(stimLabels, 'grooming');
    case '5';            stimLabels = horzcat(stimLabels, 'holding');
    case '6';            stimLabels = horzcat(stimLabels, 'following');
    case '7';            stimLabels = horzcat(stimLabels, 'observing');
    case '8';            stimLabels = horzcat(stimLabels, 'foraging');
    case '9';            stimLabels = horzcat(stimLabels, 'sitting');
  end
  if (length(stimParts) == 3)
      stimLabels = horzcat(stimLabels, 'faces','bodies','hands','background');
  end
  if (length(stimParts) == 4) && strcmp(stimParts{3}, 'Dephased')
      stimLabels = {'scramble'};
  end
  if (length(stimParts) > 3) && strncmpi(stimParts(3),'Decomp',6)
    switch stimParts{3}
        case 'DecompA'
            switch stimParts{4}
                case 'B';   stimLabels = horzcat(stimLabels, 'faces','background');
                case 'F';   stimLabels = horzcat(stimLabels, 'bodies', 'hands','background');
                case 'H';   stimLabels = horzcat(stimLabels, 'bodies', 'faces','background');
                case 'FB';  stimLabels = horzcat(stimLabels, 'background');
                case 'FH';  stimLabels = horzcat(stimLabels, 'bodies','background');
            end
        case 'DecompB'
            switch stimParts{4}
                case 'B';   stimLabels = horzcat(stimLabels, 'bodies','hands');
                case 'F';   stimLabels = horzcat(stimLabels, 'faces');
                case 'H';   stimLabels = horzcat(stimLabels, 'hands');
                case 'FB';  stimLabels = horzcat(stimLabels, 'faces', 'bodies');
                case 'FH';  stimLabels = horzcat(stimLabels, 'faces', 'hands');
            end
    end    
  end
  stimList{ii} = horzcat(stimList{ii}, stimLabels);
end

categoryLabels = {'agents', 'objects', 'scrambles', 'landscapes', 'interaction', 'GoalDirected/Moving','Still/Idle'...
    'chasing', 'fighting', 'mounting', 'grooming', 'holding', 'following', 'observing','foraging','sitting'...
    'faces','bodies','hands','background'};

paramArray = stimList;

% pictureLabels = {}    %labelscheme
%     %loop through stimList, maybe move to loop above
%     %If the name is 3 long (inc .avi), name + 4th code number
%     %If it is 4 long, name + scramble
%     %If it is 5 long, look at DecompA or B, add word "minus" or "plus",
%     %then add labels as single letters.
%     pictureLabels

save('StimParamFileSocialVids_Decomp.mat','paramArray','categoryLabels','pictureLabels')