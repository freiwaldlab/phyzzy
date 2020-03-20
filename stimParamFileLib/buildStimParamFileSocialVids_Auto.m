function [ ] = buildStimParamFileSocialVids_Auto()
%Creates a .m file containing cell arrays detailing the parameters of a
%particular stimulus set.

%File produced 2 Cell Arrays. The first is a ParamArray, with the name of
%the video. {{VideoID1; {label1}; {label2}...}{VideoID2; {label1};
%{label2}...}. The second is a TriggerLabels, containing all possible labels. 

% For experiment run 2018.10.10
% Updated 2019.13.1 to include full stim set, including August stuff.

%Find the folder with what you want
StimFolder = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories';
stimType = '.avi';

%Find every file in this folder with the extension desired.
%Add other files
tmpStimFilesStack = [];
for ii = 1:length(StimFolder)
  tmpStimFiles = dir(fullfile(StimFolder, '**', ['*' stimType]));
  tmpStimFilesStack = [tmpStimFilesStack; tmpStimFiles];
end

% Remove duplicates
stimNames = {tmpStimFilesStack.name};
[~, uniInd] = unique(stimNames);
stimList = {tmpStimFilesStack(uniInd).name}';

%% Turn that list into the appropriate file for phyzzy..
pictureLabels = cell(size(stimList));

for ii = 1:length(stimList)
  stim = stimList{ii};
  modWord = [];
  label = [];
  stimParts = split(stim,["_","."]);
  stimParts = stimParts(1:end-1); %Remove the tag
  code = stimParts{2};
  stimLabels = cell(1);
  switch code(1)
    case '1';            stimLabels{1} = 'agents'; %stimLabels = horzcat(stimLabels, 'faces','bodies','hands','background');
    case '2';            stimLabels{1} = 'objects'; %stimLabels = horzcat(stimLabels, 'goalDirectedOrMoving'); %This only works because I don't use still items.
    case '3';            stimLabels{1} = 'scramble';
    case '4';            stimLabels{1} = 'scene'; %stimLabels{2} = 'nonInteraction';
  end
  switch code(2)
    case '1';            stimLabels = horzcat(stimLabels, 'interaction');
    case '3';            stimLabels = horzcat(stimLabels, 'idle');
  end
  if strcmp(code(1:3), '110')
    stimLabels = horzcat(stimLabels ,'goalDirected');
  end
  % Make sure to not catch controls for animated stimuli
  if (((length(stimParts) > 2) && ~strncmp(stimParts(3),'C',1)) || (length(stimParts) == 2)) && ~strcmp(stimParts{1}, 'scramble')
    switch code(3)
      case '1';            stimLabels = horzcat(stimLabels ,'chasing', 'socialInteraction');
      case '2';            stimLabels = horzcat(stimLabels ,'fighting', 'socialInteraction');
      case '3';            stimLabels = horzcat(stimLabels ,'mounting', 'socialInteraction');
      case '4';            stimLabels = horzcat(stimLabels ,'grooming', 'socialInteraction');
      case '5';            stimLabels = horzcat(stimLabels ,'holding', 'socialInteraction');
      case '6';            stimLabels = horzcat(stimLabels ,'following', 'socialInteraction');
      case '7';            stimLabels = horzcat(stimLabels ,'observing', 'socialInteraction');
      case '8';            stimLabels = horzcat(stimLabels ,'foraging', 'socialInteraction');
      case '9';            stimLabels = horzcat(stimLabels ,'sitting', 'socialInteraction');
    end
  elseif ((length(stimParts) > 2) && strncmp(stimParts(3),'C',1))
    stimLabels = horzcat(stimLabels ,'animControl');
  end
  
  % Animation Processing
  if strcmp(stimParts{2}(end),'A')
    stimLabels = horzcat(stimLabels ,'animated');
    %A change needed for the naming convention on animated stimuli controls
    if length(stimParts) > 2 && strncmp(stimParts(3),'C',1)
      stimLabels = horzcat(stimLabels ,'animControl');
    end
  end
  
  % Decomposed Stimuli processing
  if any(strcmp(stimParts, 'Dephased'))
    stimLabels = {'scramble'};
    label = 'Scrambled';
  end
  
  if (length(stimParts) > 2)
    if strncmp(stimParts(3),'Decomp',5)
      modWord = 'no';
      stimLabels = horzcat(stimLabels, 'Decomp');
      switch stimParts{3}
        case 'DecompA'
          switch stimParts{4}
            case 'B';   stimLabels = horzcat(stimLabels, 'faces','background'); label = 'Bodies';
            case 'F';   stimLabels = horzcat(stimLabels, 'bodies', 'hands','background'); label = 'Faces';
            case 'H';   stimLabels = horzcat(stimLabels, 'bodies', 'faces','background'); label = 'Hands';
            case 'FB';  stimLabels = horzcat(stimLabels, 'background'); label = 'FacesBodies';
            case 'FH';  stimLabels = horzcat(stimLabels, 'bodies','background'); label = 'FaceHands';
          end
        case 'DecompB'
          switch stimParts{4}
            case 'B';   stimLabels = horzcat(stimLabels, 'bodies','hands'); label = 'FacesScene';
            case 'F';   stimLabels = horzcat(stimLabels, 'faces'); label = 'BodiesScene';
            case 'H';   stimLabels = horzcat(stimLabels, 'hands'); label = 'FacesBodiesScene';
            case 'FB';  stimLabels = horzcat(stimLabels, 'faces', 'bodies'); label = 'Scene';
            case 'FH';  stimLabels = horzcat(stimLabels, 'faces', 'hands'); label = 'BodiesScene';
          end
      end
    elseif strncmp(stimParts(3),'C',1)
      label = strjoin(stimParts(3:end));
    end
  end
  
  % Accounts for stimuli with no Labels
  stimLabels = horzcat(stimLabels, 'allStim'); %#ok<*AGROW>
  
  % Generate final structure
  stimList{ii} = horzcat(stimList{ii}, stimLabels);
  
  % Make the picture Label
  if (strcmp(code(1),'1'))
    eventTag = stimParts{1};
    cutStartInd = find(isstrprop(eventTag, 'upper'));
    speciesTag = eventTag(1);
    eventTag = eventTag(cutStartInd: end);
    if strcmp(eventTag(end-2:end), 'ing')
      eventTag = eventTag(1:end-3);
      if strcmp(eventTag, 'Chas')
        eventTag = 'Chase';
      end
    end
    pictureLabels{ii} = [eventTag '_' speciesTag '_' stimParts{2}(4:end) modWord label];
  else
    pictureLabels{ii} = [stimParts{1} '_' stimParts{2}(end) modWord label];
  end
  
end

%% Adding subEvents
% events within Phyzzy may be stimuli or subEvents within the stimuli. to
% facilitate the latter, the code below looks for an 'eventData.mat', and
% add each event within as a 'subEvent'.

eventDataPath = dir([StimFolder '\**\eventData.mat']);
load(fullfile(eventDataPath(1).folder, eventDataPath(1).name), 'eventData');
subEvents2Add = eventData.Properties.VariableNames;
[stimAdd, labelAdd] = deal(cell(length(subEvents2Add),1));
for event_i = 1:length(subEvents2Add)
  event = subEvents2Add(event_i);
  stimLabels = ['subEvents', strsplit(event{1},'_')];
  stimAdd{event_i} = [event, stimLabels];
  % Generate Label
  stimParts = strsplit(event{1},'_');
  if length(stimParts) > 1
    labelAdd{event_i} = [stimParts{1}, '_' upper(stimParts{2}(1))];
  else
    labelAdd{event_i} = stimParts{1};
  end
end
stimList = [stimList; stimAdd];
pictureLabels = [pictureLabels; labelAdd];

%% Package Outputs
%pictureLabels = pictureLabels(:,1);

categoryLabels = {'agents', 'objects', 'scramble', 'scene', 'interaction', 'nonInteraction', 'idle', 'socialInteraction','nonSocialInteraction'...
  'goalDirected', 'chasing', 'fighting', 'mounting', 'grooming', 'holding', 'following', 'observing','animated', 'animControl',...
  'foraging','sitting','faces','bodies','hands','background', 'subEvents', 'headTurn', 'bodyTurn', 'allStim'};
  
paramArray = stimList;

%% Remove replicates
eventIDs = cell(size(paramArray,1),1);
for event_i = 1:length(paramArray)
  eventIDs{event_i} = paramArray{event_i}{1}; %per the stimParamFile spec, this is the event ID
end

[~, uInd, ~] = unique(eventIDs,'stable');

paramArray = paramArray(uInd);
pictureLabels = pictureLabels(uInd);

save('StimParamFileSocialVids_Full.mat','paramArray','categoryLabels','pictureLabels')
end