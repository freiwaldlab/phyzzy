function [ ] = buildStimParamFileSocialVids_Auto()
%Creates a .m file containing cell arrays detailing the parameters of a
%particular stimulus set.

%File produced 2 Cell Arrays. The first is a ParamArray, with the name of
%the video. {{VideoID1; {label1}; {label2}...}{VideoID2; {label1};
%{label2}...}. The second is a TriggerLabels, containing all possible labels. 

% For experiment run 2018.10.10
% Updated 2019.13.1 to include full stim set, including August stuff.

%Find the folder with what you want
fprintf('Select folder containing files...\n');
%StimFolder = uigetdir();
StimFolder = {'G:\StimuliForFaridfromJulia\SocialDecomposed_renamed', 'G:\StimuliForFaridfromJulia\SocialCategories' };


%Find every file in this folder with the extension desired.
%Add other files
for ii = 1:length(StimFolder)
  tmpStimFiles = dir(StimFolder{ii});
  tmpStimFiles = tmpStimFiles(3:end); %Remove the '.' and '..'
  if ii == 1
    tmpStimFilesStack = tmpStimFiles;
  else
    tmpStimFilesStack = [tmpStimFilesStack; tmpStimFiles];
  end
end

tmpStimFiles = tmpStimFilesStack;

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

%% Turn that list into the appropriate file for phyzzy..
%pictureLabels = split(stimList, '.'); %Due to some names being
%.avi_Dephased.avi, split won't work. Going to just chop off the last 4
%characters, assuming these are a 3 character extension + dot.

%pictureLabels = cellfun(@(x) x(1:end-4), stimList, 'UniformOutput', false);
stimListTmp = stimList;
pictureLabels = cell(size(stimList));

for ii = 1:length(stimList)
  modWord = [];
  label = [];
  stim = stimList{ii};
  stimParts = split(stim,["_","."]);
  stimParts = stimParts(1:end-1);
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
  switch code(3)
    case '1';            stimLabels = horzcat(stimLabels ,'chasing');
    case '2';            stimLabels = horzcat(stimLabels ,'fighting');
    case '3';            stimLabels = horzcat(stimLabels ,'mounting');
    case '4';            stimLabels = horzcat(stimLabels ,'grooming');
    case '5';            stimLabels = horzcat(stimLabels ,'holding');
    case '6';            stimLabels = horzcat(stimLabels ,'following');
    case '7';            stimLabels = horzcat(stimLabels ,'observing');
    case '8';            stimLabels = horzcat(stimLabels ,'foraging');
    case '9';            stimLabels = horzcat(stimLabels ,'sitting');
  end
  %Due to issues w/ grouping of controls, adding this below.
  if (strcmp(code(1:2),'11'))
    if strcmp(code(3),'0')
      stimLabels = horzcat(stimLabels ,'goalDirected');
    else
      stimLabels = horzcat(stimLabels ,'socialInteraction');
    end
  end
  if any(strcmp(stimParts, 'Dephased'))
    stimLabels = {'scramble'};
    label = 'Scrambled';
  end
  if (length(stimParts) > 2) && strncmp(stimParts(3),'Decomp',5)
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
  end
  stimLabels = horzcat(stimLabels, 'allStim');
  stimList{ii} = horzcat(stimList{ii}, stimLabels);
  %Make the picture Label
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
    pictureLabels{ii} = [eventTag '_' speciesTag '_' stimParts{2}(end) modWord label];
  else
    pictureLabels{ii} = [stimParts{1} '_' stimParts{2}(end) modWord label];
  end
end

%% Package Outputs
%pictureLabels = pictureLabels(:,1);

categoryLabels = {'agents', 'objects', 'scramble', 'scene', 'interaction', 'nonInteraction', 'idle', 'socialInteraction','nonSocialInteraction'...
  'goalDirected', 'chasing', 'fighting', 'mounting', 'grooming', 'holding', 'following', 'observing',...
  'foraging','sitting','faces','bodies','hands','background','allStim'};
  
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