function [ ] = buildStimParamFileFOB( )
%
%
logFilename = '/Users/stephenserene/Desktop/Freiwald/ALAN_DATA/Visiko/160922ALAN/160922ALAN0001.log';
outFilename = '/Users/stephenserene/Desktop/Freiwald/AnalysisSerene/StimParamsFullFOB3.txt';
outFile = fopen(outFilename, 'w');
logStruct = xml2struct(logFilename);
pictureFilenames = {};
for i = 1:length(logStruct.VISIKOLOG.Object)
  stimulusStruct = logStruct.VISIKOLOG.Object{i};
  if isfield(stimulusStruct,'PictureCompleted')
    continue
  end
  pictureFilenames = vertcat(pictureFilenames, stimulusStruct.Pictures.Picture.Filename.Text);
end
pictureFilenames = unique(pictureFilenames);
pictureLabels = {};
for i = 1:length(pictureFilenames)
  tmp = regexp(pictureFilenames{i}, '\','split');
  tmp2 = regexp(tmp{end},'.bm','split');
  pictureLabels = vertcat(pictureLabels, regexprep(tmp2{1},'\.|_','')); %remove '.' and '_' characters
end
disp(pictureLabels);
categories = {'Humanhead','HumanheadoriA','HumanheadoriB','HumanheadoriC','HumanheadoriD','HumanheadoriE','Monkeyhead',...
  'MonkeyheadoriA','MonkeyheadoriB','MonkeyheadoriC','MonkeyheadoriD','MonkeyheadoriE','fruit','hand','humanbody','monkeybody',...
  'monkeybodypart','monkeybodywhole','place','techno', 'head','headoriA','headoriB','headoriC','headoriD','headoriE'};
categoryLabels = {'humanFace','humanFaceL90','humanFaceL45','humanFaceFront','humanFaceR45','humanFaceR90','monkeyFace',...
  'monkeyFaceL90','monkeyFaceL45','monkeyFaceFront','monkeyFaceR45','monkeyFaceR90','fruit','hand','humanBody','monkeyBody',...
  'monkeyBodyPart','monkeyBodyWhole','place','techno', 'face','faceL90','faceL45','faceFront','faceR45','faceR90'};
wildcardCategories = {'Humanhead*_1.','Humanhead*_07.','Humanhead*_11.','Humanhead*_12.','Humanhead*_14.','Humanhead*_19.','Humanhead*_25.','Humanhead*_26.',...
  'Monkeyhead*_1.','Monkeyhead*_2.','Monkeyhead*_3.','Monkeyhead*_4.','Monkeyhead*_5.','Monkeyhead*_6.','Monkeyhead*_7.','Monkeyhead*_8.'};
wildcardCategoryLabels = {'human01','human07','human11','human12','human14','human19','human25','human26',...
  'monkey01','monkey02','monkey03','monkey04','monkey05','monkey06','monkey07','monkey08'};
% metaCategories = {{'techno','fruit'}};
% metaCategoryLabels = {'objects'};
tuningCurveParamLabels = {'viewing angle (degrees)','viewing angle (degrees)'};
tuningCurveItems = {{'humanFaceL90','humanFaceL45','humanFaceFront','humanFaceR45','humanFaceR90'},...
  {'monkeyFaceL90','monkeyFaceL45','monkeyFaceFront','monkeyFaceR45','monkeyFaceR90'}}; %can be images or categories
tuningCurveItemType = {'category','category'}; % category or image
tuningCurveParamValues = {[-90 -45 0 45 90], [-90 -45 0 45 90]};
tuningCurveTitles = {'Human face view','Monkey face view'};
paramArray = {};
for i = 1:length(pictureFilenames)
  tmp = regexp(pictureFilenames{i},'\','split');
  fname = strtrim(tmp{end});
  paramString = fname;
  paramCell = {fname};
  for j = 1:length(categories)
    if ~isempty(regexp(fname,categories{j}, 'ONCE'))
      paramString = strcat(paramString,'\t',categoryLabels{j});
      paramCell = horzcat(paramCell, categoryLabels{j});
    end  
  end
  for j = 1:length(wildcardCategories)
    if ~isempty(regexp(fname,regexptranslate('wildcard',wildcardCategories{j}), 'ONCE'))
      paramString = strcat(paramString,'\t',wildcardCategoryLabels{j});
      paramCell = horzcat(paramCell, wildcardCategoryLabels{j});
    end 
  end
  if isempty(regexp(paramString,'face','ONCE')) && isempty(regexp(paramString,'grayback.bmp','ONCE'))
    paramString = strcat(paramString,'\t','nonface');
    paramCell = horzcat(paramCell, 'nonface');
  end
  if ~isempty(regexp(paramString,'fruit','ONCE')) || ~isempty(regexp(paramString,'techno','ONCE'))
    paramString = strcat(paramString,'\t','object');
    paramCell = horzcat(paramCell, 'object');
  end
  if ~isempty(regexp(paramString,'hand','ONCE')) || ~isempty(regexp(paramString,'Body','ONCE'))
    paramString = strcat(paramString,'\t','body');
    paramCell = horzcat(paramCell, 'body');
  end
  paramArray = vertcat(paramArray,{paramCell});
  paramString = strcat(paramString,'\n');
  fprintf(outFile,paramString);
end
% get image filenames without path, extension, or underscores
fclose(outFile);
categoryLabels = horzcat(categoryLabels,wildcardCategoryLabels,{'nonface','object','body'});
save('StimParamsFullFOB3.mat','paramArray','categoryLabels','pictureLabels', 'tuningCurveParamLabels','tuningCurveParamValues','tuningCurveItems','tuningCurveItemType','tuningCurveTitles');
end