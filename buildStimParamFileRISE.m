function [ ] = buildStimParamFileRISE(logFilename, outFilename )
%
%
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
  pictureLabels = vertcat(pictureLabels, regexprep(tmp2{1},'_',''));
end
categories = {};
categoryLabels = {};
wildcardCategories = {};
wildcardCategoryLabels = {};
% metaCategories = {{'techno','fruit'}};
% metaCategoryLabels = {'objects'};
paramArray = {};
for i = 1:length(pictureFilenames)
  tmp = regexp(pictureFilenames{i},'\','split');
  fname = strtrim(tmp{end});
  paramString = fname;
  paramCell = {fname};
  paramArray = vertcat(paramArray,{paramCell});
  paramString = strcat(paramString,'\n');
  fprintf(outFile,paramString);
end

% RISE tuning curve specification
pictureStems = {};
riseIndsStr = {};
riseIndsNum = [];
for i = 1:length(pictureLabels)
  if ~isempty(regexp(pictureLabels{i},'grey','ONCE'))
    continue
  end
  tmp = regexp(pictureLabels{i},'RISE','split');
  pictureStems = vertcat(pictureStems,tmp(1));
  riseIndsStr = vertcat(riseIndsStr, tmp(end));
  tmp = regexp(tmp{end},'-','split');
  riseIndsNum = vertcat(riseIndsNum, str2double(tmp{1}));
end
basePictures = unique(pictureStems);
riseIndsStr = unique(riseIndsStr);
riseIndsNum = unique(riseIndsNum);

% tuning curve specification
tuningCurveParamLabels = repmat({'Percent Phase Scramble'},length(basePictures),1);

tuningCurveItems = cell(length(basePictures),1);
for basePic_i = 1:length(basePictures)
  basePic = basePictures{basePic_i};
  items = cell(length(riseIndsNum),1);
  for riseInd = 1:length(riseIndsStr)
    items{riseInd} = strcat(basePic,'RISE',riseIndsStr{riseInd});
  end
  tuningCurveItems{basePic_i} = items;
end
tuningCurveItemType = repmat({'image'},length(basePictures),1); % category or image
tuningCurveParamValues = repmat({riseIndsNum*(100/max(riseIndsNum))},length(basePictures),1);
tuningCurveTitles = basePictures;
% end tuning curve specification

fclose(outFile);
save('StimParamsRISE.mat','paramArray','categoryLabels','pictureLabels', 'tuningCurveParamLabels','tuningCurveParamValues','tuningCurveItems','tuningCurveItemType','tuningCurveTitles');
end




