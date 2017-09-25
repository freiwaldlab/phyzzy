function [ ] = buildStimParamFileFOB_WZ( )
%
%
logFilename = '/Users/stephenserene/Desktop/Freiwald/YOGI_DATA/Y170828/Y1708280008.log';
outFilename = '/Users/stephenserene/Desktop/Freiwald/AnalysisSerene/StimParamsFullFOB_WZ.txt';
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

categories = {'facem','apple','btfly','mushr','fribbl','facem01','facem02','facem03','facem04','facem05','facem06','facem07','facem08',...
  'left90','left45','front00','right45','right90','neutral','threat','feargrin','lipsmack'};
categoryLabels = {'face','apple','btfly','mushr','fribbl','facem01','facem02','facem03','facem04','facem05','facem06','facem07','facem08',...
  'left90','left45','front00','right45','right90','neutral','threat','feargrin','lipsmack'};
wildcardCategories = {};  %'Humanhead*_1.'
wildcardCategoryLabels = {};


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
  if isempty(regexp(paramString,'face','ONCE')) && isempty(regexp(paramString,'grayblank.bmp','ONCE')) && ...
      isempty(regexp(paramString,'outlinem.bmp','ONCE')) && isempty(regexp(paramString,'pinknoise.bmp','ONCE'))
    paramString = strcat(paramString,'\t','nonface');
    paramCell = horzcat(paramCell, 'nonface');
  end
  if ~isempty(regexp(paramString,'mushr','ONCE')) || ~isempty(regexp(paramString,'btfly','ONCE')) || ~isempty(regexp(paramString,'apple','ONCE'))
    paramString = strcat(paramString,'\t','object');
    paramCell = horzcat(paramCell, 'object');
  end
  if ~isempty(regexp(paramString,'fribbl','ONCE'))
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
save('StimParamsFullFOB3_WZ.mat','paramArray','categoryLabels','pictureLabels');
end