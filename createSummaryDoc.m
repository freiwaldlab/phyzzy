%% Figure export script
function createSummaryDoc(layoutParamFilename, figDir)
% Function which draws together files found in the output Analysis folder
% into a single document.
%figDir = 'D:\ESIN_Ephys_Files\Analysis\Analyzed\20190930Mo\Basic\001';
%layoutParamFilename = 'buildLayoutParamFile';
%% Parse Arguments
% Want it to accept:
% a directory with .figs or .pngs
% a configuration file which denotes how to layer them, with a default
% config loaded. Default will be One page summary of available data about
% experiment (in my case, MonkeyLogic Fig, waveClus fig, spikeform fig,
% maybe eye correlogram). Following is 1 page per Unit, + Multiunit and LFP
% at the end.

%% Variables
addpath(genpath('dependencies'))
addpath(genpath('buildLayoutParamFiles'))

%% Core
%Load the parameter file
config = feval(layoutParamFilename, figDir);
load(config);

% Turn the Page switches into a cell array.
pageSwitchVal = find(cellfun(@(x)(pageSwitch.(x)),fieldnames(pageSwitch)));
configPages = fields(pageSwitch);
tagsToRecover = configPages(pageSwitchVal);

for page_ind = 1:length(tagsToRecover)
  configPages{page_ind} = eval(sprintf('pageConfig.%s.tags', tagsToRecover{page_ind}));
end

%If 'unit' is in the config file, break it up into the units found in the folder. 
if any(strcmp([configPages{:}],'unit'))
  configPages = unitExpand(configPages, figDir);
end

%Convert cell array of Tags into cell array of figure paths.
docFigPath = cell(length(configPages),1);
truePageInd = 1;
for config_ind = 1:length(configPages)
  figPath = {};
  totFigInd = 1;
  for tag_ind = 1:length(configPages{config_ind})
    figStructTmp = dir([figDir filesep '*' [configPages{config_ind}{tag_ind}] '*']);
    for fig_ind = 1:length(figStructTmp)
      if length(unique(figPath)) >= params.maxFigs
        docFigPath{truePageInd} = unique(figPath);
        clear figPath
        totFigInd = 1;
        truePageInd = truePageInd + 1;
      end
      figPath{totFigInd} = [figStructTmp(fig_ind).folder filesep figStructTmp(fig_ind).name];
      totFigInd = totFigInd + 1;
    end
  end
  docFigPath{truePageInd} = unique(figPath);
  truePageInd = truePageInd + 1;
end

%remove empty cells
nonEmptyInd = cellfun((@(cell) ~isempty(cell)), docFigPath);
docFigPath = docFigPath(nonEmptyInd);
outputDoc = cell(length(docFigPath), 1);

%Cycle through them, Run them into the subplotsFromFigs function.
for page_ind = 1:length(docFigPath)
  figsToPlot = docFigPath{page_ind};
  if ~isempty(figsToPlot)
    
    %Retrieve appropriate sizes for each figure, size them accordingly.
    pageLayout = layoutSize{length(figsToPlot)};
    figHand = gobjects(length(figsToPlot),1);
    for ii = 1:length(figsToPlot)
      figHand(ii) = openfig(figsToPlot{ii});
      figHand(ii).Units = 'inches';
      figHand(ii).Position = [24 24 pageLayout{ii}.*pageSize];
    end
    
    %Retrieve appropriate indicies for each figure, and create the correct
    %cell array for subplotsFromFigs.
    pageIndicies = layoutPosition{length(figsToPlot)};
    figIn = cell(pageIndicies{end}(1),1);
    for ii = 1:length(figHand)
      figInds = pageIndicies{ii};
      figIn{figInds(1)}{figInds(2)} = figHand(ii);
    end
    
    %subplotFigs requires an n*n input. in the case of 3 figures, find rows
    %which don't match this configuration and fix them.
    squareN = max(cellfun((@(cell) length(cell)), figIn));
    for row_ind = 1:length(figIn)
      if length(figIn{row_ind}) ~= squareN
        figIn{row_ind}{end+1:squareN} = figIn{row_ind}{end};   
      end
    end
    
    %run subplotsFromFigs.
    outputDoc{page_ind} = sprintf('%s%d.png', [finalDir, filesep 'Page'], page_ind);
    subplotsFromFigs(finalDir, outputDoc{page_ind}, figIn, params);
    close all %Incase something isn't closed in the function.
    fprintf('Done Page %d of %d...\n', page_ind, length(docFigPath))
  end
end

if strcmp(params.outputType,'pdf')
  %Name of Summary plot
  sumFile = [finalDir filesep phyzzyTag '_Summary.pdf'];
  %Cycle through figures, opening them, turn
  for page_ind = 1:length(outputDoc)
    imagesc(imread(outputDoc{page_ind}))
    currFig = gcf;
    currFig.Units = 'inches';
    currFig.Position = [0 0 8.5 11];
    axis off
    if page_ind == 1
      export_fig(sumFile,'-m1.2', '-pdf','-opengl')
    else
      export_fig(sumFile,'-m1.2', '-pdf','-opengl','-append')
    end
    delete(outputDoc{page_ind});
  end
end
close all
fprintf('Summary Complete.\n')
end

function configPages = unitExpand(configPages, figDir)
%Find all figs with unit in name.
unitFigs = dir([figDir '\*unit*']);
if ~isempty(unitFigs)
  unitFigList = {unitFigs.name}';
  %Extract all the unit numbers
  unitNum = @(str) strsplit(str,{' ','_'});
  unitFigListFrags = cellfun(unitNum, unitFigList, 'UniformOutput', 0);
  unitNumInd = strcmp(unitFigListFrags{1},'1'); %Find the unit number's location
  pullUnitNum = @(str) str2double(str{unitNumInd});
  unitNumList = cellfun(pullUnitNum, unitFigListFrags); %Recover all unit numbers
  unitCount = max(unitNumList); %Get a unit count
  configInset = cell(1, unitCount);
  for unit_ind = 1:unitCount
    configInset{unit_ind} = {sprintf('unit %d', unit_ind)}; %units to place
  end
  configInset{unit_ind+1} = {'unsorted'};
end

%Find the position in the config array they go in.
for config_ind = 1:length(configPages)
  configUnitInd(config_ind) = any(strcmp(configPages{config_ind}, {'unit'}));
end
insetStart = find(configUnitInd);
if ~isempty(unitFigs)
  %Stick in the new cells
  configTmp = cell(1, length(configPages)-1 + length(configInset));
  %before and after Inset
  configTmp(1:insetStart-1) = configPages(1:insetStart-1);
  configTmp(insetStart+length(configInset):end) = configPages(insetStart+1:end);
  configTmp(insetStart:insetStart+length(configInset)-1) = configInset;
  configPages = configTmp;
else
  configPages(insetStart) = [];
end
end