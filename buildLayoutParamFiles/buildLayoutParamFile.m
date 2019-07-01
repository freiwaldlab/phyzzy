function [layoutParamPath] = buildLayoutParamFile(varargin)    

%Use the config file to collect tags to look for. Each cell is a page. the
%'unit' page will be expanded to account for number of units.

%Directories
outputDir = 'analysisSummary';
if isempty(varargin)
  figDir = 'C:\Data 2018\Analysis_2019\180628Mo\Basic\003';
elseif length(varargin) == 1
  figDir = varargin{1};
end
%figDir = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\plotzzy\TestFigDir';

%subplotFromFigs parameters
params.interRowSpace = 10;
params.interColSpace = 10;
params.verticalMargin = 10;
params.horizontalMargin = 10;
params.pageColor = [1 1 1];
params.equalColWidths = 1;
params.equalRowHeights = 1;

%createSummaryDoc Parameters
params.outputType = 'pdf';
params.maxFigs = 2;

% Page Switches determine what pages are drawn.
pageSwitch.behavior = 1;
pageSwitch.ephys = 1;
pageSwitch.unit = 1;
pageSwitch.LFP = 1;
pageSwitch.MUA = 1;

%% Page specific Parameters - these dictate features particular page switches produce.
pageConfig.behavior.tags = {'Behav','EyeCorr_'};
pageConfig.behavior.maxFigs = 4;

pageConfig.ephys.tags = {'spike','wave'};
pageConfig.ephys.maxFigs = 4;

pageConfig.unit.tags = {'unit'};
pageConfig.unit.maxFigs = 4;

pageConfig.MUA.tags = {'MUA'};
pageConfig.MUA.maxFigs = 4;

pageConfig.LFP.tags = {'LFP'};
pageConfig.MUA.maxFigs = 4;

%Page Sizes
pageSize = [8.5 11];

%The fractions of the above page sizes each picture makes up, W*H.
layoutSize = {...
  {[1 1]},...                             % 1 figure
  {[1 .5],[1 .5]}...                      % 2 figures
  {[1 .5],[.5 .5],[.5 .5]}...             % 3 figures
  {[.5 .5],[.5 .5],[.5 .5],[.5 .5]}...    % 4 figures
  };

layoutPosition = {...
  {[1 1]},...                             % 1 figure
  {[1 1],[2 1]}...                        % 2 figures
  {[1 1],[2 1],[2 2]}...                  % 3 figures
  {[1 1],[1 2],[2 1],[2 2]}...            % 4 figures
  };

%% Setup
%Process the tag to create the appropriate directory
[A, B, ~] = fileparts(figDir);
pathSplit = split(A, ["/","\"]);
phyzzyTag = [pathSplit{end} '_' pathSplit{end-1} B];

finalDir = [outputDir filesep phyzzyTag];

%Make the tmp and output directory if they don't exist
if ~exist(finalDir, 'dir')
  mkdir(finalDir)
end

%% Save
layoutParamFilename = 'layoutParams.mat';
layoutParamPath = [outputDir filesep layoutParamFilename];
save(layoutParamPath);
end