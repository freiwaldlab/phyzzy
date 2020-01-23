function [plotMat, briefStimList, params]= plotIndex(stimuliList, params)
%% plotIndex
% use to generate an index for retrieving relevant data structures for
% plotting
% Input
% - stimuliList - a cell array of stimuli file names, should line up with a
%   data structure one wishes to index into.
% - plotIndexParam - a struct with the following fields
%   - ().stimParamFile - a full path to a phyzzy style stimParamFile.
%   - ().plotLabels - a cell array list of labels for inclusion. each list
%   should be contained in a cell array. the resulting indicies in plotMat
%   will match the index in the list. 
% Output
% - plotMat - a stimuli * length(plotLabel) matrix with 1s for stimuli to
% include 0, for others. 

% Load broad label correspondences for plotting or for reassigning labels.
tmp = load(params.stimParamsFilename);
totalEventIDs = cell(length(tmp.paramArray),1);
for event_i = 1:length(tmp.paramArray)
  totalEventIDs{event_i} = tmp.paramArray{event_i}{1}; %per the stimParamFile spec, this is the event ID
end
[~, paramSortVec] = ismember(stimuliList, totalEventIDs);
paramArray = tmp.paramArray(paramSortVec);
briefStimList = tmp.pictureLabels(paramSortVec);

% Iterate across stimuli and identify appropriate indicies
plotMat = zeros(length(stimuliList),length(params.plotLabels));

for stim_ind = 1:size(plotMat,1)
  paramStimSet = paramArray{stim_ind};
  for label_ind = 1:size(plotMat,2)
    % If single label, do simple comparison, store logical value
    if ~iscell(params.plotLabels{label_ind})
      plotMat(stim_ind, label_ind) = any(ismember(paramStimSet,params.plotLabels{label_ind}));
    else
      % If Cell, iterate through values, store the max.
      tmp = zeros(length(params.plotLabels{label_ind}),1);
      for subLabel_ind = 1:length(params.plotLabels{label_ind})
        tmp(subLabel_ind) = any(ismember(paramStimSet,params.plotLabels{label_ind}{subLabel_ind}))*subLabel_ind;
      end
      plotMat(stim_ind, label_ind) = max(tmp);
    end
  end
end

% Remove catagories containing 0.
keepInd = sum(plotMat,1) > 0;
plotMat = plotMat(:,keepInd);
params.plotLabels = params.plotLabels(keepInd);


end