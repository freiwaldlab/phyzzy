function batchAnalysis()
%% Sliding Window ANOVA
% performs an ANOVA on rates across a trial (Fixation, presentation, and
% reward), starting at trial time 0, in some predefined steps.
dataDir = 'D:\DataAnalysis\ANOVA_FullTime';
stimParamFile = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy\stimParamFileLib\StimParamFileSocialVids_Full.mat';
spikeDataBaseFilename = 'spikeDataBase.mat';
binSize = 100;
binStep = 25;
target = {'socialInteraction','agents','interaction'}; %Labels which must exist in the stimParamFile associated with the runs. 
batchRunxls = [dataDir filesep 'BatchRunResults.xlsx']; %Batch analysis xlsx produced by processRunBatch.
recordingLogxls = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Data 2018\RecordingsMoUpdated.xlsx'; %Used to exclude phase 2 to give accurate unit counts.
excludePhase2 = 0; % a switch which can be used to remove data from the same neuron collected in subsequent runs. Good for getting accurate counts.
%% Get the relevant data
spikeDataBaseFile = [dataDir filesep spikeDataBaseFilename];
if exist(spikeDataBaseFile)
  load(spikeDataBaseFile,'spikeDataBank','allStimuliVec')
else
  % Cycle through analyzedData.mat files, store and organize the relevant structures.
  preprocessedList = dir([dataDir filesep '**' filesep 'preprocessedData.mat']);
  analyzedList = dir([dataDir filesep '**' filesep 'analyzedData.mat']);
  assert(length(preprocessedList) == length(analyzedList), 'Lists arent same length, confirm every preprocessed file has an analyzed file')
  
  sessionList = cell(length(preprocessedList),1);
  spikeDataBank = struct();
  allStimuliVec = {};
  
  for ii = 1:length(preprocessedList)
    tmp = load([preprocessedList(ii).folder filesep preprocessedList(ii).name],'spikesByEvent','eventIDs','eventCategories','preAlign','postAlign');
    tmp2 = load([analyzedList(ii).folder filesep analyzedList(ii).name], 'dateSubject', 'runNum', 'groupLabelsByImage','psthByImage');
    
    sessField = sprintf('S%s%s', tmp2.dateSubject, tmp2.runNum);
    allStimuliVec = [allStimuliVec; tmp.eventIDs];
    sessionList{ii} = [tmp2.dateSubject tmp2.runNum];
    spikeDataBank.(sessField).dateSubject = tmp2.dateSubject;
    spikeDataBank.(sessField).runNum = tmp2.runNum;
    spikeDataBank.(sessField).spikesByEvent = tmp.spikesByEvent;
    spikeDataBank.(sessField).psthByImage = tmp2.psthByImage;
    spikeDataBank.(sessField).eventIDs = tmp.eventIDs;
    spikeDataBank.(sessField).eventCategories = tmp.eventCategories;
    spikeDataBank.(sessField).groupLabelsByImage = tmp2.groupLabelsByImage;
    spikeDataBank.(sessField).start = -tmp.preAlign;
    spikeDataBank.(sessField).end = tmp.postAlign;
    spikeDataBank.(sessField).groupLabel = target;
    spikeDataBank.(sessField).figDir = preprocessedList(ii).folder;
  end
  save(spikeDataBaseFile, 'spikeDataBank','allStimuliVec')
end
%% stimuli information
%Code below creates a single large vector of stimuli used, and uses this to
%create individual vectors containing which viewing of the stimulus this
%represent (i.e. 'this run represents the 10th viewing of X.avi')
allStimuliVec = unique(allStimuliVec);
runList = fields(spikeDataBank);
stimLogicalArray = zeros(length(allStimuliVec),length(runList));

for run_ind = 1:length(runList)  
  stimLogicalArray(:,run_ind) = ismember(allStimuliVec,spikeDataBank.(runList{run_ind}).eventIDs);
end

%Final product is matrix which gives 0 for non present stim, count of stim presentation otherwise.
csStimLogicalArray = cumsum(stimLogicalArray,2);
csStimLogicalArray(~stimLogicalArray) = 0; 
totalStimPresCount = sum(stimLogicalArray,2);

% When was a stimulus first seen? Index of runList where first presentation took place. 
firstStimPresInd = zeros(length(allStimuliVec),1);
for stim_ind = 1:length(allStimuliVec)
  firstStimPresInd(stim_ind) = find(stimLogicalArray(stim_ind,:),1,'first');
end

%Add this vector to each field, also construct a single dateTime vector.
allDateTimeVec = NaT(1, length(runList));
for run_ind = 1:length(runList)
  spikeDataBank.(runList{run_ind}).stimPresArray = csStimLogicalArray(:,run_ind);
  spikeDataBank.(runList{run_ind}).dateTime = datetime(extractBetween(spikeDataBank.(runList{run_ind}).dateSubject,1,8),'InputFormat','yyyyMMdd');
  if run_ind == 1
    spikeDataBank.(runList{run_ind}).daysSinceLastRec = 1000;
  else
    spikeDataBank.(runList{run_ind}).daysSinceLastRec = days(diff([spikeDataBank.(runList{run_ind-1}).dateTime, spikeDataBank.(runList{run_ind}).dateTime]));
  end
  allDateTimeVec(run_ind) = spikeDataBank.(runList{run_ind}).dateTime;
end
recordingDays = unique(allDateTimeVec');

%Distance between showings
daysSinceLastPres = zeros(size(stimLogicalArray));
for stim_ind = 1:size(stimLogicalArray,1)
  presentationInd = logical(stimLogicalArray(stim_ind,:)); %When was the stim shown
  daysSinceLastPres(stim_ind,presentationInd) = days([1000, diff(allDateTimeVec(presentationInd))]); %Duration between those dates in days
end

%Return this to spikeDataBank, arranged in way that matches stim table.
daysSinceLastRec = [1000, days(diff(allDateTimeVec))];
small2BigInd = zeros(size(stimLogicalArray)); %Used latter for PSTHs.
for run_ind = 1:size(stimLogicalArray,2)
  [~, big2SmallInd] = ismember(spikeDataBank.(runList{run_ind}).eventIDs,allStimuliVec);
  [~, small2BigInd(:,run_ind)] = ismember(allStimuliVec, spikeDataBank.(runList{run_ind}).eventIDs);
  spikeDataBank.(runList{run_ind}).stimPresCount = daysSinceLastPres(big2SmallInd,run_ind);
end

%% Exclude phase 2 recording data, making sure to remove either whole fields
if excludePhase2
  % for a session or individual channel content when appropriate.
  [trueCellInd, trueCellRun] = trueCellCount(batchRunxls, recordingLogxls);
  
  for run_ind = 1:length(sessionList)
    validInd = trueCellInd(strcmp(sessionList{run_ind},trueCellRun));
    if sum(validInd) == 0
      spikeDataBank = rmfield(spikeDataBank,(runList{run_ind}));
    else
      for event_ind = 1:length(spikeDataBank.(runList{run_ind}).spikesByEvent)
        % for every run, use the corresponding 'validInd' to remove non-valid
        % channel data from each spikesByEvent cell.
        spikeDataBank.(runList{run_ind}).spikesByEvent{event_ind} = spikeDataBank.(runList{run_ind}).spikesByEvent{event_ind}(validInd);
      end
    end
  end
  runList = fields(spikeDataBank);
  %Future Note - other variables related to stimulus must also be
  %reconstructed if this happens. Maybe move to before that step.
end
%% Combine PSTH across all runs for a particular stimulus.
% This will crash if the PSTHs aren't the same length.
plotTopStim = 1;
broadLabel = 1; %Transitions individual stimuli to broad catagory (e.g. chasing).
topStimPresThreshold = 100; %At least this many stim presentations to be plotted when plotTopStim is on.
figureTitle = 'per Stimuli';

%Label swapping before plotting
if broadLabel
  figureTitle = 'Broad Labels';
  tmp = load(stimParamFile);
  for event_i = 1:length(tmp.paramArray)
    totalEventIDs{event_i} = tmp.paramArray{event_i}{1}; %per the stimParamFile spec, this is the event ID
  end
  totalEventIDs = totalEventIDs';
  [~, paramSortVec] = ismember(allStimuliVec, totalEventIDs);
  paramArray = tmp.paramArray(paramSortVec);
  allStimuliVec = labelSwap(allStimuliVec, paramArray);
end
[uniqueStimLabels, ~, C] = unique(allStimuliVec);

c = arrayfun(@(x)length(find(C == x)), unique(C), 'Uniform', false);
uniqueStimLabelCounts = cell2mat(c);

%concatonating PSTHs.
allStimuliPSTH = zeros(length(allStimuliVec),size(spikeDataBank.(runList{1}).psthByImage{1}{end},2));
catagoryPSTHBins = cell(length(uniqueStimLabels),1);
for stim_ind = 1:length(allStimuliVec)
  stimIndex = small2BigInd(stim_ind,:);
  psthStimIndex = nonzeros(stimIndex);
  subRunList = runList(find(stimIndex));
  [cumulativeUnsortedPSTH, cumulativeUnitPSTH, cumulativeMUAPSTH] = deal([]);
  for subRun_ind = 1:length(subRunList)
    tmpPSTH = spikeDataBank.(subRunList{subRun_ind}).psthByImage;
    for chan_ind = 1:length(tmpPSTH)
      for unit_ind = 1:length(tmpPSTH{chan_ind})
        if unit_ind == length(tmpPSTH{chan_ind})
          cumulativeMUAPSTH = [cumulativeMUAPSTH; tmpPSTH{chan_ind}{unit_ind}(psthStimIndex(subRun_ind),:)];
        elseif unit_ind == 1
          cumulativeUnsortedPSTH = [cumulativeUnsortedPSTH; tmpPSTH{chan_ind}{unit_ind}(psthStimIndex(subRun_ind),:)];
        else
          cumulativeUnitPSTH = [cumulativeUnitPSTH; tmpPSTH{chan_ind}{unit_ind}(psthStimIndex(subRun_ind),:)];
        end
      end
    end
  end
  cumulativePSTHAll = {cumulativeUnsortedPSTH, cumulativeUnitPSTH, cumulativeMUAPSTH};
  if isempty(catagoryPSTHBins{C(stim_ind)})
    catagoryPSTHBins{C(stim_ind)} = [catagoryPSTHBins{C(stim_ind)}; cumulativePSTHAll];
  else
    for group_ind = 1:3
      catagoryPSTHBins{C(stim_ind)}{group_ind} = [catagoryPSTHBins{C(stim_ind)}{group_ind}; cumulativePSTHAll{group_ind}];
    end
  end
end

%Remove stimuli which aren't displayed over N times (defined at top of
%section.
if plotTopStim && ~broadLabel
  figureTitle = sprintf('stimuli with over %d runs', topStimPresThreshold);
  topIndex = totalStimPresCount > topStimPresThreshold;
  catagoryPSTHBins = catagoryPSTHBins(topIndex,:);
  uniqueStimLabels = uniqueStimLabels(topIndex);
end

%Remove low frequency examples
if broadLabel
  catagoryLowInd = uniqueStimLabelCounts < 4;
  catagoryPSTHBins(catagoryLowInd) = [];
  uniqueStimLabels(catagoryLowInd) = [];
end

groupingType = {'Unsorted','Units','MUA'};
for group_ind = 1:length(groupingType)
  figure()
  hold on
  maxVal = 0;
  minVal = 100;
  for psth_ind = 1:size(catagoryPSTHBins,1)
    meanPSTH = sum(catagoryPSTHBins{psth_ind}{group_ind})/size(catagoryPSTHBins{psth_ind}{group_ind},1);
    plot(meanPSTH);
    maxVal = ceil(max([maxVal, max(meanPSTH)]));
    minVal = floor(min([minVal, min(meanPSTH)]));
  end
  xlim([0,size(spikeDataBank.(runList{1}).psthByImage{1}{end},2)])
  %ylim([20, ceil(max(max(catagoryPSTHyBins{1}{1})))])
  legend(uniqueStimLabels,'Location','northeastoutside');
  plot([800, 800],[minVal,maxVal],'LineWidth',3,'color','k','HandleVisibility','off')
  plot([3600, 3600],[minVal,maxVal],'LineWidth',3,'color','k','HandleVisibility','off')
  title(sprintf('Mean PSTH Across all %s - %s', groupingType{group_ind}, figureTitle))
end
%% generate bin times and spike rates, and proper memberships to groups.
for run_ind = 1:length(runList)
  runStruct = spikeDataBank.(runList{run_ind});
  starts = (runStruct.start:binStep:(runStruct.end - binSize))';
  ends = (runStruct.start+binSize:binStep:(runStruct.end))';
  spikeDataBank.(runList{run_ind}).epochs = [starts,ends];
  spikeDataBank.(runList{run_ind}).epochRates = cell(length(starts),1);
  for bin_ind = 1:length(starts)
    [spikeDataBank.(runList{run_ind}).epochRates{bin_ind}, ~, ~] = spikeCounter(spikeDataBank.(runList{run_ind}).spikesByEvent, starts(bin_ind), ends(bin_ind));
  end
  for group_ind = 1:length(target)
    groupDef = target{group_ind};
    nestStrCmp = @(x) any(strcmp(x, groupDef));
    spikeDataBank.(runList{run_ind}).catagoryInd(:,group_ind) = cell2mat(cellfun(nestStrCmp, runStruct.eventCategories,'UniformOutput',0));
  end
end

%% Perform ANOVA, store values.
save('postSpikeCounter.mat');
plotCells = 0;
depth = 2;

for run_ind = 1:length(runList)
  [pVec, nullVec, errVec, omegaVec, nullO95Vec, nullO99Vec] = deal(recurCellZeros(spikeDataBank.(runList{run_ind}).epochRates{1}, length(spikeDataBank.(runList{run_ind}).epochs), length(target), depth)); % Initialize each p value array.
  [maxOmegaVec] = deal(recurCellZeros(spikeDataBank.(runList{run_ind}).epochRates{1}, 1, length(target), depth)); % Initialize each p value array.
  pStatsVec = recurDouble2ANOVATable(pVec);
  for bin_ind = 1:length(spikeDataBank.(runList{run_ind}).epochs)
    for chan_ind = 1:length(spikeDataBank.(runList{run_ind}).epochRates{bin_ind})
      for unit_ind = 1:length(spikeDataBank.(runList{run_ind}).epochRates{bin_ind}{chan_ind})
        unitResponsePerEvent = spikeDataBank.(runList{run_ind}).epochRates{bin_ind}{chan_ind}{unit_ind};
        catagortyInd = spikeDataBank.(runList{run_ind}).catagoryInd;
        %unitData{event}.rates = trial*1
        for target_ind = 1:length(target)
          [trialSpikes, trialLabels]  = deal([]);
          %grab the relevant events
          targetInd = catagortyInd(:,target_ind);
          targetSpikes = unitResponsePerEvent(targetInd);
          otherSpikes = unitResponsePerEvent(~targetInd);
          %Initialize relevant vecotrs
          spikeGroups = {targetSpikes otherSpikes};
          spikeGroupLabels ={(target{target_ind}) (['non-' target{target_ind}])};
          %Cluster and reshape the arrays properly
          for group_i = 1:length(spikeGroups)
            tmp = spikeGroups{group_i};
            tmp = [tmp{:}];
            dataVec = vertcat(tmp.rates);
            labelVec = repmat(spikeGroupLabels(group_i), length(dataVec),1);
            trialSpikes = vertcat(trialSpikes,dataVec);
            trialLabels = vertcat(trialLabels, labelVec);
          end
          % Check for social v non-social
          [pVec{chan_ind}{unit_ind}(bin_ind,target_ind), pStatsTable, ~] = anovan(trialSpikes,{trialLabels},'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
          top = pStatsTable{2,3} * (pStatsTable{2,5} - pStatsTable{3,5});
          bottom = (pStatsTable{2,3}*pStatsTable{2,5})+(pStatsTable{4,3}-pStatsTable{2,3})*pStatsTable{3,5};
          omegaVec{chan_ind}{unit_ind}(bin_ind,target_ind) = top/bottom;
          [nullPVec, nullOmegaVec] = deal(zeros(1,100));
          for rand_ind = 1:100
            trialLabels = trialLabels(randperm(length(trialLabels)));
            [nullPVec(rand_ind), pStatsTable, ~] = anovan(trialSpikes,{trialLabels},'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
            top = pStatsTable{2,3} * (pStatsTable{2,5} - pStatsTable{3,5});
            bottom = (pStatsTable{2,3}*pStatsTable{2,5})+(pStatsTable{4,3}-pStatsTable{2,3})*pStatsTable{3,5};
            nullOmegaVec(rand_ind) = top/bottom;
          end
          nullVec{chan_ind}{unit_ind}(bin_ind,target_ind) = mean(nullPVec);
          nullO95Vec{chan_ind}{unit_ind}(bin_ind,target_ind) = prctile(nullOmegaVec,95);
          nullO99Vec{chan_ind}{unit_ind}(bin_ind,target_ind) = prctile(nullOmegaVec,99);
          if bin_ind == length(spikeDataBank.(runList{run_ind}).epochs)
            maxOmegaVec{chan_ind}{unit_ind} = max(omegaVec{chan_ind}{unit_ind});
          end
        end
      end
    end
  end
  spikeDataBank.(runList{run_ind}).pVec = pVec;
  spikeDataBank.(runList{run_ind}).nullPVec = nullVec;
  spikeDataBank.(runList{run_ind}).pStatsVec = pStatsVec;
  spikeDataBank.(runList{run_ind}).maxOmegaVec = maxOmegaVec;
  spikeDataBank.(runList{run_ind}).omegaVec = omegaVec;
  spikeDataBank.(runList{run_ind}).nullO95Vec = nullO95Vec;
  spikeDataBank.(runList{run_ind}).nullO99Vec = nullO99Vec;
  %Plot the results for each unit seen
  if plotCells
    for chan_ind = 1:length(pVec)
      for unit_ind = 1:length(pVec{chan_ind})
        if unit_ind == 1
          ANOVAvarName = ['Ch' num2str(chan_ind) ' Unsorted - SocVsNonSoc'];
        elseif unit_ind == length(pVec{chan_ind})
          ANOVAvarName = ['Ch' num2str(chan_ind) ' MUA - SocVsNonSoc'];
        else
          ANOVAvarName = ['Ch' num2str(chan_ind) ' U' num2str(unit_ind-1) ' - SocVsNonSoc'];
        end
        figure()
        plot(pVec{chan_ind}{unit_ind},'linewidth',2)
        hold on
        xlim([0,size(pVec{chan_ind}{unit_ind},1)]);
        ylim([0,1]);
        plot([0,size(pVec{chan_ind}{unit_ind},1)],[0.05, 0.05],'k','linewidth',3)
        title(sprintf('Sliding Scale ANOVA - %s',ANOVAvarName))
        fracSig = round(sum(pVec{chan_ind}{unit_ind} < 0.05)/length(pVec{chan_ind}{unit_ind}), 3, 'significant');
        fragSigNull = round(sum(pVec{chan_ind}{unit_ind} < 0.4)/length(pVec{chan_ind}{unit_ind}), 3, 'significant');
        for leg_ind = 1:length(target)
          legendText{leg_ind} = sprintf('%s (%s;%s)', target{leg_ind}, num2str(fracSig(leg_ind)),num2str(fragSigNull(leg_ind)));
        end
        for null_ind = 1:length(target)
          h.(sprintf('NullLine%d', null_ind)) = mseb(1:length(nullVec{chan_ind}{unit_ind}),nullVec{chan_ind}{unit_ind}(:,null_ind)', errVec{chan_ind}{unit_ind}(:,null_ind)');
          h.(sprintf('NullLine%d', null_ind)).patch.FaceAlpha = '0.5';
        end
        legend(legendText); 
      end
    end
  end
  if mod(run_ind,10) == 0 || run_ind == length(runList)
    save('ANOVAandNull')
  end
end

%% Do Population wide ANOVA
totalUnitCount = 0;
totalChannelCount = 0;
[sigBinCountUnit, sigBinCountUnsorted, sigBinCountMUA, ...
  nulSigBinCountUnit, nulSigBinCountUnsorted, nulSigBinCountMUA] = deal(zeros(length(spikeDataBank.(runList{1}).epochs), length(target)));
sigBinCountTotal = [];
targetRunLengths = struct();
for targ_i = 1:length(target)
  targetRunLengths.(target{targ_i}) = [];
end

for run_ind = 1:length(runList)
  for chan_ind = 1:length(spikeDataBank.(runList{run_ind}).pVec)
    totalChannelCount = totalChannelCount + 1;
    for unit_ind = 1:length(spikeDataBank.(runList{run_ind}).pVec{chan_ind})
      unitLen = length(spikeDataBank.(runList{run_ind}).pVec{chan_ind});
      uPVec = spikeDataBank.(runList{run_ind}).pVec{chan_ind}{unit_ind};
      uNullVec =  spikeDataBank.(runList{run_ind}).nullPVec{chan_ind}{unit_ind};
      %Check for significant bins
      tmpSigCount = uPVec < 0.05;
      tmpNullSigCount = uNullVec < 0.05;
      %Add to relevant structures
      if unit_ind == 1
        sigBinCountUnsorted = tmpSigCount + sigBinCountUnsorted;
        nulSigBinCountUnsorted = tmpNullSigCount + nulSigBinCountUnsorted;
      elseif unit_ind == unitLen
        sigBinCountMUA = tmpSigCount + sigBinCountMUA;
        nulSigBinCountMUA = tmpNullSigCount + nulSigBinCountMUA;
      else
        totalUnitCount = totalUnitCount + 1;
        sigBinCountUnit = tmpSigCount + sigBinCountUnit;
        nulSigBinCountUnit = tmpNullSigCount +nulSigBinCountUnit;
        sigBinCountTotal = [sigBinCountTotal; sum(tmpSigCount)];
        %Keep count of consecutive bins
        for targ_i = 1:size(tmpSigCount,2)
          starts = find(diff([0; tmpSigCount(:,targ_i)]) == 1);
          ends = find(diff([tmpSigCount(:,targ_i); 0]) == -1)+1;
          runLengths = ends - starts;
          targetRunLengths.(target{targ_i}) = [targetRunLengths.(target{targ_i}); runLengths];
        end
      end
    end
  end
end

save('postSlidingANOVA')
%% Plot Population Code
%Figure 1 - Soc vs Non-Soc significance bins
target = {'socialInteraction','agents','interaction'};
allUnitTraces = [sigBinCountUnit, nulSigBinCountUnit];
allMUATraces = [sigBinCountMUA, nulSigBinCountMUA];
allUnsortedTraces = [sigBinCountUnsorted,nulSigBinCountUnsorted];
legendCells = cell(length(target)*2,1);
for legend_i = 1:length(target)
  legendCells{legend_i} = target{legend_i};
  legendCells{legend_i+length(target)} = ['Label scramble ' target{legend_i}];
end
figure
subplot(3,1,1)
plot(allUnitTraces,'Linewidth',3);
legend(legendCells);
ylim([0 totalUnitCount/5]);
xlim([0 length(allUnitTraces)])
title('Significant Units per Bin')
subplot(3,1,2)
plot(allMUATraces,'Linewidth',3);
legend(legendCells);
ylim([0 totalUnitCount/5]);
xlim([0 length(allUnitTraces)])
title('Significant MUA per Bin')
subplot(3,1,3)
plot(allUnsortedTraces,'Linewidth',3);
legend(legendCells);
ylim([0 totalUnitCount/5]);
xlim([0 length(allUnitTraces)])
title('Significant Unsorted per Bin')

% Fig 2 - Numbers of bins per stimuli, length of stretches
figure()
for ii = 1:length(target)
  %Plot significant bin count per unit
  subplot(length(target), 2, (ii*2)-1)
  histogram(sigBinCountTotal(:,ii),50)
  title(sprintf('Number of bins per unit (%s)',target{ii}))
  %Plot stretch duration per unit
  subplot(length(target), 2, (ii * 2))
  histogram(targetRunLengths.(target{ii}))
  title(sprintf('Length of significant bin runs per unit (%s)',target{ii}))
end

%% Plot Omega curves
allMaxOmegas = [];
figData.binSize = binSize; %figData required to save figure.
figData.binStep = binStep;
for run_ind = 1:length(runList)
  omegaVec = spikeDataBank.(runList{run_ind}).omegaVec;
  for chan_ind = 1:length(omegaVec)
    for unit_ind = 1:length(omegaVec{chan_ind})
      if unit_ind == 1
        ANOVAvarName = ['Ch' num2str(chan_ind) ' Unsorted'];
      elseif unit_ind == length(omegaVec{chan_ind})
        ANOVAvarName = ['Ch' num2str(chan_ind) ' MUA'];
      else
        ANOVAvarName = ['Ch' num2str(chan_ind) ' unit ' num2str(unit_ind-1)];
        allMaxOmegas = [allMaxOmegas; spikeDataBank.(runList{run_ind}).maxOmegaVec{chan_ind}{unit_ind}];
      end
      figure()
      plot(omegaVec{chan_ind}{unit_ind},'linewidth',2)
      hold on
      xlim([0,size(omegaVec{chan_ind}{unit_ind},1)]);
      ylim([0,0.2]);
      title(sprintf('Sliding Scale ANOVA Omega - %s, %d ms bin', ANOVAvarName, binSize))
      for leg_ind = 1:length(target)
        legendText{leg_ind} = sprintf('%s', target{leg_ind});
      end
      legend(legendText);
      starts = spikeDataBank.(runList{run_ind}).epochs(:,1);
      xticks(1:10:length(starts))
      xticklabels(starts(1:10:length(starts)));
      fileName = sprintf('Omega Curve - %s, %d ms bins, %d ms step', ANOVAvarName, binSize, binStep);
      savePath = [spikeDataBank.(runList{run_ind}).figDir filesep];
      saveFigure(savePath, fileName, figData, 1, 0, 0, runList{run_ind}, 'close')
    end
  end
end

%Distribution of peak sensitivities 
for ii = 1:size(allMaxOmegas,2)
  figure
  hist(allMaxOmegas(:,ii),20);
  title(target{ii});
end

save('postOmega')

end

function OutCell = recurCellZeros(inCell, N, M, depth)
%Takes a nested cell structure and initializes a blank N * M double array
%within. If the cell structure contains cells, function goes deeper.

%To-do - use this function to initialize all matching size and format
%arrays above.

if iscell(inCell) && depth ~= 0
  OutCell = cell(size(inCell));
  for inCell_ind = 1:length(inCell)
    OutCell{inCell_ind} = recurCellZeros(inCell{inCell_ind}, N, M, depth-1);
  end
else
  OutCell = zeros(N,M);
end
end

function OutCell = recurDouble2ANOVATable(inCell)
%Takes a nested cell structure and initializes a blank N * M double array
%within. If the cell structure contains cells, function goes deeper.

%To-do - use this function to initialize all matching size and format
%arrays above.

if isa(inCell, 'double') && ismatrix(inCell)
  OutCell = cell(size(inCell));
  [OutCell{:}] = deal(cell(4,7));%the size of the ANOVA table...
else
  OutCell = cell(size(inCell));
  for inCell_ind = 1:length(inCell)
    OutCell{inCell_ind} = recurDouble2ANOVATable(inCell{inCell_ind});
  end
end
end

function newLabels = labelSwap(OldLabels, paramArray)
newLabelSet = {'chasing','fighting','mounting','grooming','holding','following','observing',...
  'foraging','sitting','objects','goalDirected','idle','scramble','scene'};

for label_ind = 1:length(OldLabels)
  paramStimSet = paramArray{label_ind};
  newLabels{label_ind} = paramStimSet{ismember(paramStimSet,newLabelSet)};
end
newLabels = newLabels';
end