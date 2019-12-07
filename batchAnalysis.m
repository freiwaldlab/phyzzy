function batchAnalysis()
%% Sliding Window ANOVA
% performs an ANOVA on rates across a trial (Fixation, presentation, and
% reward), starting at trial time 0, in some predefined steps.
dataDir = 'D:\DataAnalysis\ANOVA_FullTime';

binSize = 100;
binStep = 25;
figData.binSize = binSize;
figData.binStep = binStep;

%% Get the relevant data
% Method 1: Cycle through analyzedData.mat files in a specified directory,
% store and organize the relevant structures. 
target = {'socialInteraction','agents','interaction'};

preprocessedList = dir([dataDir '\**\preprocessedData.mat']);
analyzedList = dir([dataDir '\**\analyzedData.mat']);
assert(length(preprocessedList) == length(analyzedList), 'Lists arent same length, confirm every preprocessed file has an analyzed file')
trueCellComp = cell(length(preprocessedList),1);

spikeDataBank = struct();
allStimuliVec = {};
for ii = 1:length(preprocessedList)
  tmp = load([preprocessedList(ii).folder filesep preprocessedList(ii).name],'spikesByEvent','eventIDs','eventCategories','preAlign','postAlign');
  tmp2 = load([analyzedList(ii).folder filesep analyzedList(ii).name], 'dateSubject', 'runNum', 'groupLabelsByImage');

  sessField = sprintf('S%s%s', tmp2.dateSubject, tmp2.runNum);
  
  trueCellComp{ii} = [tmp2.dateSubject tmp2.runNum];
  spikeDataBank.(sessField).dateSubject = tmp2.dateSubject;
  spikeDataBank.(sessField).runNum = tmp2.runNum;
  spikeDataBank.(sessField).spikesByEvent = tmp.spikesByEvent;
  spikeDataBank.(sessField).eventIDs = tmp.eventIDs;
  spikeDataBank.(sessField).eventCategories = tmp.eventCategories;
  spikeDataBank.(sessField).groupLabelsByImage = tmp2.groupLabelsByImage;
  spikeDataBank.(sessField).start = -tmp.preAlign;
  spikeDataBank.(sessField).end = tmp.postAlign;
  spikeDataBank.(sessField).groupLabel = target;
  spikeDataBank.(sessField).figDir = preprocessedList(ii).folder;
  allStimuliVec = [allStimuliVec; tmp.eventIDs];
end
allStimuliVec = unique(allStimuliVec);

%% New Stimuli Array
runList = fields(spikeDataBank);
stimLogicalArray = zeros(length(allStimuliVec),length(runstimLogicalArrayList));

for run_ind = 1:length(runList)  
  stimLogicalArray(:,run_ind) = ismember(allStimuliVec,spikeDataBank.(runList{run_ind}).eventIDs);
end

csStimLogicalArray = cumsum(stimLogicalArray,2);
csStimLogicalArray(~stimLogicalArray) = 0; %Final product is matrix which gives 0 for non present stim, count of stim presentation otherwise.

% First stimulus presentation
firstStimPresInd = zeros(length(allStimuliVec),1);
for stim_ind = 1:length(allStimuliVec)
  firstStimPresInd(stim_ind) = find(stimLogicalArray(stim_ind,:),1,'first');
end

%Add this vector to each field, also construct a single dateTime vector.
allDateTimeVec = [];
for run_ind = 1:length(runList)
  spikeDataBank.(runList{run_ind}).stimPresArray = csStimLogicalArray(:,run_ind);
  spikeDataBank.(runList{run_ind}).dateTime = datetime(extractBetween(spikeDataBank.(runList{run_ind}).dateSubject,1,8),'InputFormat','yyyyMMdd');
  if run_ind == 1
    spikeDataBank.(runList{run_ind}).daysSinceLastRec = 1000;
  else
    spikeDataBank.(runList{run_ind}).daysSinceLastRec = days(diff([spikeDataBank.(runList{run_ind-1}).dateTime, spikeDataBank.(runList{run_ind}).dateTime]));
  end
  allDateTimeVec = [allDateTimeVec, spikeDataBank.(runList{run_ind}).dateTime];
end

%Distance between showings
daysSinceLastPres = zeros(size(stimLogicalArray));
for stim_ind = 1:size(stimLogicalArray,1)
  presentationInd = logical(stimLogicalArray(stim_ind,:)); %When was the stim shown
  daysSinceLastPres(stim_ind,presentationInd) = days([1000, diff(allDateTimeVec(presentationInd))]); %Duration between those dates in days
end

%Return this to spikeDataBank, arranged in way that matches stim table.
daysSinceLastRec = [1000, days(diff(allDateTimeVec))];
for run_ind = 1:size(stimLogicalArray,2)
  [~, big2SmallInd] = ismember(spikeDataBank.(runList{run_ind}).eventIDs,allStimuliVec);
  spikeDataBank.(runList{run_ind}).stimPresCount = daysSinceLastPres(big2SmallInd,run_ind);
end

%% Exclude phase 2 recording data, making sure to remove either whole fields
% for a session or individual channel content when appropriate. 
batchRunxls = [dataDir filesep 'BatchRunResults.xlsx'];
recordingLogxls = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Data 2018\RecordingsMoUpdated.xlsx';
[trueCellInd, trueCellRun] = trueCellCount(batchRunxls, recordingLogxls);

for run_ind = 1:length(trueCellComp)
  validInd = trueCellInd(strcmp(trueCellComp{run_ind},trueCellRun));
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
      savePath = [dataDir filesep spikeDataBank.(runList{run_ind}).dateSubject filesep 'Basic' filesep spikeDataBank.(runList{run_ind}).runNum filesep];
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