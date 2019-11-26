%% Sliding Window ANOVA
% performs an ANOVA on rates across a trial (Fixation, presentation, and
% reward), starting at trial time 0, in some predefined steps.
dataDir = 'D:\DataAnalysis\ANOVA_FullTime';

startTime = 0;
binSize = 100;
binStep = 25;

%% Get the relevant data
% Method 1: Cycle through analyzedData.mat files in a specified directory,
% store and organize the relevant structures. 
target = {'socialInteraction','agents','interaction'};

preprocessedList = dir([dataDir '\**\preprocessedData.mat']);
analyzedList = dir([dataDir '\**\analyzedData.mat']);
assert(length(preprocessedList) == length(analyzedList), 'Lists arent same length, confirm every preprocessed file has an analyzed file')
spikeDataBank = struct();

for ii = 1:length(preprocessedList)
  tmp = load([preprocessedList(ii).folder filesep preprocessedList(ii).name],'spikesByEvent','eventIDs','eventCategories','preAlign','postAlign');
  tmp2 = load([analyzedList(ii).folder filesep analyzedList(ii).name], 'dateSubject', 'runNum', 'groupLabelsByImage');

  sessField = sprintf('S%s%s', tmp2.dateSubject, tmp2.runNum);
  
  spikeDataBank.(sessField).dateSubject = tmp2.dateSubject;
  spikeDataBank.(sessField).runNum = tmp2.runNum;
  spikeDataBank.(sessField).spikesByEvent = tmp.spikesByEvent;
  spikeDataBank.(sessField).eventIDs = tmp.eventIDs;
  spikeDataBank.(sessField).eventCategories = tmp.eventCategories;
  spikeDataBank.(sessField).groupLabelsByImage = tmp2.groupLabelsByImage;
  spikeDataBank.(sessField).start = -tmp.preAlign;
  spikeDataBank.(sessField).end = tmp.postAlign;
  spikeDataBank.(sessField).groupLabel = target;
end

runList = fields(spikeDataBank);

%% generate bin times and spike rates, and proper memberships to groups.

for run_ind = 1:length(runList)
  runStruct = spikeDataBank.(runList{run_ind});
  starts = (startTime:binStep:(runStruct.end - binSize))';
  ends = (startTime+binSize:binStep:(runStruct.end))';
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
%load('postSpikeCounter.mat');
plotCells = 0;

for run_ind = 1:length(runList)
  %spikeCountsByImageByEpoch{epoch}{channel}{unit}{event}.rates = trials*1
  [pVec, nullVec, errVec] = deal(cell(1,length(spikeDataBank.(runList{run_ind}).epochRates{1})));
  [pVec{:}, nullVec{:}, errVec{:}] = deal(cell(1,length(spikeDataBank.(runList{run_ind}).epochRates{1}{1})));
  [pVec{:}{:}, nullVec{:}{:}, errVec{:}{:}] = deal(zeros(length(spikeDataBank.(runList{run_ind}).epochs), length(target)));
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
          [pVec{chan_ind}{unit_ind}(bin_ind,target_ind), ~, ~] = anovan(trialSpikes,{trialLabels},'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
          nullPVec = zeros(1,100);
          for rand_ind = 1:100
            trialLabels = trialLabels(randperm(length(trialLabels)));
            nullPVec(rand_ind) = anovan(trialSpikes,{trialLabels},'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
            assert(~isnan(nullPVec(rand_ind)),'Nan found');
          end
          nullVec{chan_ind}{unit_ind}(bin_ind,target_ind) = mean(nullPVec);
          errVec{chan_ind}{unit_ind}(bin_ind,target_ind) = std(nullPVec)/length(nullPVec);
        end
      end
    end
    spikeDataBank.(runList{run_ind}).pVec = pVec;
    spikeDataBank.(runList{run_ind}).nullPVec = nullPVec;
    spikeDataBank.(runList{run_ind}).errVec = errVec;
  end
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
end

%% Do Population wide ANOVA
save('postSlidingANOVA')
disp('Population')




