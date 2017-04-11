function [  ] = runAnalyses( analysisParamFilename, spikesByChannel, lfpData, analogInData, taskData, taskDataAll, psthImDur, preAlign, postAlign, ...
  categoryList, pictureLabels, jumpsByImage, spikesByImage, psthEmptyByImage, spikesByCategory, psthEmptyByCategory,...
  spikesByImageForTF, spikesByCategoryForTF, lfpByImage, lfpByCategory, channelUnitNames, stimTiming)
%runAnalyses should be the main site for customization
%   - ideally, function signature should remain constant
%  this version of runAnalyses does the following:
%   - makes psth's and evoked potentials for (possibly overlapping)categories, 
%     specified at picParamsFilename
%   - RF analyses: spike rate, evoked power, mean spike latency
%   - performs time-frequency and spike-field coherence analyses for each
%     channel, using the multitaper method implemented in Chronux
%   - performs across channel coupling analyses: spike-spike, field-field,
%     and spike-field, using Chronux
%   - makes tuning curves for parameterized images or categories, as
%     specified in the stim param file
%   - todo: perform band-resolved Granger causality analysis across bands
%     - currently, assumes that 'face' and 'nonface' appear as categories in
%     picParamsFilename

load(analysisParamFilename);
pictureLabelsTmp = pictureLabels; %note: hack to avoid overwriting list of not presented stimuli
load(picParamsFilename);
pictureLabels = pictureLabelsTmp; % conclusion of hack
channelNames = ephysParams.channelNames;
spikeChannels = ephysParams.spikeChannels;
lfpChannels = ephysParams.lfpChannels;
psthPre = psthParams.psthPre;
psthPost = psthParams.psthPost;
smoothingWidth = psthParams.smoothingWidth;

lfpPreAlign = lfpAlignParams.msPreAlign; 
lfpPostAlign = lfpAlignParams.msPostAlign;

lfpPaddedBy = tfParams.movingWin(1)/2;

movingWin = tfParams.movingWin;
specgramRowAve = tfParams.specgramRowAve;
samPerMS = ephysParams.samPerMS;
if frCalcOff < frCalcOn
  frCalcOff = psthImDur+frCalcOn;
end
%spike psth color plot
if makeImPSTH
  for channel_i = 1:length(spikesByImage{1})
    for unit_i = 1:length(spikesByImage{1}{channel_i}) 
      if length(spikesByImage{1}{channel_i}) == 2 && unit_i == 1
        continue;
      end
      imagePSTH = zeros(length(pictureLabels),psthPre+1+psthImDur+psthPost);
      for image_i = 1:length(pictureLabels)
        if ~psthEmptyByImage{image_i}{channel_i}{unit_i}            
          paddedPsth = 1000*psth(spikesByImage{image_i}{channel_i}{unit_i},smoothingWidth,'n',[-preAlign postAlign],0,-preAlign:postAlign);
          imagePSTH(image_i,:) = paddedPsth(3*smoothingWidth+1:end-3*smoothingWidth);
        end
      end
      psthTitle = sprintf('%s, %s',channelNames{channel_i}, channelUnitNames{channel_i}{unit_i});
      figure('Name',psthTitle,'NumberTitle','off');
      plotPSTH(imagePSTH, [], psthPre, psthPost, psthImDur, 'color', psthTitle, pictureLabels, psthColormap );
      clear figData
      figData.z = imagePSTH;
      figData.x = -psthPre:psthImDur+psthPost;
      saveFigure(outDir, sprintf('imPSTH_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
end
if makeCatPSTH
  catMuaPsthByChannel = zeros(length(spikesByImage{1}),length(categoryList),psthPre+1+psthImDur+psthPost);
  for channel_i = 1:length(spikesByImage{1})
    for unit_i = 1:length(spikesByImage{1}{channel_i})
      if length(spikesByImage{1}{channel_i}) == 2 && unit_i == 1
        continue;
      end
      categoryPSTH = zeros(length(categoryList),psthPre+1+psthImDur+psthPost);
      for cat_i = 1:length(categoryList)
        if ~psthEmptyByCategory{cat_i}{channel_i}{unit_i}         
          paddedPsth = 1000*psth(spikesByCategory{cat_i}{channel_i}{unit_i},smoothingWidth,'n',[-preAlign postAlign],0,-preAlign:postAlign);
          categoryPSTH(cat_i,:) = paddedPsth(3*smoothingWidth+1:end-3*smoothingWidth);
        end
      end
      psthTitle = sprintf('%s, %s',channelNames{channel_i}, channelUnitNames{channel_i}{unit_i}); 
      figure('Name',psthTitle,'NumberTitle','off');
      plotPSTH(categoryPSTH, [], psthPre, psthPost, psthImDur, 'color', psthTitle, categoryList, psthColormap );
      clear figData
      figData.z = categoryPSTH;
      figData.x = -psthPre:psthImDur+psthPost;
      saveFigure(outDir, sprintf('catPSTH_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
    catMuaPsthByChannel(channel_i,:,:) = categoryPSTH;
  end
end

[imSpikeCounts, imFr, imFrErr] = spikeCounter(spikesByImage, frCalcOn, frCalcOff);
[imSpikeCountsEarly, imFrEarly, imFrErrEarly] = spikeCounter(spikesByImage, frCalcOnEarly, frCalcOffEarly);
[imSpikeCountsLate, imFrLate, imFrErrLate] = spikeCounter(spikesByImage, frCalcOnLate, frCalcOffLate);

[catSpikeCounts, catFr, catFrErr] = spikeCounter(spikesByCategory, frCalcOn, frCalcOff);
[catSpikeCountsEarly, catFrEarly, catFrErrEarly] = spikeCounter(spikesByCategory, frCalcOnEarly, frCalcOffEarly);
[catSpikeCountsLate, catFrLate, catFrErrLate] = spikeCounter(spikesByCategory, frCalcOnLate, frCalcOffLate);

faceCatNum = find(strcmp(categoryList,'face'));
nonfaceCatNum = find(strcmp(categoryList,'nonface'));
bodyCatNum = find(strcmp(categoryList, 'body'));
objectCatNum = find(strcmp(categoryList, 'object'));
categorySlimInds = zeros(length(categoryListSlim),1);
for i = 1:length(categoryListSlim)
  categorySlimInds(i) = find(strcmp(categoryList,categoryListSlim{i}));
end
colors = ['b','c','y','g','m','r','k'];
colorsCell = {'b';'c';'y';'g';'m';'r';'k'};
chColors = ['b','g','m'];

if ~taskData.RFmap
  % preferred images
  for channel_i = 1:length(spikeChannels)
    for unit_i = 1:length(channelUnitNames{channel_i})
      [imageSortedRates, imageSortOrder] = sort(imFr{channel_i}(unit_i,:),2,'descend');
      imFrErrSorted = imFrErr{channel_i}(unit_i,imageSortOrder);
      sortedImageLabels = pictureLabels(imageSortOrder);
      fprintf('\n\n\nPreferred Images: %s, %s\n\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i});
      for i = 1:10
        fprintf('%d) %s: %.2f +/- %.2f Hz\n',i,sortedImageLabels{i},imageSortedRates(i),imFrErrSorted(i));
      end
      fprintf('\nLeast Preferred Images: %s, %s\n\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i});
      for i = 1:5
        fprintf('%d) %s: %.2f +/- %.2f Hz\n',i,sortedImageLabels{end-i},imageSortedRates(end-i), imFrErrSorted(end-i));
      end
      % preferred images raster plot
      figure();
      raster(spikesByImage(imageSortOrder(1:10)), sortedImageLabels(1:10), psthPre, psthPost, psthImDur, stimTiming.ISI, channel_i, unit_i, colors);
      title(sprintf('Preferred Images, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
      saveFigure(outDir, sprintf('prefImRaster_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      % preferred images raster-evoked overlay
      figure();
      rasterEvoked(spikesByImage(imageSortOrder(1:10)), lfpByImage(imageSortOrder(1:10)), sortedImageLabels(1:10), psthPre, psthPost, psthImDur, stimTiming.ISI, lfpPaddedBy, channel_i, colors, 1)
      title(sprintf('Preferred Images, from top, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
      saveFigure(outDir, sprintf('prefImRaster-LFP_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      % preferred images raster-evoked overlay, with other channels
      figure();
      rasterEvoked(spikesByImage(imageSortOrder(1:10)), lfpByImage(imageSortOrder(1:10)), sortedImageLabels(1:10), psthPre, psthPost, psthImDur, stimTiming.ISI, lfpPaddedBy, channel_i, colors, 1)
      title(sprintf('Preferred Images, from top, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
      saveFigure(outDir, sprintf('prefImRaster-LFP_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end


  % multi-channel MUA image preference
  if length(channelNames) > 1
    multiChSpikesMin = zeros(length(pictureLabels));
    multiChSpikesMinNorm = zeros(length(pictureLabels));
    multiChSpikesMeanNorm = zeros(length(pictureLabels));
    multiChMua = zeros(length(spikeChannels),length(pictureLabels));
    multiChMuaNorm = zeros(length(spikeChannels),length(pictureLabels));
    for channel_i = 1:length(spikeChannels)
      multiChMua(channel_i,:) = imFr{channel_i}(end,:);
      multiChMuaNorm(channel_i,:) = imFr{channel_i}(end,:)/max(imFr{channel_i}(end,:));
    end
    multiChSpikesMin = min(multiChMua);
    multiChSpikesMinNorm = min(multiChMuaNorm);
    multiChSpikesMeanNorm = mean(multiChMuaNorm);
    % multi-channel preferred images
    [imageSortedRates, imageSortOrder] = sort(multiChSpikesMin,2,'descend');
    sortedImageLabels = pictureLabels(imageSortOrder);
    fprintf('\n\n\nMulti-channel Preferred Images, Maximin\n\n');
    for i = 1:10
      fprintf('%d) %s: %.2f Hz\n',i,sortedImageLabels{i},imageSortedRates(i));
    end
    fprintf('\n\nMulti-channel Least preferred Images, Maximin\n\n');
    for i = 1:5
      fprintf('%d) %s: %.2f Hz \n',i,sortedImageLabels{end-i},imageSortedRates(end-i));
    end

    [imageSortedRates, imageSortOrder] = sort(multiChSpikesMinNorm,2,'descend');
    sortedImageLabels = pictureLabels(imageSortOrder);
    fprintf('\n\n\nMulti-channel Preferred Images, Channel-Normalized Maximin\n\n');
    for i = 1:10
      fprintf('%d) %s: %.2f Hz\n',i,sortedImageLabels{i},imageSortedRates(i));
    end
    fprintf('\n\nMulti-channel Least preferred Images, Channel-Normalized Maximin\n\n');
    for i = 1:5
      fprintf('%d) %s: %.2f Hz \n',i,sortedImageLabels{end-i},imageSortedRates(end-i));
    end

    [imageSortedRates, imageSortOrder] = sort(multiChSpikesMeanNorm,2,'descend');
    sortedImageLabels = pictureLabels(imageSortOrder);
    fprintf('\n\n\nMulti-channel Preferred Images, Channel-Normalized Mean\n\n');
    for i = 1:10
      fprintf('%d) %s: %.2f Hz\n',i,sortedImageLabels{i},imageSortedRates(i));
    end
    fprintf('\n\nMulti-channel Least preferred Images, Channel-Normalized Mean\n\n');
    for i = 1:5
      fprintf('%d) %s: %.2f Hz \n',i,sortedImageLabels{end-i},imageSortedRates(end-i));
    end
  end
  % face selectivity index
  for channel_i = 1:length(spikeChannels)
    for unit_i = 1:length(channelUnitNames{channel_i})
      faceVsNonIndex = (catFr{channel_i}(unit_i,faceCatNum) - catFr{channel_i}(unit_i,nonfaceCatNum)) / (catFr{channel_i}(unit_i,faceCatNum) + catFr{channel_i}(unit_i,nonfaceCatNum));
      faceVsObjectIndex = (catFr{channel_i}(unit_i,faceCatNum) - catFr{channel_i}(unit_i,objectCatNum)) / (catFr{channel_i}(unit_i,faceCatNum) + catFr{channel_i}(unit_i,objectCatNum));
      faceVsBodyIndex = (catFr{channel_i}(unit_i,faceCatNum) - catFr{channel_i}(unit_i,bodyCatNum)) / (catFr{channel_i}(unit_i,faceCatNum) + catFr{channel_i}(unit_i,bodyCatNum));
      fprintf('\n\n%s %s Preference Indices, %dms - %d ms\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},frCalcOn,frCalcOff)
      fprintf('Face > Nonface: %.3f\n',faceVsNonIndex);
      fprintf('Face > Object: %.3f\n',faceVsObjectIndex);
      fprintf('Face > Body: %.3f\n',faceVsBodyIndex);
    end
  end
  % face selectivity index, early
  for channel_i = 1:length(spikeChannels)
    for unit_i = 1:length(channelUnitNames{channel_i})
      faceVsNonIndex = (catFrEarly{channel_i}(unit_i,faceCatNum) - catFrEarly{channel_i}(unit_i,nonfaceCatNum)) / (catFrEarly{channel_i}(unit_i,faceCatNum) + catFrEarly{channel_i}(unit_i,nonfaceCatNum));
      faceVsObjectIndex = (catFrEarly{channel_i}(unit_i,faceCatNum) - catFrEarly{channel_i}(unit_i,objectCatNum)) / (catFrEarly{channel_i}(unit_i,faceCatNum) + catFrEarly{channel_i}(unit_i,objectCatNum));
      faceVsBodyIndex = (catFrEarly{channel_i}(unit_i,faceCatNum) - catFrEarly{channel_i}(unit_i,bodyCatNum)) / (catFrEarly{channel_i}(unit_i,faceCatNum) + catFrEarly{channel_i}(unit_i,bodyCatNum));
      fprintf('\n\n%s %s Preference Indices, %dms - %d ms\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},frCalcOnEarly,frCalcOffEarly)
      fprintf('Face > Nonface: %.3f\n',faceVsNonIndex);
      fprintf('Face > Object: %.3f\n',faceVsObjectIndex);
      fprintf('Face > Body: %.3f\n',faceVsBodyIndex);
    end
  end
  % face selectivity index, late
  for channel_i = 1:length(spikeChannels)
    for unit_i = 1:length(channelUnitNames{channel_i})
      faceVsNonIndex = (catFrLate{channel_i}(unit_i,faceCatNum) - catFrLate{channel_i}(unit_i,nonfaceCatNum)) / (catFrLate{channel_i}(unit_i,faceCatNum) + catFrLate{channel_i}(unit_i,nonfaceCatNum));
      faceVsObjectIndex = (catFrLate{channel_i}(unit_i,faceCatNum) - catFrLate{channel_i}(unit_i,objectCatNum)) / (catFrLate{channel_i}(unit_i,faceCatNum) + catFrLate{channel_i}(unit_i,objectCatNum));
      faceVsBodyIndex = (catFrLate{channel_i}(unit_i,faceCatNum) - catFrLate{channel_i}(unit_i,bodyCatNum)) / (catFrLate{channel_i}(unit_i,faceCatNum) + catFrLate{channel_i}(unit_i,bodyCatNum));
      fprintf('\n\n%s %s Preference Indices, %dms - %d ms\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},frCalcOnLate,frCalcOffLate)
      fprintf('Face > Nonface: %.3f\n',faceVsNonIndex);
      fprintf('Face > Object: %.3f\n',faceVsObjectIndex);
      fprintf('Face > Body: %.3f\n',faceVsBodyIndex);
    end
  end

  % category preference bar plot
  % full
  for channel_i = 1:length(spikeChannels)
    figure();
    if length(channelUnitNames{channel_i}) == 2
      ax = subplot(1,2,1);
      superbar(ax, catFr{channel_i}(2, categorySlimInds),'E',catFrErr{channel_i}(2, categorySlimInds),'BarFaceColor',colorsCell);
      set(gca,'XTickLabel',categoryList(categorySlimInds),'XTickLabelRotation',45,'XTick',1:length(categorySlimInds),'TickDir','out');
      ylabel('Firing rate, Hz');
      ax = subplot(1,2,2);
      superbar(ax, catFr{channel_i}(2, [faceCatNum objectCatNum bodyCatNum]),'E',catFrErr{channel_i}(2, [faceCatNum objectCatNum bodyCatNum]),'BarFaceColor',{'b'; 'g'; 'r'}); %'P', [0 1 1; 1 0 0; 1 0 0]
      set(gca,'XTickLabel',{'face','object', 'body'},'XTickLabelRotation',45,'XTick',1:3,'TickDir','out');
      ylabel('Firing rate, Hz');
      saveFigure(outDir, sprintf('catPrefBars_%s_MUA_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    else
      for unit_i = 1:length(channelUnitNames{channel_i})
        ax = subplot(1,length(channelUnitNames{channel_i})+1,unit_i);
        superbar(ax,catFr{channel_i}(unit_i, categorySlimInds),'E',catFrErr{channel_i}(2, categorySlimInds),'BarFaceColor',colorsCell);
        set(gca,'XTickLabel',categoryList(categorySlimInds),'XTickLabelRotation',45,'XTick',1:length(categorySlimInds),'TickDir','out');
        ylabel('Firing rate, Hz');
        title(channelUnitNames{channel_i}{unit_i});
      end
      ax = subplot(1,length(channelUnitNames{channel_i})+1,length(channelUnitNames{channel_i})+1);
      superbar(ax, catFr{channel_i}(end, [faceCatNum objectCatNum bodyCatNum]),'E',catFrErr{channel_i}(end, [faceCatNum objectCatNum bodyCatNum]),'BarFaceColor',{'b'; 'g'; 'r'}); %'P', [0 1 1; 1 0 0; 1 0 0]
      set(gca,'XTickLabel',{'face','object','body'}, 'XTickLabelRotation',45,'XTick',1:3,'TickDir','out');
      ylabel('Firing rate, Hz');
      title('MUA');
    end
    suptitle(sprintf('%s spikes category tuning, %dms-%dms post-onset',channelNames{channel_i},frCalcOn,frCalcOff));
    saveFigure(outDir, sprintf('catPrefBars_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  end

  % early
  for channel_i = 1:length(spikeChannels)
    figure();
    if length(channelUnitNames{channel_i}) == 2
      ax = subplot(1,2,1);
      superbar(ax, catFrEarly{channel_i}(2, categorySlimInds),'E',catFrErrEarly{channel_i}(2, categorySlimInds),'BarFaceColor',colorsCell);
      set(gca,'XTickLabel',categoryList(categorySlimInds),'XTickLabelRotation',45,'XTick',1:length(categorySlimInds),'TickDir','out');
      ylabel('Firing rate, Hz');
      ax = subplot(1,2,2);
      superbar(ax, catFrEarly{channel_i}(2, [faceCatNum objectCatNum bodyCatNum]),'E',catFrErrEarly{channel_i}(2, [faceCatNum objectCatNum bodyCatNum]),'BarFaceColor',{'b'; 'g'; 'r'}); %'P', [0 1 1; 1 0 0; 1 0 0]
      set(gca,'XTickLabel',{'face','object', 'body'},'XTickLabelRotation',45,'XTick',1:3,'TickDir','out');
      ylabel('Firing rate, Hz');
      saveFigure(outDir, sprintf('catPrefBarsEarly_%s_MUA_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    else
      for unit_i = 1:length(channelUnitNames{channel_i})
        ax = subplot(1,length(channelUnitNames{channel_i})+1,unit_i);
        superbar(ax,catFrEarly{channel_i}(unit_i, categorySlimInds),'E',catFrErrEarly{channel_i}(2, categorySlimInds),'BarFaceColor',colorsCell);
        set(gca,'XTickLabel',categoryList(categorySlimInds),'XTickLabelRotation',45,'XTick',1:length(categorySlimInds),'TickDir','out');
        ylabel('Firing rate, Hz');
        title(channelUnitNames{channel_i}{unit_i});
      end
      ax = subplot(1,length(channelUnitNames{channel_i})+1,length(channelUnitNames{channel_i})+1);
      superbar(ax, catFrEarly{channel_i}(end, [faceCatNum objectCatNum bodyCatNum]),'E',catFrErrEarly{channel_i}(end, [faceCatNum objectCatNum bodyCatNum]),'BarFaceColor',{'b'; 'g'; 'r'}); %'P', [0 1 1; 1 0 0; 1 0 0]
      set(gca,'XTickLabel',{'face','object','body'}, 'XTickLabelRotation',45,'XTick',1:3,'TickDir','out');
      ylabel('Firing rate, Hz');
      title('MUA');
    end
    suptitle(sprintf('%s spikes category tuning, %dms-%dms post-onset',channelNames{channel_i},frCalcOnEarly,frCalcOffEarly));
    saveFigure(outDir, sprintf('catPrefBarsEarly_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  end

  %late
  for channel_i = 1:length(spikeChannels)
    figure();
    if length(channelUnitNames{channel_i}) == 2
      ax = subplot(1,2,1);
      superbar(ax, catFrLate{channel_i}(2, categorySlimInds),'E',catFrErrLate{channel_i}(2, categorySlimInds),'BarFaceColor',colorsCell);
      set(gca,'XTickLabel',categoryList(categorySlimInds),'XTickLabelRotation',45,'XTick',1:length(categorySlimInds),'TickDir','out');
      ylabel('Firing rate, Hz');
      ax = subplot(1,2,2);
      superbar(ax, catFrLate{channel_i}(2, [faceCatNum objectCatNum bodyCatNum]),'E',catFrErrLate{channel_i}(2, [faceCatNum objectCatNum bodyCatNum]),'BarFaceColor',{'b'; 'g'; 'r'}); %'P', [0 1 1; 1 0 0; 1 0 0]
      set(gca,'XTickLabel',{'face','object', 'body'},'XTickLabelRotation',45,'XTick',1:3,'TickDir','out');
      ylabel('Firing rate, Hz');
      saveFigure(outDir, sprintf('catPrefBarsLate_%s_MUA_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    else
      for unit_i = 1:length(channelUnitNames{channel_i})
        ax = subplot(1,length(channelUnitNames{channel_i})+1,unit_i);
        superbar(ax,catFrLate{channel_i}(unit_i, categorySlimInds),'E',catFrErrLate{channel_i}(2, categorySlimInds),'BarFaceColor',colorsCell);
        set(gca,'XTickLabel',categoryList(categorySlimInds),'XTickLabelRotation',45,'XTick',1:length(categorySlimInds),'TickDir','out');
        ylabel('Firing rate, Hz');
        title(channelUnitNames{channel_i}{unit_i});
      end
      ax = subplot(1,length(channelUnitNames{channel_i})+1,length(channelUnitNames{channel_i})+1);
      superbar(ax, catFrLate{channel_i}(end, [faceCatNum objectCatNum bodyCatNum]),'E',catFrErrLate{channel_i}(end, [faceCatNum objectCatNum bodyCatNum]),'BarFaceColor',{'b'; 'g'; 'r'}); %'P', [0 1 1; 1 0 0; 1 0 0]
      set(gca,'XTickLabel',{'face','object','body'}, 'XTickLabelRotation',45,'XTick',1:3,'TickDir','out');
      ylabel('Firing rate, Hz');
      title('MUA');
    end
    suptitle(sprintf('%s spikes category tuning, %dms-%dms post-onset',channelNames{channel_i},frCalcOnLate,frCalcOffLate));
    saveFigure(outDir, sprintf('catPrefBarsLate_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  end

  % tuning curves
  %%% variables to work with are:
  % tuningCurveParams = {'humanHeadView','monkeyHeadView'};
  % tuningCurveItems = {{'humanFaceL90','humanFaceL45','humanFaceFront','humanFaceR45','humanFaceR90'},...
  %   {'monkeyFaceL90','monkeyFaceL45','monkeyFaceFront','monkeyFaceR45','monkeyFaceR90'}}; %can be images or categories
  % tuningCurveItemType = {'category','category'}; % category or image
  % tuningCurveParamValues = {[-90 -45 0 45 90], [-90 -45 0 45 90]};

  %full
  for param_i = 1:length(tuningCurveParamLabels)
    assert(strcmp(tuningCurveItemType{param_i},'category') || strcmp(tuningCurveItemType{param_i},'image'),'Invalid tuning curve type: must be category or image');
    muaTcFig = figure();
    for channel_i = 1:length(channelNames)
      channelTcFig = figure();
      for unit_i = 1:length(channelUnitNames{channel_i})
        tcFrs = zeros(length(tuningCurveItems{param_i}),1);
        tcFrErrs = zeros(length(tuningCurveItems{param_i}),1);
        if strcmp(tuningCurveItemType{param_i},'category')
          for item_i = 1:length(tuningCurveItems{param_i})
            tcFrs(item_i) = catFr{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
            tcFrErrs(item_i) = catFrErr{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
          end
        else
          for item_i = 1:length(tuningCurveItems{param_i})
            tcFrs(item_i) = imFr{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
            tcFrErrs(item_i) = imFrErr{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
          end
        end
        if length(channelUnitNames{channel_i}) == 2  % don't make the single-channel plot if no unit defined
          close(channelTcFig);
          break
        end
        subplot(1,length(channelUnitNames{channel_i}),unit_i);
        errorbar(tuningCurveParamValues{param_i},tcFrs,tcFrErrs,'linestyle','-','linewidth',4);
        xlabel(tuningCurveParamLabels{param_i}); 
        ylabel('firing rate (Hz)');
        title(channelUnitNames{channel_i}{unit_i});
      end
      suptitle(sprintf('%s, %s, %dms - %d ms post-onset',tuningCurveTitles{param_i}, channelNames{channel_i}, frCalcOn, frCalcOff));
      saveFigure(outDir, sprintf('%s_%s_Run%s',regexprep(tuningCurveTitles{param_i},' ',''),channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );

      figure(muaTcFig);
      subplot(1,length(channelNames),channel_i);
      errorbar(tuningCurveParamValues{param_i},tcFrs,tcFrErrs,'linestyle','-','linewidth',4);
      xlabel(tuningCurveParamLabels{param_i}); 
      ylabel('firing rate (Hz)');
      title(sprintf('%s MUA',channelNames{channel_i}));
    end
    suptitle(sprintf('%s, %dms - %d ms post-onset',tuningCurveTitles{param_i}, frCalcOn, frCalcOff));
    saveFigure(outDir, sprintf('%s_Run%s',regexprep(tuningCurveTitles{param_i},' ',''),runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  end

  % early
  for param_i = 1:length(tuningCurveParamLabels)
    assert(strcmp(tuningCurveItemType{param_i},'category') || strcmp(tuningCurveItemType{param_i},'image'),'Invalid tuning curve type: must be category or image');
    muaTcFig = figure();
    for channel_i = 1:length(channelNames)
      channelTcFig = figure();
      for unit_i = 1:length(channelUnitNames{channel_i})
        tcFrs = zeros(length(tuningCurveItems{param_i}),1);
        tcFrErrs = zeros(length(tuningCurveItems{param_i}),1);
        if strcmp(tuningCurveItemType{param_i},'category')
          for item_i = 1:length(tuningCurveItems{param_i})
            tcFrs(item_i) = catFrEarly{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
            tcFrErrs(item_i) = catFrErrEarly{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
          end
        else
          for item_i = 1:length(tuningCurveItems{param_i})
            tcFrs(item_i) = imFrEarly{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
            tcFrErrs(item_i) = imFrErrEarly{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
          end
        end
        if length(channelUnitNames{channel_i}) == 2  % don't make the single-channel plot if no unit defined
          close(channelTcFig);
          break
        end
        subplot(1,length(channelUnitNames{channel_i}),unit_i);
        errorbar(tuningCurveParamValues{param_i},tcFrs,tcFrErrs,'linestyle','-','linewidth',4);
        xlabel(tuningCurveParamLabels{param_i}); 
        ylabel('firing rate (Hz)');
        title(channelUnitNames{channel_i}{unit_i});
      end
      suptitle(sprintf('%s, %s, %dms - %d ms post-onset',tuningCurveTitles{param_i}, channelNames{channel_i}, frCalcOnEarly, frCalcOffEarly));
      saveFigure(outDir, sprintf('%s_%s_Early_Run%s',regexprep(tuningCurveTitles{param_i},' ',''),channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      
      figure(muaTcFig);
      subplot(1,length(channelNames),channel_i);
      errorbar(tuningCurveParamValues{param_i},tcFrs,tcFrErrs,'linestyle','-','linewidth',4);
      xlabel(tuningCurveParamLabels{param_i}); 
      ylabel('firing rate (Hz)');
      title(sprintf('%s MUA',channelNames{channel_i}));
    end
    suptitle(sprintf('%s, %dms - %d ms post-onset',tuningCurveTitles{param_i}, frCalcOnEarly, frCalcOffEarly));
    saveFigure(outDir, sprintf('%s_Early_Run%s',regexprep(tuningCurveTitles{param_i},' ',''),runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  end

  % late
  for param_i = 1:length(tuningCurveParamLabels)
    assert(strcmp(tuningCurveItemType{param_i},'category') || strcmp(tuningCurveItemType{param_i},'image'),'Invalid tuning curve type: must be category or image');
    muaTcFig = figure();
    for channel_i = 1:length(channelNames)
      channelTcFig = figure();
      for unit_i = 1:length(channelUnitNames{channel_i})
        tcFrs = zeros(length(tuningCurveItems{param_i}),1);
        tcFrErrs = zeros(length(tuningCurveItems{param_i}),1);
        if strcmp(tuningCurveItemType{param_i},'category')
          for item_i = 1:length(tuningCurveItems{param_i})
            tcFrs(item_i) = catFrLate{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
            tcFrErrs(item_i) = catFrErrLate{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
          end
        else
          for item_i = 1:length(tuningCurveItems{param_i})
            tcFrs(item_i) = imFrLate{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
            tcFrErrs(item_i) = imFrErrLate{channel_i}(unit_i,strcmp(categoryList,tuningCurveItems{param_i}{item_i}));
          end
        end
        if length(channelUnitNames{channel_i}) == 2  % don't make the single-channel plot if no unit defined
          close(channelTcFig);
          break
        end
        subplot(1,length(channelUnitNames{channel_i}),unit_i);
        errorbar(tuningCurveParamValues{param_i},tcFrs,tcFrErrs,'linestyle','-','linewidth',4);
        xlabel(tuningCurveParamLabels{param_i}); 
        ylabel('firing rate (Hz)');
        title(channelUnitNames{channel_i}{unit_i});
      end
      suptitle(sprintf('%s, %s, %dms - %d ms post-onset',tuningCurveTitles{param_i}, channelNames{channel_i}, frCalcOnEarly, frCalcOffEarly));
      saveFigure(outDir, sprintf('%s_%s_Late_Run%s',regexprep(tuningCurveTitles{param_i},' ',''),channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      
      figure(muaTcFig);
      subplot(1,length(channelNames),channel_i);
      errorbar(tuningCurveParamValues{param_i},tcFrs,tcFrErrs,'linestyle','-','linewidth',4);
      xlabel(tuningCurveParamLabels{param_i}); 
      ylabel('firing rate (Hz)');
      title(sprintf('%s MUA',channelNames{channel_i}));
    end
    suptitle(sprintf('%s, %dms - %d ms post-onset',tuningCurveTitles{param_i}, frCalcOnLate, frCalcOffLate));
    saveFigure(outDir, sprintf('%s_Late_Run%s',regexprep(tuningCurveTitles{param_i},' ',''),runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  end
end


% firing rate RF color plot
rfGrid  = unique(taskData.pictureJumps,'rows');
gridX = unique(rfGrid(:,1));
gridsize = 2*mean(gridX(2:end,1)-gridX(1:end-1,1));
% these are the x and y values at which to interpolate RF values
xi = linspace(min(rfGrid(:,1)),max(rfGrid(:,1)),200);
yi = linspace(min(rfGrid(:,2)),max(rfGrid(:,2)),200);
if taskData.RFmap
  for channel_i = 1:length(spikesByImage{1})
    for unit_i = 1:length(spikesByImage{1}{channel_i}) %TODO: BUG! if no spikes on first stimulus, can skip unit entirely (12/12/16: still true?)
      meanRF = zeros(length(rfGrid),1);
      for image_i = 1:length(pictureLabels)
        imageRF = zeros(length(rfGrid),1);
        spikeLatencyRF = zeros(length(rfGrid),1);
        %also calc evoked peak time? power-weighted latency? peak of mean evoked correlogram?
        evokedPowerRF = zeros(length(rfGrid),1);
        for grid_i = 1:length(rfGrid)
          gridPointTrials = ismember(jumpsByImage{image_i},rfGrid(grid_i,:),'rows');
          if sum(gridPointTrials) == 0 
            imageRF(grid_i) = 0;
            disp('no data for this image and location');
            continue;
          end
          trialSpikes = spikesByImage{image_i}{channel_i}{unit_i}(gridPointTrials);
          totalSpikes = 0;
          spikeLatency = 0;
          evokedPotential = zeros(size(lfpByImage{image_i}(1,channel_i,1,samPerMS*(frCalcOn+psthPre):samPerMS*(frCalcOff+psthPre)))); 
          gridPointTrialInds = 1:length(gridPointTrials);
          gridPointTrialInds = gridPointTrialInds(gridPointTrials);
          for trial_i = 1:length(trialSpikes)
            totalSpikes = totalSpikes + sum(trialSpikes(trial_i).times > frCalcOn & trialSpikes(trial_i).times < frCalcOff);
            % bad latency measure; try weighting spike times by 1/ISI
            spikeLatency = spikeLatency + mean(trialSpikes(trial_i).times(trialSpikes(trial_i).times > frCalcOn & trialSpikes(trial_i).times < frCalcOff));
            evokedPotential = evokedPotential + lfpByImage{image_i}(1,channel_i,gridPointTrialInds(trial_i),samPerMS*(frCalcOn+psthPre):samPerMS*(frCalcOff+psthPre));
          end
          imageRF(grid_i) = (1000/(frCalcOff-frCalcOn))*totalSpikes/length(trialSpikes);
          spikeLatencyRF(grid_i) = 1/sum(gridPointTrials)*spikeLatency;
          evokedPotential = 1/sum(gridPointTrials)*(evokedPotential - mean(evokedPotential));
          evokedPowerRF(grid_i) = sum(evokedPotential.^2); %todo: move out of by-unit loop
        end
        meanRF = meanRF + imageRF;
        display_map(rfGrid(:,1),rfGrid(:,2),imageRF,xi,yi,2.2857*gridsize,0,saveFig,sprintf('%s %s, %s RF',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i}),...
          [outDir sprintf('RF_%s_%s_%s_Run%s.png',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i},runNum)]);
        if calcLatencyRF
          display_map(rfGrid(:,1),rfGrid(:,2),spikeLatencyRF,xi,yi,2.2857*gridsize,0,saveFig,sprintf('%s %s, %s Latency RF',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i}),...
            [outDir sprintf('LatencyRF_%s_%s_%s_Run%s.png',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i},runNum)]);
        end
        if calcEvokedPowerRF
          display_map(rfGrid(:,1),rfGrid(:,2),evokedPowerRF,xi,yi,2.2857*gridsize,0,saveFig,sprintf('%s %s, %s Evoked Power RF',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i}),...
            [outDir sprintf('EvokedPowerRF_%s_%s_%s_Run%s.png',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i},runNum)]);
        end
        % todo: add background subtracted version
        % todo: add subplots version
        % todo: coherency RFs (cc, cpt, ptpt)
        % todo: granger RF
        % todo: bandpassed power RFs
        % todo: stimulus decodability RF
      end
    end
  end
  display_map(rfGrid(:,1),rfGrid(:,2),meanRF,xi,yi,2.2857*gridsize,0,saveFig,sprintf('Channel %d, Unit %d, Mean RF',channel_i,unit_i),...
    [outDir sprintf('MeanRF_%s_Unit%d_Run%s.png',channelNames{channel_i},unit_i,runNum)]);
  return
end  


%%%% make evoked potential plots by category

for channel_i = 1:length(lfpChannels)
  % create evoked-psth lineplot, face vs. non, one pane
  if channel_i == 1
    psthEvokedFig = figure();
    subplot(2,1,1);
    title('Face');
    xlabel('time (ms)');
    yyaxis right
    ylabel('lfp (normalized)'); % todo: units
    yyaxis left
    ylabel('firing rate (Hz)');
    hold on
    subplot(2,1,2);
    title('Non-face');
    xlabel('time (ms)');
    yyaxis right
    ylabel('lfp (normalized)'); % todo: units
    yyaxis left
    ylabel('firing rate (Hz)');
    hold on
    handlesForLegend1 = [];
    handlesForLegend2 = [];
    forLegend = {};
  end
  times = -psthPre:psthImDur+psthPost;
  faceV = squeeze(mean(lfpByCategory{faceCatNum}(1,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
  nfaceV = squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
  
  % face vs. nonface evoked potential plot, one channel
  figure();
  plot(times,faceV,times, nfaceV, 'linewidth',3);
  hold on
  plot([0, psthImDur],[1.05*min(min(faceV),min(nfaceV)), 1.05*min(min(faceV),min(nfaceV))],'color','k','linewidth',3);
  legend('face','nonface');
  title(sprintf('%s evoked potential',channelNames{channel_i}));
  xlabel('time (ms)', 'FontSize',18);
  ylabel('lfp (uV)', 'FontSize',18);
  set(gca,'fontsize',18);
  clear figData
  figData.y = vertcat(faceV,nfaceV);
  figData.x = times;
  saveFigure( outDir,sprintf('Evoked_faceVnon_%s_Run%s',channelNames{channel_i},runNum), vertcat(faceV,nfaceV), saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) ); 
  
  % contribute to evoked-psth lineplot, face vs. non, one pane
  figure(psthEvokedFig);
  subplot(2,1,1);
  yyaxis right
  lhLFP = plot(times,faceV/max(faceV),'color',chColors(channel_i), 'linestyle','-','linewidth',2);  %todo: add errorbars?
  yyaxis left
  facePSTH = squeeze(catMuaPsthByChannel(channel_i,faceCatNum,:));
  lhMUA = plot(times,facePSTH,'color',chColors(channel_i),'linestyle','--','linewidth',2); 
  handlesForLegend1 = [handlesForLegend1,lhLFP,lhMUA];
  forLegend = [forLegend,strcat(channelNames{channel_i},' LFP'),strcat(channelNames{channel_i},' MUA')];
  if channel_i == length(lfpChannels)
    legend(handlesForLegend1,forLegend);
    h = get(gca,'ylim');
    plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3,'linestyle','-');
  end
  subplot(2,1,2);
  yyaxis right
  lhLFP = plot(times,nfaceV/max(nfaceV),'color',chColors(channel_i), 'linestyle', '-','linewidth',2);
  yyaxis left
  nfacePSTH = squeeze(catMuaPsthByChannel(channel_i,nonfaceCatNum,:));
  lhMUA = plot(times,nfacePSTH,'color',chColors(channel_i),'linestyle','--','linewidth',2);
  handlesForLegend2 = [handlesForLegend2, lhLFP, lhMUA];
  if channel_i == length(lfpChannels)
    legend(handlesForLegend2,forLegend);
    h = get(gca,'ylim');
    plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3,'linestyle','-');
    saveFigure(outDir, sprintf('faceVnonPsthEvokedOverlay_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  end
  
  
  %lfp category psth
  figure();
  channelCatEvoked = zeros(length(categoryListSlim),length(times));
  for i = 1:length(categoryListSlim)
    channelCatEvoked(i,:) = squeeze(mean(lfpByCategory{categorySlimInds(i)}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
    plot(times, channelCatEvoked(i,:), 'color',colors(i), 'linewidth',3);
    hold on
  end
  legend(categoryListSlim);
  h = get(gca,'ylim');
  plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
  hold off
  title(sprintf('%s evoked potentials',channelNames{channel_i}), 'FontSize',18);
  xlabel('time after stimulus (ms)', 'FontSize',18);
  ylabel('lfp (uV)', 'FontSize',18);
  set(gca,'fontsize',18);
  clear figData
  figData.y = squeeze(mean(lfpByCategory{categorySlimInds(i)}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
  figData.x = times;
  saveFigure(outDir,sprintf('Evoked_byCat_%s_Run%s',channelNames{channel_i},runNum), channelCatEvoked, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  
  % make lfp - psth subplot
  categoryPSTH = squeeze(catMuaPsthByChannel(channel_i,categorySlimInds,:));
  figure();
  ah1 = subplot(2,1,1);
  plotPSTH(categoryPSTH, [], psthPre, psthPost, psthImDur, 'color', psthTitle, categoryListSlim, psthColormap );
  yyaxis right
  plot(times,mean(categoryPSTH,1),'Color',[0.8,0.8,0.9],'LineWidth',4);
  ylabel('Mean PSTH (Hz)');
  hold off
  ah2 = subplot(2,1,2);
  for i = 1:length(categoryListSlim)
    plot(times,channelCatEvoked(i,:), 'color',colors(i), 'linewidth',3);
    hold on
  end
  h = get(gca,'ylim');
  plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
  hold off
  legend(categoryListSlim,'Location','northeastoutside','FontSize',10);
  xlabel('time after stimulus (ms)', 'FontSize',14);
  ylabel('lfp (uV)', 'FontSize',14);
  set(gca,'fontsize',14);
  % align x axes with colorbar; see http://stackoverflow.com/questions/5259853/matlab-how-to-align-the-axes-of-subplots-when-one-of-them-contains-a-colorbar
  % note: drawnow line synchronizes the rendering thread with the main
  drawnow;
  pos1 = get(ah1,'Position');
  pos2 = get(ah2,'Position');
  pos2(1) = max(pos1(1),pos2(1)); % right limit
  pos1(1) = max(pos1(1),pos2(1));
  pos2(3) = min(pos1(3),pos2(3)); % left limit
  pos1(3) = min(pos1(3),pos2(3));
  set(ah2,'Position',pos2);
  set(ah1,'Position',pos1);
  title(ah1,sprintf('%s PSTH and Evoked Potentials',channelNames{channel_i}));
  % put the data in a struct for saving
  clear figData
  figData.ax11.z = categoryPSTH;
  figData.ax11.x = times;
  figData.ax21.y = channelCatEvoked;
  figData.ax21.x = times;
  saveFigure(outDir,sprintf('PSTH_Evoked_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  % single trial evoked lfps; probably better as subplot
  figure();
  subplot(2,1,1)
  hold on
  ydata = zeros(50,length(times));
  for i = 1:50
    plot(times, squeeze(lfpByCategory{faceCatNum}(1,channel_i,i,lfpPaddedBy+1:end-lfpPaddedBy)));
    ydata(i,:) = squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,i,lfpPaddedBy+1:end-lfpPaddedBy));
  end;
  h = get(gca,'ylim');
  plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
  hold off
  title(sprintf('single trial face LFPs, %s', channelNames{channel_i}), 'FontSize',18);
  xlabel('time after stimulus (ms)', 'FontSize',18);
  ylabel('voltage (uV)', 'FontSize',18);
  set(gca,'fontsize',18);
  xlim([min(times) max(times)]);
  clear figData
  figData.ax11.y = ydata;
  figData.ax11.x = times;
  
  subplot(2,1,2);
  hold on
  ydata = zeros(50,length(times));
  for i = 1:50
    plot(times, squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,i,lfpPaddedBy+1:end-lfpPaddedBy)));
    ydata(i,:) = squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,i,lfpPaddedBy+1:end-lfpPaddedBy));
  end;
  h = get(gca,'ylim');
  plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
  hold off
  title(sprintf('single trial nonface LFPs, %s', channelNames{channel_i}), 'FontSize',18);
  xlabel('time after stimulus (ms)', 'FontSize',18);
  ylabel('voltage (uV)', 'FontSize',18);
  xlim([min(times) max(times)]);
  figData.ax21.y = ydata;
  figData.ax21.x = times;
  saveFigure(outDir,sprintf('Evoked_singleTrials_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  
  % spectra
  [faceS,faceF, faceE] = mtspectrumc(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))', chr_params);
  [nfaceS,nfaceF, nfaceE] = mtspectrumc(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))', chr_params);
  
  figure();
  plot(1000*faceF,log(faceS),'linewidth',3,'marker','o');
  hold on
  plot(1000*nfaceF,log(nfaceS),'linewidth',3,'color','r','marker','o');
  hold off
  title(sprintf('%s evoked power spectrum',channelNames{channel_i}), 'FontSize',18);
  xlabel('frequency (Hz)', 'FontSize',18);
  ylabel('voltage, log(uV)', 'FontSize',18);
  legend('face','non');
  set(gca,'fontsize',18);
  clear figData
  figData.y = vertcat(faceS,nfaceS);
  figData.x = vertcat(faceF,nfaceF);
  saveFigure(outDir,sprintf('spectrum_faceVnon_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  
  figure();
  specModel = fitlm(log(1000*faceF(6:end)), log(.5*(faceS(6:end)+nfaceS(6:end))));
  m = specModel.Coefficients.Estimate(2);
  y0 = specModel.Coefficients.Estimate(1);
  h1 = loglog(1000*faceF,exp(m*log(1000*faceF) + y0), 'k--');
  hold on
  h2 = loglog(1000*faceF,faceS,'linewidth',3,'marker','o');
  h3 = loglog(1000*nfaceF,nfaceS,'linewidth',3,'color','r','marker','o');
  hold off
  title(sprintf('%s evoked power spectrum, loglog',channelNames{channel_i}), 'FontSize',18);
  xlim([1000*min(faceF) 1000*max(faceF)]); %todo: fix error
  xlabel('frequency (Hz)', 'FontSize',18);
  ylabel('voltage (uV)', 'FontSize',18);
  legend([h2;h3;h1],{'face','nonface',sprintf('fit,m = %.3f',m)});
  set(gca,'fontsize',18);
  clear figData
  figData.y = vertcat(faceS,nfaceS);
  figData.x = vertcat(faceF,nfaceF);
  saveFigure(outDir,sprintf('spectrum_log_faceVnon_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
  
  
  %%%%%% TODO %%%%%%
  % todo: multichannel evoked lineplot, categories, subplots
  % todo: find local maxes and mins in evoked and psth
  % psth category lineplot (probably in subplot with evoked by category)

  % image-wise time-frequency plots, spikes and lfp
  if imageTF
    for image_i = 1:length(pictureLabels)
      % todo: put in dB conversion and specgramrowave option
      [S,t,f,R]=mtspecgrampt(spikesByImage{image_i}{channel_i}{end},movingWin,chr_params); %last optional param not used: fscorr
      t = t - lfpAlignParams.msPreAlign;
      f = 1000*f;
      figure();
      imagesc(t,f,S'); %TODO: fix to match cat version, which is correct!!
      set(gca,'Ydir','normal'); % also need 'Yscale,'log'??
      xlabel('Time'); % todo: units
      ylabel('Frequency (Hz)');
      c = colorbar();
      ylabel(c,'Power'); % todo: check units
      hold on
      draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
      draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
      title(sprintf('%s MUA Time-Frequency, %s',channelNames{channel_i},pictureLabels{image_i}));
      clear FigData
      figData.x = t;
      figData.y = f;
      figData.z = S';
      saveFigure(outDir,sprintf('TF_MUA_%s_%s.fig',channelNames{channel_i},pictureLabels{image_i}), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      
      [S,t,f]=mtspecgramc(squeeze(lfpByImage{image_i}(1,channel_i,:,:))',movingWin,chr_params);
      t = t - lfpAlignParams.msPreAlign;
      f = 1000*f;
      figure();
      imagesc(t,f,S');
      set(gca,'Ydir','normal'); % also need 'Yscale,'log'??
      xlabel('Time'); % todo: units
      ylabel('Frequency (Hz)');
      c = colorbar();
      ylabel(c,'Power'); % todo: check units
      hold on
      draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
      draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
      title(sprintf('%s LFP Time-Frequency, %s',channelNames{channel_i},pictureLabels{image_i}));
      clear figData
      figData.x = t;
      figData.y = f;
      figData.z = S;
      saveFigure(outDir,sprintf('TF_LFP_%s_%s',channelNames{channel_i},pictureLabels{image_i}), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
  % category-wise time-frequency plots, spikes and lfp
  
  if catTF
    if channel_i == 1
      lfpTfFig = figure();
      muaTfFig = figure();
    end
    for cat_i = 1:length(categoryList)
      if cat_i ~= faceCatNum && cat_i ~= nonfaceCatNum
        continue
      end
      [S,t,f,R]=mtspecgrampt(spikesByCategoryForTF{cat_i}{channel_i}{end},movingWin,chr_params); %last optional param not used: fscorr
      t = t - lfpAlignParams.msPreAlign;
      f = 1000*f;
      totalPower = sum(S,2)';
      figure();
      if specgramRowAve
        for i = 1:size(S,2)
          S(:,i) = S(:,i)/mean(S(:,i)); 
        end
        imagesc(t,f,S'); axis xy; c = colorbar();
        ylabel(c,'Row-Normalized Power');
      else
        S = 10*log10(S);
        imagesc(t,f,S'); axis xy; c = colorbar();
        ylabel(c,'Power (dB)');
      end
      xlabel('Time (ms)'); 
      ylabel('Frequency (Hz)');
      hold on
      draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
      draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
      title(sprintf('%s LFP Time-Frequency, %s',channelNames{channel_i},pictureLabels{image_i}));
      yyaxis right
      plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
      ylabel('Integrated Power');
      hold off
      title(sprintf('%s MUA Time-Frequency, %s',channelNames{channel_i},categoryList{cat_i}));
      clear figData
      figData.x = t;
      figData.y = f;
      figData.z = S;
      saveFigure(outDir,sprintf('TF_MUA_%s_%s',channelNames{channel_i},categoryList{cat_i}), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      % contribute to shared figure
      figure(muaTfFig);
      if cat_i == faceCatNum
        subplot(3,2,2*channel_i-1);
        forTitle = sprintf('Face MUA TF, %s',channelNames{channel_i});
      else
        subplot(3,2,2*channel_i);
        forTitle = sprintf('Nonface MUA TF, %s',channelNames{channel_i});
      end
      imagesc(t,f,S'); axis xy; c = colorbar(); ylabel(c,'Power (dB)'); %todo: fix unit when row-normalizing
      hold on
      draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
      draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
      title(sprintf('%s LFP Time-Frequency, %s',channelNames{channel_i},pictureLabels{image_i}));
      xlabel('Time (ms)'); 
      ylabel('Frequency (Hz)');
      title(forTitle);
      
      [S,t,f]=mtspecgramc(squeeze(lfpByCategory{cat_i}(1,channel_i,:,:))',movingWin,chr_params);
      t = t - lfpAlignParams.msPreAlign;
      f = 1000*f;
      totalPower = sum(S,2)';
      
      figure();
      if specgramRowAve
        for i = 1:size(S,2)
          S(:,i) = S(:,i)/mean(S(:,i));
        end
        imagesc(t,f,S'); axis xy; c = colorbar();
        ylabel(c,'Row-Normalized Power');
      else
        S = 10*log10(S);
        imagesc(t,f,S'); axis xy; c = colorbar();
        ylabel(c,'Power (dB)');
      end
      xlabel('Time (ms)'); 
      ylabel('Frequency (Hz)');
      hold on
      draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
      draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
      yyaxis right
      plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
      ylabel('Integrated Power');
      hold off
      title(sprintf('%s LFP Time-Frequency, %s',channelNames{channel_i},categoryList{cat_i}));
      clear figData
      figData.x = t;
      figData.y = f;
      figData.z = S';
      saveFigure(outDir,sprintf('TF_LFP_%s_%s.fig',channelNames{channel_i},categoryList{cat_i}), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      
      % contribute to shared figure
      figure(lfpTfFig);
      if cat_i == faceCatNum
        subplot(3,2,2*channel_i-1);
        forTitle = sprintf('Face LFP TF, %s',channelNames{channel_i});
      else
        subplot(3,2,2*channel_i);
        forTitle = sprintf('Nonface LFP TF, %s',channelNames{channel_i});
      end
      imagesc(t,f,S'); axis xy; c = colorbar(); ylabel(c,'Power (dB)'); %todo: fix unit when row-normalizing
      hold on
      draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
      draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
      xlabel('Time (ms)'); 
      ylabel('Frequency (Hz)');
      title(forTitle);
    end
    
    % todo: convert to saveFigure and build its fig data as we go
    if saveFig && channel_i == length(lfpChannels) && cat_i == max(faceCatNum, nonfaceCatNum)
      figure(muaTfFig);
      savefig(strcat(outDir,sprintf('TF_MUA_%s_%s.fig',channelNames{channel_i},categoryList{cat_i})));
      export_fig([outDir sprintf('TF_MUA_%s_%s.png',channelNames{channel_i},categoryList{cat_i})],'-m1.2','-transparent','-opengl');
      figure(lfpTfFig);
      savefig(lfpTfFig,strcat(outDir,sprintf('TF_LFP_%s_%s.fig',channelNames{channel_i},categoryList{cat_i})));
      export_fig([outDir sprintf('TF_LFP_%s_%s.png',channelNames{channel_i},categoryList{cat_i})],'-m1.2','-transparent','-opengl');
    end
  end
  
  drawnow;
  
  %%% across channels
  if crossTF
    for channel2_i = channel_i:length(lfpChannels)
      % face vs. nonface spike field coherence, within and across channels
      % channel_i spike -- channel_i field
      useJacknife = 0;
      if useJacknife
        chr_params.err = [2 .05];
        % todo: eliminate 'faceMuaML' etc.
        [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
          spikesByCategoryForTF{faceCatNum}{channel_i}{end},chr_params);
        [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
          spikesByCategory{nonfaceCatNum}{channel_i}{end},chr_params);
        errs = vertcat(faceErrs,nfaceErrs);
      else
        [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
          spikesByCategoryForTF{faceCatNum}{channel_i}{end},chr_params);
        [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
          spikesByCategoryForTF{nonfaceCatNum}{channel_i}{end},chr_params);
        errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
      end
      figure();
      mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
      legend('face', 'nonface');
      xlabel('frequency (Hz)');
      ylabel('coherency');
      title(sprintf('%s spike - %s field coherence',channelNames{channel_i},channelNames{channel_i}),'FontSize',18);
      clear figData
      figData.x = vertcat(fface, fnface);
      figData.y = [Cface, Cnface];
      saveFigure(outDir,sprintf('coh_MUA-LFP_%s_%s_FACEvsNON_Run%s',channelNames{channel_i},channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      if channel2_i > channel_i 
        % channel2_i spike -- channel_i field
        useJacknife = 0;
        if useJacknife
          chr_params.err = [2 .05];
          [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
            spikesByCategoryForTF{faceCatNum}{channel2_i}{end},chr_params);
          [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
            spikesByCategoryForTF{nonfaceCatNum}{channel2_i}{end},chr_params);
          errs = vertcat(faceErrs,nfaceErrs);
        else
          [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
            spikesByCategoryForTF{faceCatNum}{channel2_i}{end},chr_params);
          [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
            spikesByCategoryForTF{nonfaceCatNum}{channel2_i}{end},chr_params);
          errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
        end
        figure();
        mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
        legend('face', 'nonface');
        xlabel('frequency (Hz)');
        ylabel('coherency');
        title(sprintf('%s spike - %s field coherence',channelNames{channel2_i},channelNames{channel_i}),'FontSize',18);
        clear figData
        figData.x = vertcat(fface, fnface);
        figData.y = [Cface, Cnface];
        saveFigure(outDir,sprintf('coh_MUA-LFP_%s_%s_FACEvsNON_Run%s',channelNames{channel2_i},channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        % channel_i spike -- channel2_i field
        useJacknife = 0;
        if useJacknife
          chr_params.err = [2 .05];
          [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
            spikesByCategoryForTF{faceCatNum}{channel_i}{end},chr_params);
          [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
            spikesByCategory{nonfaceCatNum}{channel_i}{end},chr_params);
          errs = vertcat(faceErrs,nfaceErrs);
        else
          [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
            spikesByCategoryForTF{faceCatNum}{channel_i}{end},chr_params);
          [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
            spikesByCategory{nonfaceCatNum}{channel_i}{end},chr_params);
          errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
        end
        figure();
        mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
        legend('face', 'nonface');
        xlabel('frequency (Hz)');
        ylabel('coherency');
        title(sprintf('%s spike - %s field coherence',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = vertcat(fface, fnface);
        figData.y = [Cface, Cnface];
        saveFigure(outDir,sprintf('coh_MUA-LFP_%s_%s_FACEvsNON_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        

        % field-field
        [Cface,phiface,S12,S1,S2,fface,confCface,phistdface]=coherencyc(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
          squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',chr_params);
        [Cnface,phinface,S12,S1,S2,fnface,confCnface,phistdnface]=coherencyc(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
          squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',chr_params);
        errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
        figure();
        mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
        legend('face', 'nonface');
        xlabel('frequency (Hz)');
        ylabel('coherency');
        title(sprintf('%s field - %s field coherence',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = vertcat(fface, fnface);
        figData.y = [Cface, Cnface];
        saveFigure(outDir,sprintf('coh_LFP-LFP_%s_%s_FACEvsNON_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        
        
        % field-field time-frequency coherency, face
        [Cface,phiface,S12,S1,S2,tface,fface,confCface,phistdface]=cohgramc(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
          squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',movingWin,chr_params);
        t = tface - lfpAlignParams.msPreAlign;
        f = 1000*fface;
        
        figure(); 
        imagesc(t,f,Cface'); axis xy
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'Coherency');
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s field - %s field coherence, face',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = t;
        figData.y = f;
        figData.z = Cface';
        saveFigure(outDir,sprintf('coh_TF_LFP-LFP_%s_%s_FACE_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        
        figure(); 
        imagesc(t,f,phiface'); axis xy
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'phase');
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s field - %s field phase, face',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = t;
        figData.y = fface;
        figData.z = phiface';
        saveFigure(outDir,sprintf('phase_TF_LFP-LFP_%s_%s_FACE_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        
        % lfp-lfp time-frequency coherency, nonface
        [Cnface,phinface,S12,S1,S2,tnface,fnface,confCnface,phistdnface]=cohgramc(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
          squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',movingWin,chr_params);
        t = tnface - lfpAlignParams.msPreAlign;
        f = 1000*fnface;
        figure();
        imagesc(t,f,Cnface'); axis xy
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'Coherency'); 
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s field - %s field coherence, nonface',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = t;
        figData.y = f;
        figData.z = Cnface';
        saveFigure(outDir,sprintf('coh_TF_LFP-LFP_%s_%s_NONFACE_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        
        figure(); 
        imagesc(t,f,phinface'); axis xy
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'phase');
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s field - %s field phase, nonface',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = t;
        figData.y = f;
        figData.z = phinface';
        saveFigure(outDir,sprintf('phase_TF_LFP-LFP_%s_%s_NONFACE_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        
        % spike-spike coherency, face vs. non
        [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypt(spikesByCategory{faceCatNum}{channel_i}{end},...
          spikesByCategory{faceCatNum}{channel2_i}{end},chr_params,0,1);
        [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(spikesByCategory{faceCatNum}{channel_i}{end},...
          spikesByCategory{nonfaceCatNum}{channel2_i}{end},chr_params,0,1);
        errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
        figure();
        mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
        legend('face', 'nonface');
        xlabel('frequency (Hz)');
        ylabel('coherency');
        title(sprintf('%s Spike - %s Spike coherence',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = vertcat(fface, fnface); %todo: fix; add freq cutoff array slice [ applies to all coherence figures ]
        figData.y = [Cface, Cnface];
        saveFigure(outDir,sprintf('coh_MUA-MUA_%s_%s_FACEvsNON_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        
        % spike-spike time-frequency coherency, face vs. non
        
        [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgrampt(spikesByCategoryForTF{faceCatNum}{channel_i}{end},...
          spikesByCategoryForTF{faceCatNum}{channel2_i}{end}, movingWin,chr_params);
        
        t = tface - lfpAlignParams.msPreAlign;
        f = 1000*fface;
        
        figure(); 
        imagesc(t,f,Cface'); axis xy
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'Coherency');
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s spike - %s spike coherence, face',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = t;
        figData.y = f;
        figData.z = Cface';
        saveFigure(outDir,sprintf('coh_TF_MUA-MUA_%s_%s_FACE_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        
        figure(); 
        imagesc(t,f,phiface'); axis xy
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'phase');
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s spike - %s spike phase, face',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = t;
        figData.y = fface;
        figData.z = phiface';
        saveFigure(outDir,sprintf('phase_TF_MUA-MUA_%s_%s_FACE_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        
        % nonface
        [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCnface,phistdnface]=cohgrampt(spikesByCategoryForTF{nonfaceCatNum}{channel_i}{end},...
          spikesByCategoryForTF{nonfaceCatNum}{channel2_i}{end}, movingWin,chr_params);
        
        t = tnface - lfpAlignParams.msPreAlign;
        f = 1000*fnface;
        
        figure(); 
        imagesc(t,f,Cnface'); axis xy
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'Coherency');
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s spike - %s spike coherence, nonface',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = t;
        figData.y = f;
        figData.z = Cnface';
        saveFigure(outDir,sprintf('coh_TF_MUA-MUA_%s_%s_NONFACE_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        
        figure(); 
        imagesc(t,f,phinface'); axis xy
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'phase');
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s spike - %s spike phase, nonface',channelNames{channel_i},channelNames{channel2_i}),'FontSize',18);
        clear figData
        figData.x = t;
        figData.y = fNface;
        figData.z = phinface';
        saveFigure(outDir,sprintf('phase_TF_MUA-MUA_%s_%s_NONFACE_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        % todo: granger (full, passbands)
      end
    end
  end 
end
end

