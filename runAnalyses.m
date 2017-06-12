function [  ] = runAnalyses( analysisParamFilename, spikesByChannel, lfpData, analogInData, taskData, taskDataAll, psthImDur, preAlign, postAlign, ...
  categoryList, pictureLabels, jumpsByImage, spikesByImage, psthEmptyByImage, spikesByCategory, psthEmptyByCategory,...
  spikesByImageForTF, spikesByCategoryForTF, lfpByImage, lfpByCategory, channelUnitNames, stimTiming, picCategories, onsetsByImage, OnsetsByCategory)
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

if ~calcSwitch.spikeTimes %use 1 ms bins for spikes
  spikesByCategoryBinned = cell(size(spikesByCategory));
  assert((movingWin(1)/2 > 3*smoothingWidth),'Error: current implementation assumes that movingWin/2 > 3*psthSmoothingWidth. Not true here');
  for cat_i = 1:length(spikesByCategory)
    spikesByCategoryBinned{cat_i} = cell(length(channelNames));
    for channel_i = 1:length(channelNames)
      spikesByCategoryBinned{cat_i}{channel_i} = cell(length(channelUnitNames{channel_i}));
      for unit_i = 1:length(channelUnitNames{channel_i})
          spikesByCategoryBinned{cat_i}{channel_i}{unit_i} = binspikes(spikesByCategory{cat_i}{channel_i}{unit_i},1,[-1*(psthPre+movingWin(1)/2), psthImDur+psthPost+movingWin(1)/2])';
      end
    end
  end
  spikesByImageBinned = cell(size(spikesByCategory));
  for image_i = 1:length(spikesByImage)
    spikesByImageBinned{image_i} = cell(length(channelNames));
    for channel_i = 1:length(channelNames)
      spikesByImageBinned{image_i}{channel_i} = cell(length(channelUnitNames{channel_i}));
      for unit_i = 1:length(channelUnitNames{channel_i})
          spikesByImageBinned{image_i}{channel_i}{unit_i} = binspikes(spikesByImage{image_i}{channel_i}{unit_i},1,[-1*(psthPre+movingWin(1)/2), psthImDur+psthPost+movingWin(1)/2])';
      end
    end
  end
end


%spike psth color plot
if makeImPSTH
  if ~calcSwitch.spikeTimes
    filterPoints = -3*smoothingWidth:3*smoothingWidth;
    smoothingFilter = exp(-1*filterPoints.^2/(2*smoothingWidth^2));
    smoothingFilter = smoothingFilter/sum(smoothingFilter);
  end
  
  for channel_i = 1:length(spikesByImage{1})
    for unit_i = 1:length(spikesByImage{1}{channel_i}) 
      if length(spikesByImage{1}{channel_i}) == 2 && unit_i == 1 %if no isolated unit defined, plot just MUA, not also unsorted (since it's identical)
        continue;
      end
      imagePSTH = zeros(length(pictureLabels),psthPre+1+psthImDur+psthPost);
      for image_i = 1:length(pictureLabels)
        if ~psthEmptyByImage{image_i}{channel_i}{unit_i}            
          if calcSwitch.spikeTimes
            paddedPsth = 1000*psth(spikesByImage{image_i}{channel_i}{unit_i},smoothingWidth,'n',[-preAlign postAlign],0,-preAlign:postAlign);
            imagePSTH(image_i,:) = paddedPsth(3*smoothingWidth+1:end-3*smoothingWidth);
          else %use spike bins
            paddedPsth = 1000*conv(mean(spikesByImageBinned{image_i}{channel_i}{unit_i},1),smoothingFilter,'same');
            imagePSTH(image_i,:) = paddedPsth(movingWin(1)/2+1:end-movingWin(1)/2);
          end
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
          if calcSwitch.spikeTimes
            paddedPsth = 1000*psth(spikesByCategory{cat_i}{channel_i}{unit_i},smoothingWidth,'n',[-preAlign postAlign],0,-preAlign:postAlign);
            categoryPSTH(cat_i,:) = paddedPsth(3*smoothingWidth+1:end-3*smoothingWidth);
          else
            paddedPsth= 1000*conv(mean(spikesByCategoryBinned{cat_i}{channel_i}{unit_i},1),smoothingFilter,'same');
            categoryPSTH(cat_i,:) = paddedPsth(movingWin(1)/2+1:end-movingWin(1)/2);
          end
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

colors = ['b','c','y','g','m','r','k'];
colorsCell = {'b';'c';'y';'g';'m';'r';'k'};
chColors = ['b','g','m'];

% todo: make loop?
[imSpikeCounts, imFr, imFrErr] = spikeCounter(spikesByImage, frCalcOn, frCalcOff);
[imSpikeCountsEarly, imFrEarly, imFrErrEarly] = spikeCounter(spikesByImage, frCalcOnEarly, frCalcOffEarly);
[imSpikeCountsLate, imFrLate, imFrErrLate] = spikeCounter(spikesByImage, frCalcOnLate, frCalcOffLate);

if ~isempty(spikesByCategory)
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
  imageSlimCats = zeros(length(pictureLabels),1);
  imageSlimCatColors = {};
  %categoryListFOB = {'face','object','body'};
  %imageFatCats = zeros(length(pictureLabels),1);
  for image_i = 1:length(pictureLabels)
    for slimCat_i = 1:length(categorySlimInds)
      if any(strcmp(picCategories{image_i},categoryListSlim{slimCat_i}))
        imageSlimCats(image_i) = slimCat_i;
        imageSlimCatColors = vertcat(imageSlimCatColors,colors(slimCat_i));
        break
      end
    end
    if length(imageSlimCatColors) < image_i
      imageSlimCats(image_i) = 0;
      imageSlimCatColors = vertcat(imageSlimCatColors,'w');
      fprintf('no slim category match found for %s\n',pictureLabels{image_i});
    end
  end
end

if ~taskData.RFmap
  % preferred images
  for channel_i = 1:length(spikeChannels)
    for unit_i = 1:length(channelUnitNames{channel_i})
      [imageSortedRates, imageSortOrder] = sort(imFr{channel_i}(unit_i,:),2,'descend');
      imFrErrSorted = imFrErr{channel_i}(unit_i,imageSortOrder);
      sortedImageLabels = pictureLabels(imageSortOrder);
      sortedImageSlimCatColors = imageSlimCatColors(imageSortOrder);
      fprintf('\n\n\nPreferred Images: %s, %s\n\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i});
      for i = 1:min(10,length(pictureLabels))
        fprintf('%d) %s: %.2f +/- %.2f Hz\n',i,sortedImageLabels{i},imageSortedRates(i),imFrErrSorted(i));
      end
      fprintf('\nLeast Preferred Images: %s, %s\n\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i});
      for i = 1:min(5,length(pictureLabels))
        fprintf('%d) %s: %.2f +/- %.2f Hz\n',i,sortedImageLabels{end-i},imageSortedRates(end-i), imFrErrSorted(end-i));
      end
      % preferred images raster plot
      if plotSwitch.prefImRaster
        figure();
        raster(spikesByImage(imageSortOrder(1:10)), sortedImageLabels(1:10), psthPre, psthPost, psthImDur, stimTiming.ISI, channel_i, unit_i, colors);
        title(sprintf('Preferred Images, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
        saveFigure(outDir, sprintf('prefImRaster_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      end
      % preferred images raster-evoked overlay
      if plotSwitch.prefImRasterEvokedOverlay
        figure();
        rasterEvoked(spikesByImage(imageSortOrder(1:10)), lfpByImage(imageSortOrder(1:10)), sortedImageLabels(1:10), psthPre, psthPost, psthImDur, stimTiming.ISI, lfpPaddedBy, channel_i, colors, 1)
        title(sprintf('Preferred Images, from top, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
        saveFigure(outDir, sprintf('prefImRaster-LFP_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      end
      % preferred images raster-evoked overlay, with other channels
      if plotSwitch.prefImMultiChRasterEvokedOverlay
        figure();
        rasterEvokedMultiCh(spikesByImage(imageSortOrder(1:10)), lfpByImage(imageSortOrder(1:10)), sortedImageLabels(1:10), psthPre, psthPost, psthImDur, stimTiming.ISI, lfpPaddedBy, 1:length(lfpChannels), channelNames, colors)
        title(sprintf('Preferred Images, from top, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
        saveFigure(outDir, sprintf('prefImRaster-LFP-MultiChannel_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      end
      % image preference barplot
      if plotSwitch.imageTuningSorted
        figure();
        superbar(imageSortedRates,'E',imFrErrSorted,'BarFaceColor',sortedImageSlimCatColors);
        set(gca,'XTickLabel',sortedImageLabels,'XTickLabelRotation',45,'XTick',1:length(pictureLabels),'TickDir','out');
        ylabel('Firing rate, Hz');
        title(sprintf('Image tuning, sorted, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
        saveFigure(outDir, sprintf('imageTuningSorted_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        % todo: add legend for category-color map
      end
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
    for i = 1:min(10,length(pictureLabels))
      fprintf('%d) %s: %.2f Hz\n',i,sortedImageLabels{i},imageSortedRates(i));
    end
    fprintf('\n\nMulti-channel Least preferred Images, Maximin\n\n');
    for i = 1:min(5,length(pictureLabels))
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
    for i = 1:min(10,length(pictureLabels))
      fprintf('%d) %s: %.2f Hz\n',i,sortedImageLabels{i},imageSortedRates(i));
    end
    fprintf('\n\nMulti-channel Least preferred Images, Channel-Normalized Mean\n\n');
    for i = 1:min(5,length(pictureLabels))
      fprintf('%d) %s: %.2f Hz \n',i,sortedImageLabels{end-i},imageSortedRates(end-i));
    end
  end
  % face selectivity index
  if calcSwitch.faceSelectIndex
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
  end
  % face selectivity index, early
  if calcSwitch.faceSelectIndexEarly
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
  end
  % face selectivity index, late
  if calcSwitch.faceSelectIndexLate
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
  end

  % category preference bar plot
  % full
  if plotSwitch.categoryPrefBarPlot
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
  end

  % early
  if plotSwitch.categoryPrefBarPlotEarly
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
  end

  %late
  if plotSwitch.categoryPrefBarPlotLate
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
  end

  % tuning curves
  %%% variables to work with are:
  % tuningCurveParams = {'humanHeadView','monkeyHeadView'};
  % tuningCurveItems = {{'humanFaceL90','humanFaceL45','humanFaceFront','humanFaceR45','humanFaceR90'},...
  %   {'monkeyFaceL90','monkeyFaceL45','monkeyFaceFront','monkeyFaceR45','monkeyFaceR90'}}; %can be images or categories
  % tuningCurveItemType = {'category','category'}; % category or image
  % tuningCurveParamValues = {[-90 -45 0 45 90], [-90 -45 0 45 90]};

  %full
  if plotSwitch.tuningCurves
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
  end

  % early
  if plotSwitch.tuningCurvesEarly
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
  end

  % late
  if plotSwitch.tuningCurvesLate
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
end


% firing rate RF color plot
if taskData.RFmap
  rfGrid  = unique(taskData.pictureJumps,'rows');
  gridX = unique(rfGrid(:,1));
  gridsize = 2*mean(gridX(2:end,1)-gridX(1:end-1,1));
  % these are the x and y values at which to interpolate RF values
  xi = linspace(min(rfGrid(:,1)),max(rfGrid(:,1)),200);
  yi = linspace(min(rfGrid(:,2)),max(rfGrid(:,2)),200);
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
        if plotSwitch.LatencyRF
          display_map(rfGrid(:,1),rfGrid(:,2),spikeLatencyRF,xi,yi,2.2857*gridsize,0,saveFig,sprintf('%s %s, %s Latency RF',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i}),...
            [outDir sprintf('LatencyRF_%s_%s_%s_Run%s.png',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i},runNum)]);
        end
        if plotSwtich.calcEvokedPowerRF
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

% evokedTF says to do TF analysis on the trial average psth and evoked potentials
if calcSwitch.meanEvokedTF && calcSwitch.spikeTimes
  % note: because I'm currently working with spike times rather than
  % binned spikes, I need to do this loop to make a struct whose 'times'
  % field has all the spikes from a category
  allSpikesByCategoryForTF = cell(length(spikesByCategory),1);
  for cat_i = 1:length(spikesByCategory)
    catChannelSpikes = cell(length(channelNames),1);
    for channel_i = 1:length(channelNames)
      channelUnitSpikes = cell(length(channelUnitNames{channel_i}),1);
      for unit_i = 1:length(channelUnitNames{channel_i})
        s.times = [];
        for trial_i = 1:length(spikesByCategoryForTF{cat_i}{channel_i}{unit_i})
          s.times = vertcat(s.times,spikesByCategoryForTF{cat_i}{channel_i}{unit_i}(trial_i).times);
        end
        channelUnitSpikes{unit_i} = s;
      end
      catChannelSpikes{channel_i} = channelUnitSpikes;
    end
    allSpikesByCategoryForTF{cat_i} = catChannelSpikes;
  end
end

% if we're going to need induced psth's and lfp's, build them now
if (calcSwitch.inducedSpectra || calcSwitch.inducedCatTF) && ~calcSwitch.spikeTimes && ~calcSwitch.inducedTrialMagnitudeCorrection
  % subtract the category mean psth from every trial to obtain induced psth
  spikesByCategoryBinnedInduced = spikesByCategoryBinned;
  for cat_i = 1:length(spikesByCategoryBinned)
    for channel_i = 1:length(channelNames)
      for unit_i = 1:length(channelUnitNames{channel_i})
        spikesByCategoryBinnedInduced{cat_i}{channel_i}{unit_i} = spikesByCategoryBinned{cat_i}{channel_i}{unit_i} - mean(spikesByCategoryBinned{cat_i}{channel_i}{unit_i},1); %note: matlab automatically applies subtraction to every row
      end
    end
  end
  lfpByCategoryInduced = lfpByCategory;
  for cat_i = 1:length(lfpByCategory)
    lfpByCategoryInduced{cat_i} = lfpByCategory{cat_i} - mean(lfpByCategory{cat_i},3); %note: matlab automatically applies subtraction to every trial/channel appropriately
  end
end

% if we're going to need induced, response-magnitude-corrected psth's and lfp's, build them now
if (calcSwitch.inducedSpectra || calcSwitch.inducedCatTF) && ~calcSwitch.spikeTimes && calcSwitch.inducedTrialMagnitudeCorrection
  % subtract the category mean psth from every trial to obtain induced psth
  spikesByCategoryBinnedInduced = spikesByCategoryBinned;
  for cat_i = 1:length(spikesByCategoryBinned)
    for channel_i = 1:length(channelNames)
      for unit_i = 1:length(channelUnitNames{channel_i})
        unitCatMean = mean(spikesByCategoryBinned{cat_i}{channel_i}{unit_i},1);
        unitCatMeanSpikeBkgnd = sum(unitCatMean(1:50)); %todo: set this according to psth pre etc. to get best background estimate
        unitCatMeanBkgndSub = unitCatMean - unitCatMeanSpikeBkgnd;
        unitCatMeanTotalEvokedSpikes = sum(unitCatMeanBkgndSub);
        for trial_i = 1:size(spikesByCategoryBinned{cat_i}{channel_i}{unit_i},1)
          %magnitudeCorrector = mean(spikesByCategoryBinned{cat_i}{channel_i}{unit_i},1)*sum(spikesByCategoryBinned{cat_i}{channel_i}{unit_i});
          trialSpikes = spikesByCategoryBinned{cat_i}{channel_i}{unit_i}(trial_i,:);
          trialSpikesBkgndSub = trialSpikes - mean(trialSpikes(1:50));
          evokedMagnitudeCorrector =  sum(trialSpikesBkgndSub)/unitCatMeanTotalEvokedSpikes;
          spikesByCategoryBinnedInduced{cat_i}{channel_i}{unit_i} = trialSpikesBkgndSub - evokedMagnitudeCorrector*unitCatMeanBkgndSub;
        end
      end
    end
  end
  lfpByCategoryInduced = lfpByCategory;
  for cat_i = 1:length(lfpByCategory)
    for channel_i = 1:length(channelNames)
      channelCatMean = squeeze(mean(lfpByCategory{cat_i}(1,channel_i,:,:),3));
      channelCatMeanMagnitude = mean(abs(channelCatMean));
      for trial_i = 1:size(lfpByCategory{cat_i},3)
        magnitudeCorrector = mean(squeeze(abs(lfpByCategory{cat_i}(1,channel_i,trial_i,:))))/channelCatMeanMagnitude;
        lfpByCategoryInduced{cat_i}(1,channel_i,trial_i,:) = lfpByCategory{cat_i}(1,channel_i,trial_i,:) - magnitudeCorrector*reshape(channelCatMean,1,1,1,[]);
      end
    end
  end
end

%%%% make evoked potential plots by category
times = -psthPre:psthImDur+psthPost; %todo: rename var to eg psthTimes or lfpTimes
for channel_i = 1:length(lfpChannels)
  
  % todo: programatize on category
  % todo: add isolated units
  % too: link axes on subplots
  % todo: psth lineplot with errorbars
  if plotSwitch.faceVnonEvokedPotential || plotSwitch.faceVnonEvokedMuaMultiCh
    faceV = squeeze(mean(lfpByCategory{faceCatNum}(1,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
    nfaceV = squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
  end
 
  % create evoked-psth lineplot, face vs. non, one pane
  if plotSwitch.faceVnonEvokedMuaMultiCh && channel_i == 1
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
  
  
  % face vs. nonface evoked potential plot, one channel
  if plotSwitch.faceVnonEvokedPotential
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
  end
  % contribute to evoked-psth lineplot, face vs. non, one pane
  if plotSwitch.faceVnonEvokedMuaMultiCh
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
  end
  
  
  %lfp category psth
  if plotSwitch.evokedByCategory
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
  end
  
  % make lfp - psth subplot
  if plotSwitch.psthEvokedByCategory
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
  end
end

if plotSwitch.runSummary
  for channel_i = 1:length(channelNames)
    f = figure();
    numSubplots = length(imSpikeCounts{channel_i})-1+4;
    axisHandles = [];
    %juice
    ah = subplot(numSubplots,1,numSubplots);
    ylabel('juice');
    xlabel('time (ms)');
    axisHandles = horzcat(axisHandles,ah);
    hold on
    assert(length(taskDataAll.juiceOnTimes) == length(taskDataAll.juiceOffTimes),'different number of juice on and off triggers; do not have defense');
    assert(min(taskDataAll.juiceOnTimes) < min(taskDataAll.juiceOffTimes),'juice on at beginning of trial; do not have defense');
    for juice_i = 1:length(taskDataAll.juiceOnTimes)
      plot([taskDataAll.juiceOnTimes(juice_i) taskDataAll.juiceOffTimes(juice_i)],[0 0]);
    end
    % fix spot flash
    ah = subplot(numSubplots,1,numSubplots-1);
    ylabel('flash');
    axisHandles = horzcat(axisHandles,ah);
    hold on
    assert(min(taskDataAll.fixSpotFlashStartTimes) < min(taskDataAll.fixSpotFlashEndTimes),'fix spot flash off trigger before on; do not have defense');
    %handle case where flash on at end of run with indexing on end times
    for flash_i = 1:length(taskDataAll.fixSpotFlashEndTimes)
      plot([taskDataAll.fixSpotFlashStartTimes(flash_i) taskDataAll.fixSpotFlashEndTimes(flash_i)],[0 0]);
    end
    %fixation
    ah = subplot(numSubplots,1,numSubplots-2);
    ylabel('fix');
    axisHandles = horzcat(axisHandles,ah);
    hold on
    assert(min(taskDataAll.fixationInTimes) < min(taskDataAll.fixationOutTimes),'fixation out trigger before in; do not have defense');
    %handle case where fix in at end of run
    if length(taskDataAll.fixationInTimes) == length(taskDataAll.fixationOutTimes)+1
      fixationOutTimesTmp = vertcat(taskDataAll.fixationOutTimes,taskDataAll.pictureEndTimes(end));
    else
      fixationOutTimesTmp = taskDataAll.fixationOutTimes;
    end
    for fix_i = 1:length(taskDataAll.fixationInTimes)
      plot([taskDataAll.fixationInTimes(fix_i) fixationOutTimesTmp(fix_i)],[0 0]);
    end
    %spikes and lfp
    for image_i = 1:length(onsetsByImage)
      onsets = onsetsByImage{image_i};
      lfpByTrial = squeeze(lfpByImage{image_i}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
      lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
      lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
      
      for trial_i = 1:length(onsets)
        %trialSubplot = ceil((onsets(trial_i)-firstOnset)/(trialsPerRow*(psthImDur+stimTiming.ISI)));
        %subplot(trialSubplot,1,1);
        
        %lfp
        ah = subplot(numSubplots,1,numSubplots-3);
        if image_i == trial_i && trial_i == 1
          axisHandles = horzcat(axisHandles,ah);
          ylabel('lfp pwr');
        end
        hold on
        plot([onsets(trial_i) onsets(trial_i)+psthImDur], [lfpPowerByTrial(trial_i) lfpPowerByTrial(trial_i)], 'color','blue');
        %spikes
        for unit_i = 1:length(imSpikeCounts{channel_i})-1 %for now, don't show MUA in addition to units and hash
          ah = subplot(numSubplots,1,unit_i);
          if image_i == trial_i && trial_i == 1
            axisHandles = horzcat(axisHandles,ah);
            ylabel(sprintf('fr, u%d',unit_i-1));
          end
          hold on
          plot([onsets(trial_i) onsets(trial_i)+psthImDur], [imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i) imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i)], 'color','blue')
        end
      end
    end
    linkaxes(axisHandles,'x');
    suptitle(channelNames{channel_i});
    drawnow;
    figData = 'none';
    saveFigure(outDir,sprintf('runSummary_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    close(f);
  end
end
  
if plotSwitch.runSummaryImMeanSub
  for channel_i = 1:length(channelNames)
    f = figure();
    numSubplots = length(imSpikeCounts{channel_i})-1+4;
    axisHandles = [];
    %juice
    ah = subplot(numSubplots,1,numSubplots);
    ylabel('juice');
    xlabel('time (ms)');
    axisHandles = horzcat(axisHandles,ah);
    hold on
    assert(length(taskDataAll.juiceOnTimes) == length(taskDataAll.juiceOffTimes),'different number of juice on and off triggers; do not have defense');
    assert(min(taskDataAll.juiceOnTimes) < min(taskDataAll.juiceOffTimes),'juice on at beginning of trial; do not have defense');
    for juice_i = 1:length(taskDataAll.juiceOnTimes)
      plot([taskDataAll.juiceOnTimes(juice_i) taskDataAll.juiceOffTimes(juice_i)],[0 0]);
    end
    % fix spot flash
    ah = subplot(numSubplots,1,numSubplots-1);
    ylabel('flash');
    axisHandles = horzcat(axisHandles,ah);
    hold on
    assert(min(taskDataAll.fixSpotFlashStartTimes) < min(taskDataAll.fixSpotFlashEndTimes),'fix spot flash off trigger before on; do not have defense');
    %handle case where flash on at end of run with indexing on end times
    for flash_i = 1:length(taskDataAll.fixSpotFlashEndTimes)
      plot([taskDataAll.fixSpotFlashStartTimes(flash_i) taskDataAll.fixSpotFlashEndTimes(flash_i)],[0 0]);
    end
    %fixation
    ah = subplot(numSubplots,1,numSubplots-2);
    ylabel('fix');
    axisHandles = horzcat(axisHandles,ah);
    hold on
    assert(min(taskDataAll.fixationInTimes) < min(taskDataAll.fixationOutTimes),'fixation out trigger before in; do not have defense');
    %handle case where fix in at end of run
    if length(taskDataAll.fixationInTimes) == length(taskDataAll.fixationOutTimes)+1
      fixationOutTimesTmp = vertcat(taskDataAll.fixationOutTimes,taskDataAll.pictureEndTimes(end));
    end
    for fix_i = 1:length(taskDataAll.fixationInTimes)
      plot([taskDataAll.fixationInTimes(fix_i) fixationOutTimesTmp(fix_i)],[0 0]);
    end
    %spikes and lfp
    for image_i = 1:length(onsetsByImage)
      onsets = onsetsByImage{image_i};
      lfpByTrial = squeeze(lfpByImage{image_i}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
      lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
      lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
      lfpPowerByTrialImMeanSub = lfpPowerByTrial - mean(lfpPowerByTrial);
      
      for trial_i = 1:length(onsets)
        %trialSubplot = ceil((onsets(trial_i)-firstOnset)/(trialsPerRow*(psthImDur+stimTiming.ISI)));
        %subplot(trialSubplot,1,1);
        
        %lfp
        ah = subplot(numSubplots,1,numSubplots-3);
        if image_i == trial_i && trial_i == 1
          axisHandles = horzcat(axisHandles,ah);
          ylabel('lfp pwr fluct.');
        end
        hold on
        plot([onsets(trial_i) onsets(trial_i)+psthImDur], [lfpPowerByTrialImMeanSub(trial_i) lfpPowerByTrialImMeanSub(trial_i)], 'color','blue');
        %spikes
        for unit_i = 1:length(imSpikeCounts{channel_i})-1 %for now, don't show MUA in addition to units and hash
          ah = subplot(numSubplots,1,unit_i);
          if image_i == trial_i && trial_i == 1
            axisHandles = horzcat(axisHandles,ah);
            ylabel(sprintf('fr fluct., u%d',unit_i-1));
          end
          hold on
          imageUnitMeanFr = mean(imSpikeCounts{channel_i}{unit_i}{image_i}.counts);
          plot([onsets(trial_i) onsets(trial_i)+psthImDur], [imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i) imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i)] - imageUnitMeanFr, 'color','blue')
        end
      end
    end
    linkaxes(axisHandles,'x');
    suptitle(sprintf('LFP and FR trial fluctuations, %s',channelNames{channel_i}));
    drawnow;
    figData = 'none';
    saveFigure(outDir,sprintf('runSummary_fluct_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    close(f);
  end
end
  
if plotSwitch.runSummaryImMeanSubDiv
  for channel_i = 1:length(channelNames)
    f = figure();
    numSubplots = length(imSpikeCounts{channel_i})-1+4;
    axisHandles = [];
    %juice
    ah = subplot(numSubplots,1,numSubplots);
    ylabel('juice');
    xlabel('time (ms)');
    axisHandles = horzcat(axisHandles,ah);
    hold on
    assert(length(taskDataAll.juiceOnTimes) == length(taskDataAll.juiceOffTimes),'different number of juice on and off triggers; do not have defense');
    assert(min(taskDataAll.juiceOnTimes) < min(taskDataAll.juiceOffTimes),'juice on at beginning of trial; do not have defense');
    for juice_i = 1:length(taskDataAll.juiceOnTimes)
      plot([taskDataAll.juiceOnTimes(juice_i) taskDataAll.juiceOffTimes(juice_i)],[0 0]);
    end
    % fix spot flash
    ah = subplot(numSubplots,1,numSubplots-1);
    ylabel('flash');
    axisHandles = horzcat(axisHandles,ah);
    hold on
    assert(min(taskDataAll.fixSpotFlashStartTimes) < min(taskDataAll.fixSpotFlashEndTimes),'fix spot flash off trigger before on; do not have defense');
    %handle case where flash on at end of run with indexing on end times
    for flash_i = 1:length(taskDataAll.fixSpotFlashEndTimes)
      plot([taskDataAll.fixSpotFlashStartTimes(flash_i) taskDataAll.fixSpotFlashEndTimes(flash_i)],[0 0]);
    end
    %fixation
    ah = subplot(numSubplots,1,numSubplots-2);
    ylabel('fix');
    axisHandles = horzcat(axisHandles,ah);
    hold on
    assert(min(taskDataAll.fixationInTimes) < min(taskDataAll.fixationOutTimes),'fixation out trigger before in; do not have defense');
    %handle case where fix in at end of run
    if length(taskDataAll.fixationInTimes) == length(taskDataAll.fixationOutTimes)+1
      fixationOutTimesTmp = vertcat(taskDataAll.fixationOutTimes,taskDataAll.pictureEndTimes(end));
    end
    for fix_i = 1:length(taskDataAll.fixationInTimes)
      plot([taskDataAll.fixationInTimes(fix_i) fixationOutTimesTmp(fix_i)],[0 0]);
    end
    %spikes and lfp
    for image_i = 1:length(onsetsByImage)
      onsets = onsetsByImage{image_i};
      lfpByTrial = squeeze(lfpByImage{image_i}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
      lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
      lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
      lfpPowerByTrialImMeanSubDiv = (lfpPowerByTrial - mean(lfpPowerByTrial))/mean(lfpPowerByTrial);
      
      for trial_i = 1:length(onsets)
        %trialSubplot = ceil((onsets(trial_i)-firstOnset)/(trialsPerRow*(psthImDur+stimTiming.ISI)));
        %subplot(trialSubplot,1,1);
        
        %lfp
        ah = subplot(numSubplots,1,numSubplots-3);
        if image_i == trial_i && trial_i == 1
          axisHandles = horzcat(axisHandles,ah);
          ylabel('lfp pwr fluct.');
        end
        hold on
        plot([onsets(trial_i) onsets(trial_i)+psthImDur], [lfpPowerByTrialImMeanSubDiv(trial_i) lfpPowerByTrialImMeanSubDiv(trial_i)], 'color','blue');
        %spikes
        for unit_i = 1:length(imSpikeCounts{channel_i})-1 %for now, don't show MUA in addition to units and hash
          ah = subplot(numSubplots,1,unit_i);
          if image_i == trial_i && trial_i == 1
            axisHandles = horzcat(axisHandles,ah);
            ylabel(sprintf('fr fluct., u%d',unit_i-1));
          end
          hold on
          imageUnitMeanFr = mean(imSpikeCounts{channel_i}{unit_i}{image_i}.counts);
          plot([onsets(trial_i) onsets(trial_i)+psthImDur], ([imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i) imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i)] - imageUnitMeanFr)/imageUnitMeanFr, 'color','blue');
        end
      end
    end
    linkaxes(axisHandles,'x');
    suptitle(sprintf('LFP and FR trial fractional fluctuations, %s',channelNames{channel_i}));
    drawnow;
    figData = 'none';
    saveFigure(outDir,sprintf('runSummary_fracfluct_%s_Run%s',channelNames{channel_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    close(f);
  end
end
  
  % make lfp power-fr scatter plot
if plotSwitch.lfpPowerMuaScatterAll
  for channel_i = 1:length(channelNames)
    for unit_i = 1:length(channelUnitNames{channel_i})
      figure();
      hold on;
      for cat_i = 1:length(categoryListSlim)
        lfpByTrial = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
        lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
        lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
        scatter(lfpPowerByTrial,catSpikeCounts{channel_i}{unit_i}{categorySlimInds(cat_i)}.counts,36,colors(cat_i),'filled')
      end
      xlabel('total lfp power (uV)');
      ylabel('firing rate (Hz)');
      title(sprintf('LFP power vs Firing rate, %s %s, %d ms - %d ms',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}, frCalcOn, frCalcOff))
      saveFigure(outDir,sprintf('scatter_lfpPwrVsFr_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
end
  
if plotSwitch.lfpPeakToPeakMuaScatterAll
  for channel_i = 1:length(channelNames)
    for unit_i = 1:length(channelUnitNames{channel_i})
      figure();
      hold on;
      for cat_i = 1:length(categoryListSlim)
        lfpByTrial = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
        scatter(max(lfpByTrial,[],2)-min(lfpByTrial,[],2),catSpikeCounts{channel_i}{unit_i}{categorySlimInds(cat_i)}.counts,36,colors(cat_i),'filled')
      end
      xlabel('Peak-to-Trough LFP Amplitude (uV)');
      ylabel('firing rate (Hz)');
      title(sprintf('LFP Amplitude vs Firing rate, %s %s, %d ms - %d ms',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}, frCalcOn, frCalcOff))
      saveFigure(outDir,sprintf('scatter_lfpP2PVsFr_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
end
  
if plotSwitch.lfpLatencyMuaLatency
  for channel_i = 1:length(channelNames)
  end
end
  
% make lfp power-fr scatter plot, early
if plotSwitch.lfpPowerMuaScatterAllEarly
  for channel_i = 1:length(channelNames)
    for unit_i = 1:length(channelUnitNames{channel_i})
      figure();
      hold on;
      for cat_i = 1:length(categoryListSlim)
        lfpByTrial = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOnEarly:lfpPaddedBy+1+psthPre+frCalcOffEarly));
        lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
        lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
        scatter(lfpPowerByTrial,catSpikeCountsEarly{channel_i}{unit_i}{categorySlimInds(cat_i)}.counts,36,colors(cat_i),'filled')
      end
      xlabel('total lfp power (uV)');
      ylabel('firing rate (Hz)');
      title(sprintf('LFP power vs Firing rate, %s %s, %d ms - %d ms',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}, frCalcOnEarly, frCalcOffEarly))
      saveFigure(outDir,sprintf('scatter_lfpPwrVsFr_early_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
end
  
if plotSwitch.lfpPeakToPeakMuaScatterAllEarly
  for channel_i = 1:length(channelNames)
    for unit_i = 1:length(channelUnitNames{channel_i})
      figure();
      hold on;
      for cat_i = 1:length(categoryListSlim)
        lfpByTrial = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOnEarly:lfpPaddedBy+1+psthPre+frCalcOffEarly));
        scatter(max(lfpByTrial,[],2)-min(lfpByTrial,[],2),catSpikeCountsEarly{channel_i}{unit_i}{categorySlimInds(cat_i)}.counts,36,colors(cat_i),'filled')
      end
      xlabel('Peak-to-Trough LFP Amplitude (uV)');
      ylabel('firing rate (Hz)');
      title(sprintf('LFP Amplitude vs Firing rate, %s %s, %d ms - %d ms',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}, frCalcOnEarly, frCalcOffEarly))
      saveFigure(outDir,sprintf('scatter_lfpP2PVsFr_early_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
  
if plotSwitch.lfpLatencyMuaLatencyEarly
  for channel_i = 1:length(channelNames)
  end
end
  
% make lfp power-fr scatter plot, late
if plotSwitch.lfpPowerMuaScatterAllLate
  for channel_i = 1:length(channelNames)
    for unit_i = 1:length(channelUnitNames{channel_i})
      figure();
      hold on;
      for cat_i = 1:length(categoryListSlim)
        lfpByTrial = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOnLate:lfpPaddedBy+1+psthPre+frCalcOffLate));
        lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
        lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
        scatter(lfpPowerByTrial,catSpikeCountsLate{channel_i}{unit_i}{categorySlimInds(cat_i)}.counts,36,colors(cat_i),'filled')
      end
      xlabel('total lfp power (uV)');
      ylabel('firing rate (Hz)');
      title(sprintf('LFP power vs Firing rate, %s %s, %d ms - %d ms',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}, frCalcOnLate, frCalcOffLate))
      saveFigure(outDir,sprintf('scatter_lfpPwrVsFr_late_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
end
  
if plotSwitch.lfpPeakToPeakMuaScatterAllLate
  for channel_i = 1:length(channelNames)
    for unit_i = 1:length(channelUnitNames{channel_i})
      figure();
      hold on;
      for cat_i = 1:length(categoryListSlim)
        lfpByTrial = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOnLate:lfpPaddedBy+1+psthPre+frCalcOffLate));
        scatter(max(lfpByTrial,[],2)-min(lfpByTrial,[],2),catSpikeCountsLate{channel_i}{unit_i}{categorySlimInds(cat_i)}.counts,36,colors(cat_i),'filled')
      end
      xlabel('Peak-to-Trough LFP Amplitude (uV)');
      ylabel('firing rate (Hz)');
      title(sprintf('LFP Amplitude vs Firing rate, %s %s, %d ms - %d ms',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}, frCalcOnLate, frCalcOffLate))
      saveFigure(outDir,sprintf('scatter_lfpP2PVsFr_late_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
end
  
if plotSwitch.lfpLatencyMuaLatencyLate
  for channel_i = 1:length(channelNames)
  end
end
  
% lfp power vs lfp power, across channels
if plotSwitch.lfpPowerAcrossChannels && channel_i < length(lfpChannels)
  for channel_i = 1:length(channelNames)
    for channel2_i = channel_i+1:length(lfpChannels)
      figure();
      hold on;
      for cat_i = 1:length(categoryListSlim)
        lfpByTrialCh1 = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
        lfpMeanSubByTrialCh1 = lfpByTrialCh1 - mean(lfpByTrialCh1,2);
        lfpPowerByTrialCh1 = sqrt(sum(lfpMeanSubByTrialCh1.^2,2));
        lfpByTrialCh2 = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel2_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
        lfpMeanSubByTrialCh2 = lfpByTrialCh2 - mean(lfpByTrialCh2,2);
        lfpPowerByTrialCh2 = sqrt(sum(lfpMeanSubByTrialCh2.^2,2));
        scatter(lfpPowerByTrialCh1,lfpPowerByTrialCh2,36,colors(cat_i),'filled')
      end
      xlabel(sprintf('total lfp power, %s (uV)',channelNames{channel_i}));
      ylabel(sprintf('total lfp power, %s (uV)',channelNames{channel2_i}));
      title(sprintf('LFP Power, %s vs. %s, %d ms - %d ms',channelNames{channel_i},channelNames{channel2_i}, frCalcOn, frCalcOff))
      saveFigure(outDir,sprintf('scatter_lfpPwrVSlfpPwr_%s_%s_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
end
  
% lfp p2p amplitude vs lfp p2p amplitude, across channels
if plotSwitch.lfpPeakToPeakAcrossChannels && channel_i < length(lfpChannels)
  for channel_i = 1:length(channelNames)
    for channel2_i = channel_i+1:length(lfpChannels)
      figure();
      hold on;
      for cat_i = 1:length(categoryListSlim)
        lfpByTrialCh1 = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
        lfpByTrialCh2 = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel2_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
        scatter(max(lfpByTrialCh1,[],2)-min(lfpByTrialCh1,[],2),max(lfpByTrialCh2,[],2)-min(lfpByTrialCh2,[],2),36,colors(cat_i),'filled')
      end
      xlabel(sprintf('LFP Amplitude, %s (uV)',channelNames{channel_i}));
      ylabel(sprintf('LFP Amplitude, %s (uV)',channelNames{channel2_i}));
      title(sprintf('LFP Amplitude, %s vs. %s, %d ms - %d ms',channelNames{channel_i},channelNames{channel2_i}, frCalcOn, frCalcOff))
      saveFigure(outDir,sprintf('scatter_lfpP2PVSlfpP2P_%s_%s_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
end
  
% lfp latency vs lfp latency, across channels, defined as max dot product with mean
if plotSwitch.lfpLatencyShiftAcrossChannels && channel_i < length(lfpChannels)
  for channel_i = 1:length(channelNames)
    for channel2_i = channel_i+1:length(lfpChannels)
      figure();
      hold on;
      for cat_i = 1:length(categoryListSlim)
        lfpByTrialCh1 = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
        lfpMeanSubByTrialCh1 = lfpByTrialCh1 - mean(lfpByTrialCh1,2);
        lfpTrialAveCh1 = mean(lfpMeanSubByTrialCh1,1);
        
        lfpByTrialCh2 = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel2_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
        lfpMeanSubByTrialCh2 = lfpByTrialCh2 - mean(lfpByTrialCh2,2);
        lfpTrialAveCh2 = mean(lfpMeanSubByTrialCh2,1);
        shiftsCh1 = zeros(size(lfpByTrialCh1,1),1);
        shiftsCh2 = zeros(size(lfpByTrialCh1,1),1);
        for trial_i = 1:size(lfpByTrialCh1,1)
          bestShiftCh1 = -50;
          bestMatchCh1 = 0;
          for shift = -50:50 %magic number; todo replace with variable
            shiftedTrialLfp = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,trial_i,lfpPaddedBy+1+psthPre+frCalcOn+shift:lfpPaddedBy+1+psthPre+frCalcOff+shift));
            shiftedTrialLfp = shiftedTrialLfp - mean(shiftedTrialLfp);
            match = lfpTrialAveCh1*shiftedTrialLfp;
            if match > bestMatchCh1
              bestMatchCh1 = match;
              bestShiftCh1 = shift;
            end
            shiftsCh1(trial_i) = bestShiftCh1;
          end
          bestShiftCh2 = -50;
          bestMatchCh2 = 0;
          for shift = -50:50 %magic number; todo replace with variable
            shiftedTrialLfp = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel2_i,trial_i,lfpPaddedBy+1+psthPre+frCalcOn+shift:lfpPaddedBy+1+psthPre+frCalcOff+shift));
            shiftedTrialLfp = shiftedTrialLfp - mean(shiftedTrialLfp);
            match = lfpTrialAveCh2*shiftedTrialLfp;
            if match > bestMatchCh2
              bestMatchCh2 = match;
              bestShiftCh2 = shift;
            end
            shiftsCh2(trial_i) = bestShiftCh2;
          end
        end
        scatter(shiftsCh1,shiftsCh2,36,colors(cat_i),'filled')
      end
      xlabel(sprintf('LFP latency shift from category mean, %s (ms)',channelNames{channel_i}));
      ylabel(sprintf('LFP latency shift from category mean, %s (ms)',channelNames{channel2_i}));
      title(sprintf('LFP Latency shift from category mean, %s vs. %s, %d ms - %d ms',channelNames{channel_i},channelNames{channel2_i}, frCalcOn, frCalcOff))
      saveFigure(outDir,sprintf('scatter_lfpShiftVSlfpShift_%s_%s_Run%s',channelNames{channel_i},channelNames{channel2_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
    end
  end
end
  
% other to-do: one epoch for each peak of trial ave evoked potential, ideally set automatically
% todo: across-peak latency within and across channels
% todo: spike burstiness vs. lfp bumpiness; can do as p2p evoked/psth; or (low freq power)/(high freq power)

calcSwitchNames = {'evoked', 'induced'};
tfCalcSwitchTitleSuffixes = {'',', induced'}; % appended to titles
tfCalcSwitchFnameSuffixes = {'','_induced'}; % appended to filenames
calcSwitches = [calcSwitch.evokedSpectra, calcSwitch.inducedSpectra];
for calc_i = 1:length(calcSwitches)
  if calcSwitches(calc_i)
    if strcmp(calcSwitchNames{calc_i},'induced')
      if calcSwitch.spikeTimes
        disp('induced TF not yet implemented for spike times; change to spike bins using calcSwitch.spikeTimes = 0');
        continue
      end
      disp('switching spikes and lfps to induced for spectrum computation');
      spikesByCategoryBinnedEvokedTmp = spikesByCategoryBinned;
      spikesByCategoryBinned = spikesByCategoryBinnedInduced;
      lfpByCategoryEvokedTmp = lfpByCategory;
      lfpByCategory = lfpByCategoryInduced;
    end
    % single trial evoked lfps; better as multi-channel subplot?
    % todo: programatize on category
    if plotSwitch.singleTrialLfpFaceVnon
      for channel_i = 1:length(channelNames)
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
        title(sprintf('single trial face LFPs, %s%s', channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
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
        title(sprintf('single trial nonface LFPs, %s%s', channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
        xlabel('time after stimulus (ms)', 'FontSize',18);
        ylabel('voltage (uV)', 'FontSize',18);
        xlim([min(times) max(times)]);
        figData.ax21.y = ydata;
        figData.ax21.x = times;
        saveFigure(outDir,sprintf('Evoked_singleTrials_%s%s_Run%s',channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      end
    end
      
      
    % lfp spectra, face v non
    if plotSwitch.lfpSpectraFaceVnon
      % todo: programtize on category
      for channel_i = 1:length(channelNames)
        [faceS,faceF, faceE] = mtspectrumc(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))', chr_params);
        [nfaceS,nfaceF, nfaceE] = mtspectrumc(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))', chr_params);

        figure();
        plot(1000*faceF,log(faceS),'linewidth',3,'marker','o');
        hold on
        plot(1000*nfaceF,log(nfaceS),'linewidth',3,'color','r','marker','o');
        hold off
        title(sprintf('%s evoked power spectrum%s',channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
        xlabel('frequency (Hz)', 'FontSize',18);
        ylabel('voltage, log(uV)', 'FontSize',18);
        legend('face','non');
        set(gca,'fontsize',18);
        clear figData
        figData.y = vertcat(faceS,nfaceS);
        figData.x = vertcat(faceF,nfaceF);
        saveFigure(outDir,sprintf('spectrum_faceVnon_%s%s_Run%s',channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );

        figure();
        specModel = fitlm(log(1000*faceF(6:end)), log(.5*(faceS(6:end)+nfaceS(6:end))));
        m = specModel.Coefficients.Estimate(2);
        y0 = specModel.Coefficients.Estimate(1);
        h1 = loglog(1000*faceF,exp(m*log(1000*faceF) + y0), 'k--');
        hold on
        h2 = loglog(1000*faceF,faceS,'linewidth',3,'marker','o');
        h3 = loglog(1000*nfaceF,nfaceS,'linewidth',3,'color','r','marker','o');
        hold off
        title(sprintf('%s evoked power spectrum%s',channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
        xlim([1000*min(faceF) 1000*max(faceF)]); %todo: fix error
        xlabel('frequency (Hz)', 'FontSize',18);
        ylabel('voltage (uV)', 'FontSize',18);
        legend([h2;h3;h1],{'face','nonface',sprintf('fit,m = %.3f',m)});
        set(gca,'fontsize',18);
        clear figData
        figData.y = vertcat(faceS,nfaceS);
        figData.x = vertcat(faceF,nfaceF);
        saveFigure(outDir,sprintf('spectrum_log_faceVnon_%s%s_Run%s',channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
      end
    end

    if plotSwitch.spikeSpectraFacevNon
      for channel_i = 1:length(channelNames)
        for unit_i = 1:length(channelUnitNames{channel_i})
          if length(channelUnitNames{channel_i}) == 2 && unit_i == 1 %skip unsorted if no isolated units defined
            continue
          end
          if calcSwitch.spikeTimes
            [faceS,faceF,faceR,faceE] = mtspectrumpt(spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i}, chr_params);
            [nfaceS,nfaceF,nfaceR nfaceE] = mtspectrumpt(spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i}, chr_params);
          else
            [faceS,faceF,nfaceR,faceE] = mtspectrumpb(spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}', chr_params);
            [nfaceS,nface,nfaceRF,nfaceE] = mtspectrumpb(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}', chr_params);
          end
          figure();
          plot(1000*faceF,log(faceS),'linewidth',3,'marker','o');
          hold on
          plot(1000*nfaceF,log(nfaceS),'linewidth',3,'color','r','marker','o');
          hold off
          title(sprintf('%s %s evoked power spectrum%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
          xlabel('frequency (Hz)', 'FontSize',18);
          ylabel('spike power, log(spks/s)', 'FontSize',18);
          legend('face','non');
          set(gca,'fontsize',18);
          clear figData
          figData.y = vertcat(faceS,nfaceS);
          figData.x = vertcat(faceF,nfaceF);
          saveFigure(outDir,sprintf('spectrum_faceVnon_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );

          figure();
          specModel = fitlm(log(1000*faceF(6:end)), log(.5*(faceS(6:end)+nfaceS(6:end))));
          m = specModel.Coefficients.Estimate(2);
          y0 = specModel.Coefficients.Estimate(1);
          h1 = loglog(1000*faceF,exp(m*log(1000*faceF) + y0), 'k--');
          hold on
          h2 = loglog(1000*faceF,faceS,'linewidth',3,'marker','o');
          h3 = loglog(1000*nfaceF,nfaceS,'linewidth',3,'color','r','marker','o');
          hold off
          title(sprintf('%s %s evoked power spectrum%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
          xlim([1000*min(faceF) 1000*max(faceF)]); %todo: fix error
          xlabel('frequency (Hz)', 'FontSize',18);
          ylabel('power log(spks/s)', 'FontSize',18);
          legend([h2;h3;h1],{'face','nonface',sprintf('fit,m = %.3f',m)});
          set(gca,'fontsize',18);
          clear figData
          figData.y = vertcat(faceS,nfaceS);
          figData.x = vertcat(faceF,nfaceF);
          saveFigure(outDir,sprintf('spectrum_log_faceVnon_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        end
      end
    end
    % clean up temporary variables and restore stable variables
    if strcmp(calcSwitchNames{calc_i},'induced')
      spikesByCategoryBinned = spikesByCategoryBinnedEvokedTmp;
      lfpByCategory =  lfpByCategoryEvokedTmp;
    end
  end
end
  
% image-wise time-frequency plots, spikes and lfp
% todo: put in option for just preferred image
tfCalcSwitchNames = {'evokedImageTF', 'inducedImageTF'};
tfCalcSwitchTitleSuffixes = {'',', induced'}; % appended to titles
tfCalcSwitchFnameSuffixes = {'','_induced'}; % appended to filenames
tfCalcSwitches = [calcSwitch.evokedImageTF, calcSwitch.inducedImageTF];
for calc_i = 1:length(tfCalcSwitches)
  if tfCalcSwitches(calc_i)
    if strcmp(tfCalcSwitchNames{calc_i},'inducedImageTF')
      if calcSwitch.spikeTimes
        disp('induced TF not yet implemented for spike times; change to spike bins using calcSwitch.spikeTimes = 0');
        continue
      end
      disp('switching spikes and lfps to induced for imagewise tf computation');
      spikesByImageBinnedEvokedTmp = spikesByImageBinned;
      spikesByImageBinned = spikesByImageBinnedInduced;
      lfpByImageEvokedTmp = lfpByImage;
      lfpByImage = lfpByImageInduced;
    end
    if plotSwitch.spikeSpectraTfByImage
      for channel_i = 1:length(channelNames)
        for image_i = 1:length(pictureLabels)
          % todo: put in dB conversion and specgramrowave option
          for unit_i =1:length(channelUnitNames{channel_i})
            if length(channelUnitNames{channel_i}) == 1 && unit_i ==1 %if no isolated spike, don't separate unsorted and MUA
              continue
            end
            if calcSwitch.spikeTimes
              [S,t,f,R]=mtspecgrampt(spikesByImageForTF{image_i}{channel_i}{unit_i},movingWin,chr_params); %last optional param not used: fscorr
            else
              [S,t,f,R]=mtspecgrampb(spikesByImageBinned{image_i}{channel_i}{unit_i}',movingWin,chr_params); %last optional param not used: fscorr
            end
            t = t - lfpAlignParams.msPreAlign;
            f = 1000*f;
            fh = figure();
            imagesc(t,f,S'); %TODO: fix to match cat version, which is correct!!
            set(gca,'Ydir','normal'); % also need 'Yscale,'log'??
            xlabel('Time'); % todo: units
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'Power'); % todo: check units
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s %s Time-Frequency, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i},tfCalcSwitchTitleSuffixes{calc_i}));
            clear FigData
            figData.x = t;
            figData.y = f;
            figData.z = S';
            drawnow;
            saveFigure(outDir,sprintf('TF_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},pictureLabels{image_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);
          end
        end
      end
    end
    if plotSwitch.lfpSpectraTfByImage
      for channel_i = 1:length(channelNames)
        [S,t,f]=mtspecgramc(squeeze(lfpByImage{image_i}(1,channel_i,:,:))',movingWin,chr_params);
        t = t - lfpAlignParams.msPreAlign;
        f = 1000*f;
        fh = figure();
        imagesc(t,f,S');
        set(gca,'Ydir','normal'); % also need 'Yscale,'log'??
        xlabel('Time'); % todo: units
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'Power'); % todo: check units
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s LFP Time-Frequency, %s%s',channelNames{channel_i},pictureLabels{image_i},tfCalcSwitchTitleSuffixes{calc_i}));
        clear figData
        figData.x = t;
        figData.y = f;
        figData.z = S;
        drawnow;
        saveFigure(outDir,sprintf('TF_LFP_%s_%s%s_Run%s',channelNames{channel_i},pictureLabels{image_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        close(fh);
      end
    end
    if strcmp(tfCalcSwitchNames{calc_i},'inducedImageTF')
      spikesByImageBinned = spikesByImageBinnedEvokedTmp;
      lfpByImage = lfpByImageEvokedTmp;
    end
  end
end
  
  
% category-wise time-frequency plots, spikes and lfp
tfCalcSwitchNames = {'evokedCatTF', 'inducedCatTF'};
tfCalcSwitchTitleSuffixes = {'',', induced'}; % appended to titles
tfCalcSwitchFnameSuffixes = {'','_induced'}; % appended to filenames
tfCalcSwitches = [calcSwitch.evokedCatTF, calcSwitch.inducedCatTF];
for calc_i = 1:length(tfCalcSwitches)
  if tfCalcSwitches(calc_i)
    if strcmp(tfCalcSwitchNames{calc_i},'inducedCatTF')
      if calcSwitch.spikeTimes
        disp('induced TF not yet implemented for spike times; change to spike bins using calcSwitch.spikeTimes = 0');
        continue
      end
      disp('switching spikes and lfps to induced for categorywise coupling computation');
      spikesByCategoryBinnedEvokedTmp = spikesByCategoryBinned;
      spikesByCategoryBinned = spikesByCategoryBinnedInduced;
      lfpByCategoryEvokedTmp = lfpByCategory;
      lfpByCategory = lfpByCategoryInduced;
    end
    for channel_i = 1:length(channelNames)
      if channel_i == 1
        lfpTfFig = figure();
        muaTfFig = figure();
      end
      % where I stopped
      for cat_i = 1:length(categoryList) %todo: change this to programatize on category
        if cat_i ~= faceCatNum && cat_i ~= nonfaceCatNum
          continue
        end
        % spikes
        for unit_i = 1:length(channelUnitNames{channel_i})
          if length(channelUnitNames{channel_i}) == 1 && unit_i == 1
            continue
          end
          if calcSwitch.spikeTimes
            [S,t,f,R]=mtspecgrampt(spikesByCategoryForTF{cat_i}{channel_i}{unit_i},movingWin,chr_params); %last optional param not used: fscorr
          else
            [S,t,f,R]=mtspecgrampb(spikesByCategoryBinned{cat_i}{channel_i}{unit_i}',movingWin,chr_params); %last optional param not used: fscorr
          end
          t = t - lfpAlignParams.msPreAlign;
          f = 1000*f;
          totalPower = sum(S,2)';
          fh = figure();
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
          title(sprintf('%s %s Time-Frequency, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},categoryList{cat_i},tfCalcSwitchTitleSuffixes{calc_i}));
          clear figData
          figData.x = t;
          figData.y = f;
          figData.z = S;
          drawnow;
          saveFigure(outDir,sprintf('TF_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},categoryList{cat_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
          close(fh);
          
%           % autocorrellogram, TF TODO finish
%           if ~calcSwitch.spikeTimes
%             autocorrel = zeros(201,1);
%             spikeMat = spikesByCategory{cat_i}{channel_i}{unit_i};
%             shift_i = 0;
%             for shift = 100:-1:1
%               shift_i = shift_i + 1;
%               autocorrel(shift_i) = mean(mean(spikeMat(:,1:end-shift).*spikeMat(:,shift:end)));
%             end
%             shift_i = shift_i + 1;
%             autocorrel(shift_i) = 1;
%             for shift = 1:100
%               shift_i = shift_i + 1;
%               autocorrel(shift_i) = mean(mean(spikeMat(:,1:end-shift).*spikeMat(:,shift:end)));
%             end
%           end
%           figure();
%           plot(-100:1:100,autocorrel);
%           % todo: autocorrellogram with spike times
          
          
          % contribute mua to shared figure
          if unit_i == length(channelUnitNames{channel_i})
            figure(muaTfFig);
            if cat_i == faceCatNum
              subplot(3,2,2*channel_i-1);
              forTitle = sprintf('Face MUA TF, %s%s',channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i});
            else
              subplot(3,2,2*channel_i);
              forTitle = sprintf('Nonface MUA TF, %s%s',channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i});
            end
            imagesc(t,f,S'); axis xy; c = colorbar(); ylabel(c,'Power (dB)'); %todo: fix unit when row-normalizing
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s LFP Time-Frequency, %s%s',channelNames{channel_i},pictureLabels{image_i},tfCalcSwitchTitleSuffixes{calc_i}));
            xlabel('Time (ms)'); 
            ylabel('Frequency (Hz)');
            title(forTitle);
          end
        end
        % lfp
        [S,t,f]=mtspecgramc(squeeze(lfpByCategory{cat_i}(1,channel_i,:,:))',movingWin,chr_params);
        t = t - lfpAlignParams.msPreAlign;
        f = 1000*f;
        totalPower = sum(S,2)';

        fh = figure();
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
        title(sprintf('%s LFP Time-Frequency, %s%s',channelNames{channel_i},categoryList{cat_i},tfCalcSwitchTitleSuffixes{calc_i}));
        clear figData
        figData.x = t;
        figData.y = f;
        figData.z = S';
        drawnow;
        saveFigure(outDir,sprintf('TF_%s_LFP_%s%s_Run%s',channelNames{channel_i},categoryList{cat_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
        close(fh);
        % contribute to shared figure
        figure(lfpTfFig);
        if cat_i == faceCatNum
          subplot(3,2,2*channel_i-1);
          forTitle = sprintf('Face LFP TF, %s%s',channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i});
        else
          subplot(3,2,2*channel_i);
          forTitle = sprintf('Nonface LFP TF, %s%s',channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i});
        end
        imagesc(t,f,S'); axis xy; c = colorbar(); ylabel(c,'Power (dB)'); %todo: fix unit when row-normalizing
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        xlabel('Time (ms)'); 
        ylabel('Frequency (Hz)');
        title(forTitle);
      end

      drawnow;
      % todo: convert to saveFigure and build its fig data as we go
      if saveFig && channel_i == length(lfpChannels) && cat_i == max(faceCatNum, nonfaceCatNum)
        figure(muaTfFig);
        savefig(strcat(outDir,sprintf('TF_MUA_%s_%s%s.fig',channelNames{channel_i},categoryList{cat_i},tfCalcSwitchFnameSuffixes{calc_i})));
        export_fig([outDir sprintf('TF_MUA_%s_%s%s.png',channelNames{channel_i},categoryList{cat_i},tfCalcSwitchFnameSuffixes{calc_i})],'-m1.2','-transparent','-opengl');
        close(muaTfFig);
        figure(lfpTfFig);
        savefig(lfpTfFig,strcat(outDir,sprintf('TF_LFP_%s_%s%s.fig',channelNames{channel_i},categoryList{cat_i},tfCalcSwitchFnameSuffixes{calc_i})));
        export_fig([outDir sprintf('TF_LFP_%s_%s%s.png',channelNames{channel_i},categoryList{cat_i},tfCalcSwitchFnameSuffixes{calc_i})],'-m1.2','-transparent','-opengl');
        close(lfpTfFig);
      end
    
  


      %%% coupling across channels and modalities (units/unsorted/mua, fields)
      if calcSwitch.crossTF
        % face vs. nonface spike field coherence, within and across channels
        % channel_i spike -- channel_i field
        for unit_i = 1:length(channelUnitNames{channel_i})
          if length(channelUnitNames{channel_i}) == 2 && unit_i == 1
            continue
          end
          
          %todo: STA LFP
          
          if calcSwitch.useJacknife
            chr_params.err = [2 .05];
            if calcSwitch.spikeTimes
              [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},chr_params);
              [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},chr_params);
            else
              [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpb(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',chr_params);
              [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpb(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}',chr_params);
            end
            errs = vertcat(faceErrs,nfaceErrs);
          else
            if calcSwitch.spikeTimes
              [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},chr_params);
              [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},chr_params);
            else
              [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpb(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',chr_params);
              [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpb(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}',chr_params);
            end
            errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
          end
          % note: this correction maintains consistency with X-Y coherence, where X comes first in the call to coherency, but so we can say X=spike,
          % Y=field, while chronux requires us to put field first in the function call
          phiface = -1*phiface;
          phinface = -1*phinface;

          fh = figure();
          mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
          legend('face', 'nonface');
          xlabel('frequency (Hz)');
          ylabel('coherency');
          title(sprintf('%s %s - %s field coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
          clear figData
          figData.x = vertcat(fface, fnface);
          figData.y = [Cface, Cnface];
          drawnow;
          saveFigure(outDir,sprintf('coh_%s_%s_%s_LFP_FACEvsNON%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
          close(fh);
        end
        if calcSwitch.meanEvokedTF  
          for unit_i = 1:length(channelUnitNames{channel_i})
            if length(channelUnitNames{channel_i}) == 2 && unit_i == 1
              continue
            end
            if calcSwitch.spikeTimes
              [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(mean(lfpByCategory{faceCatNum}(1,channel_i,:,:),3)),...
                allSpikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},chr_params);
              [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:),3)),...
                allSpikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},chr_params);
            else
              [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpb(squeeze(mean(lfpByCategory{faceCatNum}(1,channel_i,:,:),3)),...
                mean(spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i},1)',chr_params);
              [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpb(squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:),3)),...
                mean(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i},1)',chr_params);
            end
            errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));

            phiface = -1*phiface;
            phinface = -1*phinface;

            fh = figure();
            mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
            legend('face', 'nonface');
            xlabel('frequency (Hz)');
            ylabel('coherency');
            title(sprintf('%s %s psth - %s evoked potential coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = vertcat(fface, fnface);
            figData.y = [Cface, Cnface];
            drawnow;
            saveFigure(outDir,sprintf('coh_%s_%s_PSTH-%s_EVOKED_FACEvsNON%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);
          end
        end

        %todo: intra-channel spike-field "JPSTH"
        for unit_i = 1:length(channelUnitNames{channel_i})
          if length(channelUnitNames{channel_i}) == 2 && unit_i == 1 %skip unsorted if no isolated units; just do mua
            continue
          end
          if calcSwitch.spikeTimes
            [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgramcpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
              spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i}',movingWin,chr_params);
          else
            [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgramcpb(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
              spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',movingWin,chr_params);
          end

          t = tface - lfpAlignParams.msPreAlign;
          f = 1000*fface;

          fh = figure(); 
          imagesc(t,f,Cface'); axis xy
          xlabel('Time (ms)'); 
          ylabel('Frequency (Hz)');
          c = colorbar();
          ylabel(c,'Coherency');
          hold on
          draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
          draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
          title(sprintf('%s %s - %s field coherence, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
          clear figData
          figData.x = t;
          figData.y = f;
          figData.z = Cface';
          drawnow;
          saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_LFP_FACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
          close(fh);

          fh = figure(); 
          imagesc(t,f,phiface'); axis xy
          xlabel('Time (ms)'); 
          ylabel('Frequency (Hz)');
          c = colorbar();
          ylabel(c,'phase');
          hold on
          draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
          draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
          title(sprintf('%s %s - %s phase, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
          clear figData
          figData.x = t;
          figData.y = fface;
          figData.z = phiface';
          drawnow;
          saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_LFP_FACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
          close(fh);

          % nonface
          if calcSwitch.spikeTimes
            [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgramcpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
              spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i}', movingWin,chr_params);
          else
            [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgramcpb(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
              spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}', movingWin,chr_params);
          end

          t = tnface - lfpAlignParams.msPreAlign;
          f = 1000*fnface;

          fh = figure(); 
          imagesc(t,f,Cnface'); axis xy
          xlabel('Time (ms)'); 
          ylabel('Frequency (Hz)');
          c = colorbar();
          ylabel(c,'Coherency');
          hold on
          draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
          draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
          title(sprintf('%s %s - %s LFP coherence, nonface%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
          clear figData
          figData.x = t;
          figData.y = f;
          figData.z = Cnface';
          drawnow;
          saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_LFP_NONFACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
          close(fh);

          fh = figure(); 
          imagesc(t,f,phinface'); axis xy
          xlabel('Time (ms)'); 
          ylabel('Frequency (Hz)');
          c = colorbar();
          ylabel(c,'phase');
          hold on
          draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
          draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
          title(sprintf('%s %s - %s LFP phase, nonface%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
          clear figData
          figData.x = t;
          figData.y = fnface;
          figData.z = phinface';
          drawnow;
          saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_LFP_NONFACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
          close(fh);
        end

        %intra-channel tf, spike-spike (e.g. isolated vs. unsorted, u1 v u2, etc.)

        for unit_i = 1:length(channelUnitNames{channel_i})-1 %note: don't consider mua, b/c would have contribution of spike coupling to itself
          for unit2_i = unit_i+1:length(channelUnitNames{channel_i})-1
            %todo: JPSTH, cross-correllogram/covariogram
            if calcSwitch.spikeTimes
              [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypt(spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},...
                spikesByCategoryForTF{faceCatNum}{channel_i}{unit2_i},chr_params); %can have time grid as additional arg
              [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencypt(spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},...
                spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit2_i},chr_params); %can have time grid as additional arg
            else
              [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypb(spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',...
                spikesByCategoryBinned{faceCatNum}{channel_i}{unit2_i}',chr_params); 
              [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencypb(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}',...
                spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit2_i}',chr_params); 
            end
            errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
            fh = figure();
            mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
            legend('face', 'nonface');
            xlabel('frequency (Hz)');
            ylabel('coherency');
            title(sprintf('%s %s - %s %s coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = vertcat(fface, fnface); %todo: fix; add freq cutoff array slice [ applies to all coherence figures ]
            figData.y = [Cface, Cnface];
            drawnow;
            saveFigure(outDir,sprintf('coh_%s_%s-%s_%s_FACEvsNON%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);

            if calcSwitch.meanEvokedTF
              if calcSwitch.spikeTimes
                [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypt(allSpikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},...
                  allSpikesByCategoryForTF{faceCatNum}{channel_i}{unit2_i},chr_params); 
                [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencypt(allSpikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},...
                  allSpikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit2_i},chr_params);
              else
                [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypb(mean(spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i},3)',...
                  mean(spikesByCategoryBinned{faceCatNum}{channel_i}{unit2_i},3)',chr_params); 
                [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencypb(mean(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i},3)',...
                  mean(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit2_i},3)',chr_params);
              end
              errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
              fh = figure();
              mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
              legend('face', 'nonface');
              xlabel('frequency (Hz)');
              ylabel('coherency');
              title(sprintf('%s %s PSTH - %s %s PSTH coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = vertcat(fface, fnface); %todo: fix; add freq cutoff array slice [ applies to all coherence figures ]
              figData.y = [Cface, Cnface];
              drawnow;
              saveFigure(outDir,sprintf('coh_%s_%s_PSTH-%s_%s_PSTH_FACEvsNON%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum));
              close(fh);
            end

            % spike-spike time-frequency coherency, face vs. non
            % face
            if calcSwitch.spikeTimes
              [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgrampt(spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},...
                spikesByCategoryForTF{faceCatNum}{channel_i}{unit2_i}, movingWin,chr_params);
            else
              [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgrampb(spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',...
                spikesByCategoryBinned{faceCatNum}{channel_i}{unit2_i}', movingWin,chr_params);
            end

            t = tface - lfpAlignParams.msPreAlign;
            f = 1000*fface;

            fh = figure(); 
            imagesc(t,f,Cface'); axis xy
            xlabel('Time (ms)'); 
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'Coherency');
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s %s - %s %s coherence, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = t;
            figData.y = f;
            figData.z = Cface';
            drawnow;
            saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_%s_FACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);

            fh = figure(); 
            imagesc(t,f,phiface'); axis xy
            xlabel('Time (ms)'); 
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'phase');
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s %s - %s %s phase, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = t;
            figData.y = fface;
            figData.z = phiface';
            drawnow;
            saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_%s_FACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);

            % nonface
            if calcSwitch.spikeTimes
              [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgrampt(spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},...
                spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit2_i}, movingWin,chr_params);
            else
              [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgrampb(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}',...
                spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit2_i}', movingWin,chr_params);
            end

            t = tnface - lfpAlignParams.msPreAlign;
            f = 1000*fnface;

            fh = figure(); 
            imagesc(t,f,Cnface'); axis xy
            xlabel('Time (ms)'); 
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'Coherency');
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s %s - %s %s coherence, nonface%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i}),'FontSize',18);
            clear figData
            figData.x = t;
            figData.y = f;
            figData.z = Cnface';
            drawnow;
            saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_%s_NONFACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);

            fh = figure(); 
            imagesc(t,f,phinface'); axis xy
            xlabel('Time (ms)'); 
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'phase');
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s %s - %s %s phase, nonface%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = t;
            figData.y = fnface;
            figData.z = phinface';
            drawnow;
            saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_%s_NONFACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);
          end
        end
        %%%%%%
        %%%%%%
        % across channel
        for channel2_i = channel_i:length(lfpChannels) %these two lines make sure we only calculate once for each pair
          if channel2_i > channel_i 
            % channel2_i spike -- channel_i field
            for unit_i = 1:length(channelUnitNames{channel2_i})
              if length(channelUnitNames{channel2_i}) == 2 && unit_i == 1 %skip unsorted if no isolated units; just do mua
                continue
              end

              if calcSwitch.useJacknife
                chr_params.err = [2 .05];
                if calcSwitch.spikeTimes
                  [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                    spikesByCategoryForTF{faceCatNum}{channel2_i}{unit_i},chr_params);
                  [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                    spikesByCategoryForTF{nonfaceCatNum}{channel2_i}{unit_i},chr_params);
                else
                  [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpb(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                    spikesByCategoryBinned{faceCatNum}{channel2_i}{unit_i}',chr_params);
                  [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpb(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                    spikesByCategoryBinned{nonfaceCatNum}{channel2_i}{unit_i}',chr_params); 
                end
                errs = vertcat(faceErrs,nfaceErrs);
              else
                if calcSwitch.spikeTimes
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                    spikesByCategoryForTF{faceCatNum}{channel2_i}{unit_i},chr_params);
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                    spikesByCategoryForTF{nonfaceCatNum}{channel2_i}{unit_i},chr_params);
                else
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpb(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                    spikesByCategoryBinned{faceCatNum}{channel2_i}{unit_i}',chr_params);
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpb(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                    spikesByCategoryBinned{nonfaceCatNum}{channel2_i}{unit_i}',chr_params);
                end
                errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
              end
              phiface = -1*phiface;
              phinface = -1*phinface;

              fh = figure();
              mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
              legend('face', 'nonface');
              xlabel('frequency (Hz)');
              ylabel('coherency');
              title(sprintf('%s %s - %s field coherence%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = vertcat(fface, fnface);
              figData.y = [Cface, Cnface];
              drawnow;
              saveFigure(outDir,sprintf('coh_%s_%s-%sLFP_FACEvsNON%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);

              if calcSwitch.meanEvokedTF
                if calcSwitch.spikeTimes
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(mean(lfpByCategory{faceCatNum}(1,channel_i,:,:),3)),...
                    allSpikesByCategoryForTF{faceCatNum}{channel2_i}{unit_i},chr_params);
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:),3)),...
                    allSpikesByCategoryForTF{nonfaceCatNum}{channel2_i}{unit_i},chr_params);
                else
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpb(squeeze(mean(lfpByCategory{faceCatNum}(1,channel_i,:,:),3)),...
                    mean(spikesByCategoryBinned{faceCatNum}{channel2_i}{unit_i},1)',chr_params);
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpb(squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:),3)),...
                    mean(spikesByCategoryBinned{nonfaceCatNum}{channel2_i}{unit_i},1)',chr_params);
                end
                errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));

                phiface = -1*phiface;
                phinface = -1*phinface;

                fh = figure();
                mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
                legend('face', 'nonface');
                xlabel('frequency (Hz)');
                ylabel('coherency');
                title(sprintf('%s %s psth - %s evoked potential coherence%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = vertcat(fface, fnface);
                figData.y = [Cface, Cnface];
                drawnow;
                saveFigure(outDir,sprintf('coh_%s_%s_PSTH-%s_EVOKED_FACEvsNON%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
                close(fh);
              end

              %time-frequency coherency, channel2_i spike channel_i field tf
              if calcSwitch.spikeTimes
                [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgramcpt(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                  spikesByCategoryForTF{faceCatNum}{channel2_i}{unit_i},movingWin,chr_params);
              else
                [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgramcpb(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
                  spikesByCategoryBinned{faceCatNum}{channel2_i}{unit_i}',movingWin,chr_params);
              end

              t = tface - lfpAlignParams.msPreAlign;
              f = 1000*fface;

              fh = figure(); 
              imagesc(t,f,Cface'); axis xy
              xlabel('Time (ms)'); 
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'Coherency');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s %s - %s field coherence, face%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = Cface';
              drawnow;
              saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_LFP_FACE%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);

              fh = figure(); 
              imagesc(t,f,phiface'); axis xy
              xlabel('Time (ms)'); 
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'phase');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s %s - %s phase, face%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = fface;
              figData.z = phiface';
              drawnow;
              saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_LFP_FACE%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);

              % nonface
              if calcSwitch.spikeTimes
                [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgramcpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                  spikesByCategoryForTF{nonfaceCatNum}{channel2_i}{unit_i}, movingWin,chr_params);
              else
                [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgramcpb(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
                  spikesByCategoryBinned{nonfaceCatNum}{channel2_i}{unit_i}', movingWin,chr_params);
              end

              t = tnface - lfpAlignParams.msPreAlign;
              f = 1000*fnface;

              fh = figure(); 
              imagesc(t,f,Cnface'); axis xy
              xlabel('Time (ms)'); 
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'Coherency');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s %s - %s LFP coherence, nonface%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = Cnface';
              drawnow;
              saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_LFP_NONFACE%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);

              fh = figure(); 
              imagesc(t,f,phinface'); axis xy
              xlabel('Time (ms)'); 
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'phase');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s %s - %s LFP phase, nonface%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = fnface;
              figData.z = phinface';
              drawnow;
              saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_LFP_NONFACE%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);
            end


            % channel_i spike -- channel2_i field
            for unit_i = 1:length(channelUnitNames{channel_i})
              if length(channelUnitNames{channel_i}) == 2 && unit_i == 1 %skip unsorted if no isolated units; just do mua
                continue
              end

              if calcSwitch.useJacknife
                chr_params.err = [2 .05];
                if calcSwitch.spikeTimes
                  [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
                    spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},chr_params);
                  [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
                    spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},chr_params);
                else
                  [Cface,phi,S12,S1,S2,fface,zerosp,confCface,phistd, faceErrs]=coherencycpb(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
                    spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',chr_params);
                  [Cnface,phi,S12,S1,S2,fnface,zerosp,confCnface,phistd, nfaceErrs]= coherencycpb(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
                    spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}',chr_params); 
                end
                errs = vertcat(faceErrs,nfaceErrs);
              else
                if calcSwitch.spikeTimes
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
                    spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},chr_params);
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
                    spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},chr_params);
                else
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpb(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
                    spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',chr_params);
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpb(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
                    spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}',chr_params);
                end
                errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
              end
              phiface = -1*phiface;
              phinface = -1*phinface;

              fh = figure();
              mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
              legend('face', 'nonface');
              xlabel('frequency (Hz)');
              ylabel('coherency');
              title(sprintf('%s %s - %s field coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = vertcat(fface, fnface);
              figData.y = [Cface, Cnface];
              drawnow;
              saveFigure(outDir,sprintf('coh_%s_%s-%sLFP_FACEvsNON%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);

              if calcSwitch.meanEvokedTF
                if calcSwitch.spikeTimes
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpt(squeeze(mean(lfpByCategory{faceCatNum}(1,channel2_i,:,:),3)),...
                    allSpikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},chr_params);
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpt(squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:),3)),...
                    allSpikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},chr_params);
                else
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencycpb(squeeze(mean(lfpByCategory{faceCatNum}(1,channel2_i,:,:),3)),...
                    mean(spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i},1)',chr_params);
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencycpb(squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:),3)),...
                    mean(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i},1)',chr_params);
                end
                errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));

                phiface = -1*phiface;
                phinface = -1*phinface;

                fh = figure();
                mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
                legend('face', 'nonface');
                xlabel('frequency (Hz)');
                ylabel('coherency');
                title(sprintf('%s %s psth - %s evoked potential coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = vertcat(fface, fnface);
                figData.y = [Cface, Cnface];
                drawnow;
                saveFigure(outDir,sprintf('coh_%s_%s_PSTH-%s_EVOKED_FACEvsNON%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
                close(fh);
              end

              %time-frequency coherency, channel_i spike channel2_i field
              if calcSwitch.spikeTimes
                [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgramcpt(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
                  spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},movingWin,chr_params);
              else
                [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgramcpb(squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',...
                  spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',movingWin,chr_params);
              end

              t = tface - lfpAlignParams.msPreAlign;
              f = 1000*fface;

              fh = figure(); 
              imagesc(t,f,Cface'); axis xy
              xlabel('Time (ms)'); 
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'Coherency');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s %s - %s field coherence, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = Cface';
              drawnow;
              saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_LFP_FACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);

              fh = figure(); 
              imagesc(t,f,phiface'); axis xy
              xlabel('Time (ms)'); 
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'phase');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s %s - %s phase, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = fface;
              figData.z = phiface';
              drawnow;
              saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_LFP_FACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);

              % nonface
              if calcSwitch.spikeTimes
                [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgramcpt(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
                  spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i}, movingWin,chr_params);
              else
                [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgramcpb(squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',...
                  spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}', movingWin,chr_params);
              end

              t = tnface - lfpAlignParams.msPreAlign;
              f = 1000*fnface;

              fh = figure(); 
              imagesc(t,f,Cnface'); axis xy
              xlabel('Time (ms)'); 
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'Coherency');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s %s - %s LFP coherence, nonface%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = Cnface';
              drawnow;
              saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_LFP_NONFACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);

              fh = figure(); 
              imagesc(t,f,phinface'); axis xy
              xlabel('Time (ms)'); 
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'phase');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s %s - %s LFP phase, nonface%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = fnface;
              figData.z = phinface';
              drawnow;
              saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_LFP_NONFACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);
            end

            % field-field
            [Cface,phiface,S12,S1,S2,fface,confCface,phistdface]=coherencyc(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
              squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',chr_params);
            [Cnface,phinface,S12,S1,S2,fnface,confCnface,phistdnface]=coherencyc(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
              squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',chr_params);
            errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
            fh = figure();
            mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
            legend('face', 'nonface');
            xlabel('frequency (Hz)');
            ylabel('coherency');
            title(sprintf('%s field - %s field coherence%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = vertcat(fface, fnface);
            figData.y = [Cface, Cnface];
            drawnow;
            saveFigure(outDir,sprintf('coh_%s_LFP-%s_LFP_FACEvsNON%s_Run%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);

            if calcSwitch.meanEvokedTF
              [Cface,phiface,S12,S1,S2,fface,confCface,phistdface]=coherencyc(squeeze(mean(lfpByCategory{faceCatNum}(1,channel_i,:,:),3)),...
              squeeze(mean(lfpByCategory{faceCatNum}(1,channel2_i,:,:),3)),chr_params);
              [Cnface,phinface,S12,S1,S2,fnface,confCnface,phistdnface]=coherencyc(squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:),3)),...
                squeeze(mean(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:),3)),chr_params);
              errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
              fh = figure();
              mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
              legend('face', 'nonface');
              xlabel('frequency (Hz)');
              ylabel('coherency');
              title(sprintf('%s evoked potential - %s evoked potential coherence%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = vertcat(fface, fnface);
              figData.y = [Cface, Cnface];
              drawnow;
              saveFigure(outDir,sprintf('coh_%s_EVOKED-%s_EVOKED_FACEvsNON%s_Run%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
              close(fh);
            end

            % field-field time-frequency coherency, face
            [Cface,phiface,S12,S1,S2,tface,fface,confCface,phistdface]=cohgramc(squeeze(lfpByCategory{faceCatNum}(1,channel_i,:,:))',...
              squeeze(lfpByCategory{faceCatNum}(1,channel2_i,:,:))',movingWin,chr_params);
            t = tface - lfpAlignParams.msPreAlign;
            f = 1000*fface;

            fh = figure(); 
            imagesc(t,f,Cface'); axis xy
            xlabel('Time (ms)'); 
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'Coherency');
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s field - %s field coherence, face%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = t;
            figData.y = f;
            figData.z = Cface';
            drawnow;
            saveFigure(outDir,sprintf('coh_TF_%s_LFP-%s_LFP_FACE%s_Run%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);

            fh = figure(); 
            imagesc(t,f,phiface'); axis xy
            xlabel('Time (ms)'); 
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'phase');
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s field - %s field phase, face%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = t;
            figData.y = fface;
            figData.z = phiface';
            drawnow;
            saveFigure(outDir,sprintf('phase_TF_%s_LFP-%s_LFP_FACE%s_Run%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);

            % lfp-lfp time-frequency coherency, nonface
            [Cnface,phinface,S12,S1,S2,tnface,fnface,confCnface,phistdnface]=cohgramc(squeeze(lfpByCategory{nonfaceCatNum}(1,channel_i,:,:))',...
              squeeze(lfpByCategory{nonfaceCatNum}(1,channel2_i,:,:))',movingWin,chr_params);
            t = tnface - lfpAlignParams.msPreAlign;
            f = 1000*fnface;
            fh = figure();
            imagesc(t,f,Cnface'); axis xy
            xlabel('Time (ms)');
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'Coherency'); 
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s field - %s field coherence, nonface%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = t;
            figData.y = f;
            figData.z = Cnface';
            drawnow;
            saveFigure(outDir,sprintf('coh_TF_%s_LFP-%s_LFP_NONFACE%s_Run%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);

            fh = figure(); 
            imagesc(t,f,phinface'); axis xy
            xlabel('Time (ms)'); 
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'phase');
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s field - %s field phase, nonface%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = t;
            figData.y = f;
            figData.z = phinface';
            drawnow;
            saveFigure(outDir,sprintf('phase_TF_%s_LFP-%s_LFP_NONFACE%s_Run%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
            close(fh);

            % spike-spike coherency, face vs. non

            for unit_i = 1:length(channelUnitNames{channel_i})
              if length(channelUnitNames{channel_i}) == 1 && unit_i == 1
                continue
              end
              for unit2_i = 1:length(channelUnitNames{channel2_i})
                if length(channelUnitNames{channel2_i}) == 1 && unit2_i == 1
                  continue
                end
                if calcSwitch.spikeTimes
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypt(spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},...
                    spikesByCategoryForTF{faceCatNum}{channel2_i}{unit2_i},chr_params); %can have time grid as additional arg
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencypt(spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},...
                    spikesByCategoryForTF{nonfaceCatNum}{channel2_i}{unit2_i},chr_params); %can have time grid as additional arg
                else
                  [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypb(spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',...
                    spikesByCategoryBinned{faceCatNum}{channel2_i}{unit2_i}',chr_params); 
                  [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencypb(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}',...
                    spikesByCategoryBinned{nonfaceCatNum}{channel2_i}{unit2_i}',chr_params); 
                end
                errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
                fh = figure();
                mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
                legend('face', 'nonface');
                xlabel('frequency (Hz)');
                ylabel('coherency');
                title(sprintf('%s %s - %s %s coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = vertcat(fface, fnface); %todo: fix; add freq cutoff array slice [ applies to all coherence figures ]
                figData.y = [Cface, Cnface];
                drawnow;
                saveFigure(outDir,sprintf('coh_%s_%s-%s_%s_FACEvsNON%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
                close(fh);

                if calcSwitch.meanEvokedTF
                  if calcSwitch.spikeTimes
                    [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypt(allSpikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},...
                      allSpikesByCategoryForTF{faceCatNum}{channel2_i}{unit2_i},chr_params); 
                    [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencypt(allSpikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},...
                      allSpikesByCategoryForTF{nonfaceCatNum}{channel2_i}{unit2_i},chr_params);
                  else
                    [Cface,phiface,S12,S1,S2,fface,zerosp,confCface,phistdface]=coherencypb(mean(spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i},3)',...
                      mean(spikesByCategoryBinned{faceCatNum}{channel2_i}{unit2_i},3)',chr_params); 
                    [Cnface,phinface,S12,S1,S2,fnface,zerosp,confCnface,phistdnface]=coherencypb(mean(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i},3)',...
                      mean(spikesByCategoryBinned{nonfaceCatNum}{channel2_i}{unit2_i},3)',chr_params);
                  end
                  errs = vertcat(confCface*ones(1,length(Cface(fface < 0.1))),confCnface*ones(1,length(Cface(fface < 0.1))));
                  fh = figure();
                  mseb(1000*repmat(fface(fface < 0.1),2,1),vertcat((Cface(fface < 0.1))',(Cnface(fface < 0.1))'), errs);
                  legend('face', 'nonface');
                  xlabel('frequency (Hz)');
                  ylabel('coherency');
                  title(sprintf('%s %s PSTH - %s %s PSTH coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                  clear figData
                  figData.x = vertcat(fface, fnface); %todo: fix; add freq cutoff array slice [ applies to all coherence figures ]
                  figData.y = [Cface, Cnface];
                  drawnow;
                  saveFigure(outDir,sprintf('coh_%s_%s_PSTH-%s_%s_PSTH_FACEvsNON%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum));
                  close(fh);
                end

                % spike-spike time-frequency coherency, face vs. non
                % face
                if calcSwitch.spikeTimes
                  [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgrampt(spikesByCategoryForTF{faceCatNum}{channel_i}{unit_i},...
                    spikesByCategoryForTF{faceCatNum}{channel2_i}{unit2_i}, movingWin,chr_params);
                else
                  [Cface,phiface,S12,S1,S2,tface,fface,zerosp,confCface,phistdface]=cohgrampb(spikesByCategoryBinned{faceCatNum}{channel_i}{unit_i}',...
                    spikesByCategoryBinned{faceCatNum}{channel2_i}{unit2_i}', movingWin,chr_params);
                end

                t = tface - lfpAlignParams.msPreAlign;
                f = 1000*fface;

                fh = figure(); 
                imagesc(t,f,Cface'); axis xy
                xlabel('Time (ms)'); 
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'Coherency');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s %s - %s %s coherence, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = f;
                figData.z = Cface';
                drawnow;
                saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_%s_FACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
                close(fh);

                fh = figure(); 
                imagesc(t,f,phiface'); axis xy
                xlabel('Time (ms)'); 
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'phase');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s %s - %s %s phase, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = fface;
                figData.z = phiface';
                drawnow;
                saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_%s_FACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
                close(fh);

                % nonface
                if calcSwitch.spikeTimes
                  [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgrampt(spikesByCategoryForTF{nonfaceCatNum}{channel_i}{unit_i},...
                    spikesByCategoryForTF{nonfaceCatNum}{channel2_i}{unit2_i}, movingWin,chr_params);
                else
                  [Cnface,phinface,S12,S1,S2,tnface,fnface,zerosp,confCface,phistdface]=cohgrampb(spikesByCategoryBinned{nonfaceCatNum}{channel_i}{unit_i}',...
                    spikesByCategoryBinned{nonfaceCatNum}{channel2_i}{unit2_i}', movingWin,chr_params);
                end

                t = tnface - lfpAlignParams.msPreAlign;
                f = 1000*fnface;

                fh = figure(); 
                imagesc(t,f,Cnface'); axis xy
                xlabel('Time (ms)'); 
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'Coherency');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s %s - %s %s coherence, nonface%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = f;
                figData.z = Cnface';
                drawnow;
                saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_%s_NONFACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
                close(fh);

                fh = figure(); 
                imagesc(t,f,phinface'); axis xy
                xlabel('Time (ms)'); 
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'phase');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s %s - %s %s phase, nonface%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = fnface;
                figData.z = phinface';
                drawnow;
                saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_%s_NONFACE%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum) );
                close(fh);
              end
            end
          end
        end
      end
      % clean up temporary variables and restore stable variables
      if strcmp(tfCalcSwitchNames{calc_i},'inducedCatTF')
        spikesByCategoryBinned = spikesByCategoryBinnedEvokedTmp;
        lfpByCategory =  lfpByCategoryEvokedTmp;
      end
    end
  end
end
end

