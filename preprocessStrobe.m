function [ strobeTriggers ] = preprocessStrobe( inputDataSource, params )
%preprocessPhotodiodeStrobe extracts frame times from a photodiode signal 
%   Inputs:
%   - inputDataSource: may be:
%       - a filename, if inputDataType is either blackrockFilename or
%         that takes filenames
%       - a 1xN array of data samples
%   - params: struct with the fields
%     - needStrobe: return immediately if 0 (type: logical)
%     - inputDataType: options are 'array','blackrockFilename', or
%                      'otherFilename' (type: string)
%     - dataChannel: this is the channel number of the data to preprocess, if inputDataType is blackrockFilename or a custom handle (type: int)
%     - dataLoader: handle to function with signature dataTrace,samplingFreq] = f(inputDataSource, params); params passed
%                   unabridged (required only if inputDataSource is otherFilename)
%     - peakTimeOffset: this is the offset, in ms, of the peak from the event it's coupled to
%                       (note: will be subtracted, so should be > 0 if
%                           peak follows event)
%                           (type: numeric)
%     - strobeTroughs: set to 1 if strobe causes troughs, 0 else (type: logical)
%     - levelCalibType: (type: string) options are:
%       - 'hardcode'
%       - 'hardcodeAndPlot'
%       - 'hardcodeAndCheck'
%       - 'auto'
%       - 'autoAndPlot'
%       - 'autoAndCheck'
%       - 'manual'
%     - cleanPeaks: if true, keep peaks only if the signal dropped below the low-peak threshold between them;
%                   for a series of peaks without a sub-threshold trough, keep only the highest (type: logical)
%     - useRisingEdge: define threshold using peaks, then use rising edge
%       threshold crossing for trigger times; optional, default 0 (type: logical)
%     for trigger times
%     - checkHighLowAlternation: (optional, default zero, only used if numLevels == 2)
%     - numLevels: (type: int)
%     - minPeakNumInLevel: min number of peaks required to define a level (type: int)
%     - samplingFreq: not needed if inputDataType is blackrockFilename, or
%                     a custom handle that sets sampling freq (numeric)
%     - peaksToPlot: number of peaks to show in calibration plots; suggested value > 10 (int)
%     - peakFreq: approximate number of peaks per second
%     - noisePeaksAtPeak: number of peaks above threshold per strobe flash; default 0 (type: int)
%     - hardcodeFromFile: (type: logical) (only needed for hardcode calib types)
%     - displayStats: print info on strobe rate and variability. Optional, default 1. (type: logical)
%     - saveCalibFile: (type: logical)
%     - inputCalibrationFile: filename ending in .mat, contains field lowThreshold (numeric)
%                             if numLevels >= 2, also contains field highThreshold (numeric)
%                             if numLevels == 3, also contains field midThreshold (numeric)
%                             (required only if hardcodeFromFile): (type: string)    
%     - outputCalibrationFile: filename ending in .mat (required only if saveCalibFile): (type:string)
%     - saveFigures: (type: logical)
%     - calibFigFname: figure filename, without .fig stem; required only if saveFigures (type:string)
%     - triggersFigFname: figure filename, without .fig stem; required only if saveFigures (type:string)
%     - dateSubject: MMDDYYNAME or similar; required only if saveFigures (type: string)
%     - runNum: eg 002; required only if saveFigures (type: string)
%     - outDir: directory for figure and calib. output; include the final '/' in the path (type: string) 
%               (required only if saveFigures or saveCalibFile is 1)
%
%   Outputs:
%     - strobeTriggers: a struct of peak times in ms, with fields:
%         - all
%         - if numLevels == 2, all, low, high, lowToHigh, highToLow
%         - if numLevels == 3, all, low, mid, high, lowToHigh, highToLow,
%           lowToMid, midToLow, midToHigh, highToMid
%   Possible extensions:
%     - make calibration plot smarter, so that it plots the region of the signal with the greatest possible number of digital levels 
%     - apply parameter defaults, check parameter validity (could also be handled in checkAnalysisParamFile

if ~params.needStrobe
  strobeTriggers = struct();
  return
end

% unpack parameters
levelCalibType = params.levelCalibType;
peaksToPlot = params.peaksToPlot;
cleanPeaks = params.cleanPeaks;
numLevels = params.numLevels;
inputDataType = params.inputDataType;
peakFreq = params.peakFreq;
minPeakNumInLevel = params.minPeakNumInLevel;
saveFigures = params.saveFigures;
saveCalibFile = params.saveCalibFile;
strobeTroughs = params.strobeTroughs;
if isfield(params, 'peakTimeOffset')
  peakTimeOffset = params.peakTimeOffset;
else
  peakTimeOffset = 0;
end
if params.numLevels == 2
  checkHighLowAlternation = params.checkHighLowAlternation;
end
if params.saveFigures || params.saveCalibFile
  outDir = params.outDir;
  runNum = params.runNum;
  dateSubject = params.dateSubject;
end
if saveFigures
  calibFigFname = params.calibFigFname;
  triggersFigFname = params.triggersFigFname;
end
if ~strcmp(inputDataType, 'array')
  dataChannel = params.dataChannel;
else
  samplingFreq = params.samplingFreq;
end
if strcmp(inputDataType, 'otherFilename')
  dataLoader = params.dataLoader;
end
if params.saveCalibFile
  outputCalibrationFile = params.outputCalibrationFile;
end
if isfield(params,'noisePeaksAtPeak')
  noisePeaksAtPeak = params.noisePeaksAtPeak;
else
  noisePeaksAtPeak = 0;
end
if isfield(params,'useRisingEdge')
  useRisingEdge = params.useRisingEdge;
else
  useRisingEdge = 0;
end
if isfield(params,'displayStats')
  displayStats = params.displayStats;
else
  displayStats = 1;
end

% load photodiode data, if needed
if strcmp(inputDataType, 'blackrockFilename')
  assert(logical(exist(inputDataSource,'file')),'The data file you requested does not exist.');
  tmp = openNSx(inputDataSource,'report','read');
  samplingFreq = tmp.MetaTags.SamplingFreq;
  dataChannelIndex = find(tmp.MetaTags.ChannelID == dataChannel);
  dataTrace = tmp.Data(dataChannelIndex,:); %#ok
  assert(size(dataTrace,1) == 1, 'Invalid data channel; requested channel not recorded in file provided');
  clear tmp
elseif strcmp(inputDataType, 'otherFilename')
  [dataTrace, samplingFreq] = dataLoader(inputDataSource, params);
else
  dataTrace = inputDataSource;
end

frameNumCeiling = 2*(noisePeaksAtPeak+1)*ceil(length(dataTrace)*peakFreq/samplingFreq); %keep some non-frame peaks for numLevels=1 auto calib 

        
% Preprocessing logic:
% 1) Find peaks in the inverted diode trace. A peak occurs when the cathode ray passes the diode,
%    even on black pixels. However, additional small-prominence peaks occur due to noise
% 2) Sort peaks by their height: since both the the signal and its
%    second derivative are high at true peaks, they should be the largest
if strobeTroughs
  dataTrace = -1*(dataTrace - mean(dataTrace)); %this makes frames peaks
else
  dataTrace = dataTrace - mean(dataTrace);
end
rawPeaks = findpeaks(dataTrace);
rawPeaks.loc = rawPeaks.loc(rawPeaks.loc > 1); % don't count the first sample as a peak
rawPeaks.loc = rawPeaks.loc(rawPeaks.loc < length(dataTrace)); % don't count the last sample as a peak
rawPeaksSorted = sort(dataTrace(rawPeaks.loc),'descend');
rawPeaksSorted = rawPeaksSorted(1:frameNumCeiling);

if any(strcmp(levelCalibType, {'hardcodeAndPlot','hardcodeAndCheck','hardcode'}))  
  if params.hardcodeFromFile
    switch numLevels
      case 1
        load(params.inputCalibrationFile,'lowThreshold');
      case 2
        load(params.inputCalibrationFile,'highThreshold','lowThreshold');
      case 3
        load(params.inputCalibrationFile,'highThreshold','midThreshold','lowThreshold');
    end
  else
    switch numLevels
      case 1
        highThreshold = params.levelHigh;
      case 2
        highThreshold = params.levelHigh;
        lowThreshold = params.levelLow;
      case 3
        highThreshold = params.levelHigh;
        midThreshold = params.levelMid;
        lowThreshold = params.levelLow;
    end
  end
end

if any(strcmp(levelCalibType, {'auto','autoAndPlot','autoAndCheck'}))
  % Auto-calibration logic: 
  %  The largest jumps in the plot of sorted peak heights should occur between true peak types,  
  %  or between the smallest peak type and the largest noise peak. So find and sort those jumps, 
  %  then place level cut-offs at the mid-point of the largest ones. 

  peakDiffs = diff(rawPeaksSorted);
  [~,peakDiffSortInds] = sort(peakDiffs);
  switch numLevels
    case 1
      lowThresholdDiffInd = peakDiffSortInds(1);
      lowThreshold = rawPeaksSorted(lowThresholdDiffInd)*.5 + rawPeaksSorted(lowThresholdDiffInd+1)/2;
    case 2
      % defense against noise peak at intermediate point between levels
      tmpThreshold2Ind = 2;
      while abs(peakDiffSortInds(tmpThreshold2Ind) - peakDiffSortInds(1)) < minPeakNumInLevel
        tmpThreshold2Ind = tmpThreshold2Ind + 1;
      end
      % end defense
      tmpThreshold1 = rawPeaksSorted(peakDiffSortInds(1))*.5 + rawPeaksSorted(peakDiffSortInds(1)+1)/2;
      tmpThreshold2 = rawPeaksSorted(peakDiffSortInds(tmpThreshold2Ind))*.5 + rawPeaksSorted(peakDiffSortInds(tmpThreshold2Ind)+1)/2; 
      highThreshold = max(tmpThreshold1,tmpThreshold2);
      lowThreshold = min(tmpThreshold1,tmpThreshold2);
    case 3
      % defense against noise peak at intermediate point between levels
      tmpThreshold2Ind = 2;
      while abs(peakDiffSortInds(tmpThreshold2Ind) - peakDiffSortInds(1)) < minPeakNumInLevel
        tmpThreshold2Ind = tmpThreshold2Ind + 1;
      end
      tmpThreshold3Ind = tmpThreshold2Ind+1;
      while min(abs(peakDiffSortInds(tmpThreshold3Ind) - peakDiffSortInds(1)), abs(peakDiffSortInds(tmpThreshold3Ind) - peakDiffSortInds(tmpThreshold2Ind))) < minPeakNumInLevel
        tmpThreshold3Ind = tmpThreshold3Ind + 1;
      end
      %end defense
      tmpThreshold1 = rawPeaksSorted(peakDiffSortInds(1))*.5 + rawPeaksSorted(peakDiffSortInds(1)+1)/2;
      tmpThreshold2 = rawPeaksSorted(peakDiffSortInds(tmpThreshold2Ind))*.5 + rawPeaksSorted(peakDiffSortInds(tmpThreshold2Ind)+1)/2; 
      tmpThreshold3 = rawPeaksSorted(peakDiffSortInds(tmpThreshold3Ind))*.5 + rawPeaksSorted(peakDiffSortInds(tmpThreshold3Ind)+1)/2; 
      highThreshold = max([tmpThreshold1,tmpThreshold2,tmpThreshold3]);
      lowThreshold = min([tmpThreshold1,tmpThreshold2,tmpThreshold3]);
      midThreshold = setdiff([tmpThreshold1,tmpThreshold2,tmpThreshold3],[highThreshold,lowThreshold]);
  end
end

if ~strcmp(levelCalibType,'manual')  %todo: fix plottng, make alternation check optional
  peakSampleInds = rawPeaks.loc(dataTrace(rawPeaks.loc) > lowThreshold);
  if cleanPeaks
    peakSampleIndsClean = zeros(size(peakSampleInds));
    cleanPeaks_i = 1;
    lookAhead = 1;
    for peak_i = 1:length(peakSampleInds)
      if lookAhead > 1
        lookAhead = lookAhead - 1;
        continue
      end
      peakToKeep = peak_i;
      peakToKeepLevel = dataTrace(peakSampleInds(peak_i));
      while peak_i+lookAhead < length(peakSampleInds) && min(dataTrace(peakSampleInds(peak_i):peakSampleInds(peak_i+lookAhead))) > lowThreshold
        newPeakToKeepLevel = max(peakToKeepLevel,dataTrace(peakSampleInds(peak_i+lookAhead)));
        peakToKeep = peakToKeep*(newPeakToKeepLevel == peakToKeepLevel) + (peak_i+lookAhead)*(newPeakToKeepLevel ~= peakToKeepLevel);
        lookAhead = lookAhead + 1;
      end
      peakSampleIndsClean(cleanPeaks_i) = peakSampleInds(peakToKeep);
      cleanPeaks_i = cleanPeaks_i + 1;
    end
    peakSampleInds = peakSampleIndsClean(peakSampleIndsClean > 0);
  end

  if numLevels == 2
    highPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > highThreshold);
    lowPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > lowThreshold & dataTrace(peakSampleInds) < highThreshold);
  end

  if numLevels == 3
    highPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > highThreshold);
    midPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > midThreshold & dataTrace(peakSampleInds) < highThreshold); 
    lowPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > lowThreshold & dataTrace(peakSampleInds) < midThreshold); 
  end

  % if two levels, check whether frame levels alternate
  if numLevels == 2
    frameCountDiff = length(highPeakSampleInds) - length(lowPeakSampleInds);
    % with correct alternation, smallest interval between subsequent high peaks must be greater than longest frame; same for low peaks 
    frameAlternationError = (min(min(diff(highPeakSampleInds)),min(diff(lowPeakSampleInds))) > max(diff(peakSampleInds)));
  end

  % make plots, if requested, or if auto calib failed
  if any(strcmp(levelCalibType,{'autoAndPlot','autoAndCheck','hardcodeAndPlot','hardcodeAndCheck'})) || (numLevels == 2 && (abs(frameCountDiff) > 1 || frameAlternationError))
    if numLevels == 2 && abs(frameCountDiff) > 1 && checkHighLowAlternation
      fprintf('Frame count mismatch: Found %d bright-square frames and %d dark-square frames. Starting manual check.\n',...
        length(highPeakSampleInds),length(lowPeakSampleInds));
    end
    if numLevels == 2 && any(frameAlternationError) && checkHighLowAlternation
      fprintf('Found light-dark frame alternation error. Starting manual check.\n');
    end
    if any(strcmp(levelCalibType, {'hardcodeAndPlot','hardcodeAndCheck','hardcode'}))
      calibSuffix = 'hard code';
    else
      calibSuffix = 'auto';
    end

    figure('Name','Strobe Levels','NumberTitle','off');
    subplot(1,2,1);
    sampleRange = ceil(peakSampleInds(1) - 5*samplingFreq/peakFreq):ceil(peakSampleInds(1) + (peaksToPlot+0.5)*samplingFreq/peakFreq);
    sampleRange = sampleRange(sampleRange > 0);
    plot(sampleRange,dataTrace(sampleRange));
    hold on
    % plot threshold lines
    h = get(gca,'xlim');
    plot([h(1) h(2)], [lowThreshold lowThreshold]);
    forLegend121 = {'raw trace',sprintf('low threshold %s',calibSuffix)}; %note: 121 refers to the subplot
    if numLevels == 3
      plot([h(1) h(2)], [midThreshold midThreshold]);
      forLegend121 = horzcat(forLegend121, {sprintf('mid threshold %s',calibSuffix)});
    end
    if numLevels > 1
      plot([h(1) h(2)], [highThreshold highThreshold]);
      forLegend121 = horzcat(forLegend121, {sprintf('high threshold %s',calibSuffix)});
    end

    % plot peaks, labeled by high-mid-low if appropriate
    switch numLevels
      case 1
        plot(peakSampleInds(1:peaksToPlot+1), dataTrace(peakSampleInds(1:peaksToPlot+1)),'linestyle','none','marker','o','color','r');
        forLegend121 = horzcat(forLegend121, {'peaks'});
      case 2
        % high peaks
        highPeaksToPlot = highPeakSampleInds(highPeakSampleInds < ceil(peakSampleInds(1) + (peaksToPlot+0.5)*samplingFreq/peakFreq));
        plot(highPeaksToPlot, dataTrace(highPeaksToPlot),'linestyle','none','marker','o','color','r');
        forLegend121 = horzcat(forLegend121, {'high peaks'});
        % low peaks
        lowPeaksToPlot = lowPeakSampleInds(lowPeakSampleInds < ceil(peakSampleInds(1) + (peaksToPlot+0.5)*samplingFreq/peakFreq));
        plot(lowPeaksToPlot, dataTrace(lowPeaksToPlot),'linestyle','none','marker','o','color','g');
        forLegend121 = horzcat(forLegend121, {'low peaks'}); 
      case 3
        % high peaks
        highPeaksToPlot = highPeakSampleInds(highPeakSampleInds < ceil(peakSampleInds(1) + (peaksToPlot+0.5)*samplingFreq/peakFreq));
        plot(highPeaksToPlot, dataTrace(highPeaksToPlot),'linestyle','none','marker','o','color','r');
        forLegend121 = horzcat(forLegend121, {'high peaks'});
        % mid peaks
        midPeaksToPlot = midPeakSampleInds(midPeakSampleInds < ceil(peakSampleInds(1) + (peaksToPlot+0.5)*samplingFreq/peakFreq));
        plot(midPeaksToPlot, dataTrace(midPeaksToPlot),'linestyle','none','marker','o','color','k');
        forLegend121 = horzcat(forLegend121, {'mid peaks'});
        % low peaks
        lowPeaksToPlot = lowPeakSampleInds(lowPeakSampleInds < ceil(peakSampleInds(1) + (peaksToPlot+0.5)*samplingFreq/peakFreq));
        plot(lowPeaksToPlot, dataTrace(lowPeaksToPlot),'linestyle','none','marker','o','color','g');
        forLegend121 = horzcat(forLegend121, {'low peaks'});

    end
    legend(forLegend121,'Location','NorthOutside');  

    % show thresholds on ranked peaks plot
    subplot(1,2,2);
    plot(rawPeaksSorted, 'marker','o');
    hold on
    h = get(gca,'xlim');
    plot([h(1) h(2)], [lowThreshold lowThreshold]);
    forLegend122 = {'peak levels',sprintf('low threshold %s',calibSuffix)}; %note: 122 refers to the subplot
    if numLevels == 3
      plot([h(1) h(2)], [midThreshold midThreshold]);
      forLegend122 = horzcat(forLegend122, {sprintf('mid threshold %s',calibSuffix)});
    end
    if numLevels > 1
      plot([h(1) h(2)], [highThreshold highThreshold]);
      forLegend122 = horzcat(forLegend122, {sprintf('high threshold %s',calibSuffix)});
    end
    legend(forLegend122,'Location','NorthOutside');
    drawnow();

    if any(strcmp(levelCalibType,{'autoAndCheck','hardcodeAndCheck'}))
      tmp = input('Accept low threshold? Enter y or n: ','s');
      if strcmp(tmp,'n')
        tmp = input('Enter new low threshold, or ''quit'' to continue without photodiode sync: ','s');
        while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
          tmp = input('Invalid input: enter new low threshold, or ''quit'' to continue without photodiode sync: ','s');
        end
        if strcmp(tmp,'quit')
          return %todo: need to handle this better
        else
          lowThreshold = str2double(tmp);
          subplot(1,2,1);
          h = get(gca,'xlim');
          plot([h(1) h(2)], [lowThreshold lowThreshold]);
          %forLegend121 = horzcat(forLegend121,{'low threshold manual'});
          %legend(forLegend121);
          subplot(1,2,2);
          h = get(gca,'xlim');
          plot([h(1) h(2)], [lowThreshold lowThreshold]);
          %forLegend122 = horzcat(forLegend121,{'low threshold manual'});
          %legend(forLegend122);
        end
      end
      if numLevels == 3
        tmp = input('Accept mid threshold? Enter y or n: ','s');
        if strcmp(tmp,'n')
          tmp = input('Enter new mid threshold, or ''quit'' to continue without photodiode sync: ','s');
          while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
            tmp = input('Invalid input: enter new mid threshold, or ''quit'' to continue without photodiode sync: ','s');
          end
          if strcmp(tmp,'quit')
            return %todo: need to handle this better
          else
            midThreshold = str2double(tmp);
            subplot(1,2,1);
            h = get(gca,'xlim');
            plot([h(1) h(2)], [midThreshold midThreshold]);
            %forLegend121 = horzcat(forLegend121,{'mid threshold manual'});
            %legend(forLegend121);
            subplot(1,2,2);
            h = get(gca,'xlim');
            plot([h(1) h(2)], [midThreshold midThreshold]);
            %forLegend122 = horzcat(forLegend121,{'mid threshold manual'});
            %legend(forLegend122);
          end
        end
      end
      if numLevels > 1
        tmp = input('Accept high threshold? Enter y or n: ','s');
        if strcmp(tmp,'n')
          tmp = input('Enter new high threshold, or ''quit'' to continue without photodiode sync: ','s');
          while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
            tmp = input('Invalid input: enter new high threshold, or ''quit'' to continue without photodiode sync: ','s');
          end
          if strcmp(tmp,'quit')
            return %todo: need to handle this better
          else
            highThreshold = str2double(tmp);
            subplot(1,2,1);
            h = get(gca,'xlim');
            plot([h(1) h(2)], [highThreshold highThreshold]);
            forLegend121 = horzcat(forLegend121,{'high threshold manual'});
            legend(forLegend121);
            subplot(1,2,2);
            h = get(gca,'xlim');
            plot([h(1) h(2)], [highThreshold highThreshold]);
            forLegend122 = horzcat(forLegend121,{'high threshold manual'});
            legend(forLegend122);
          end
        end
      end
    end
  end

  % if two levels, check that frames alternate appropriately, if requested
  if numLevels == 2 && checkHighLowAlternation
    peakSampleInds = rawPeaks.loc(dataTrace(rawPeaks.loc) > lowThreshold);
    highPeakSampleInds = rawPeaks.loc(dataTrace(rawPeaks.loc) > highThreshold);
    lowPeakSampleInds = rawPeaks.loc(dataTrace(rawPeaks.loc) > lowThreshold & dataTrace(rawPeaks.loc) < highThreshold);
    frameCountDiff = length(highPeakSampleInds) - length(lowPeakSampleInds);
    % with correct alternation, smallest interval between subsequent high peaks must tbe greater than longest frame; same for low peaks 
    frameAlternationError = (min(min(diff(highPeakSampleInds)),min(diff(lowPeakSampleInds))) > max(diff(peakSampleInds)));
    notValid = (abs(frameCountDiff) > 1 || frameAlternationError);
    if notValid
      if abs(frameCountDiff) > 1
        tmp = input(sprintf('Frame count mismatch: Found %d bright-square frames and %d dark-square frames. Press ENTER to start manual check, or q to quit: ',...
          length(highPeakSampleInds),length(lowPeakSampleInds)));
      end
      if frameAlternationError
        tmp = input(sprintf('Found light-dark frame alternation error. Press ENTER to start manual check, or q to quit: '));
      end
      while ~strcmp('tmp',{'','q'})
        tmp = input('Invalid input. Press ENTER to re-start manual check, or q to quit: ');
      end
      if isempty(tmp)
        notValid = 1;
      else
        notValid = 0;
      end
    end
  else
    notValid = 0;
  end
end

if strcmp(levelCalibType,'manual') || (notValid && checkHighLowAlternation)
  notValid = 1;
  while notValid
    figure('Name','Strobe Calibration','NumberTitle','off');
    subplot(1,2,1);
    sampleRange = ceil(length(dataTrace)/2):ceil(length(dataTrace)/2 + peaksToPlot*samplingFreq/peakFreq);
    plot(sampleRange,dataTrace(sampleRange));
    hold on
    forLegend121 = {'raw trace'};
    legend(forLegend121,'Location','NorthOutside');
    subplot(1,2,2);
    plot(rawPeaksSorted, 'marker','o');
    hold on
    forLegend122 = {'peak levels'};
    legend(forLegend122,'Location','NorthOutside');
    drawnow();
    tmp = input('Enter the low threshold, or ''quit'' to continue without photodiode sync: ','s');
    while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
      tmp = input('Invalid input: enter the low threshold, or ''quit'' to continue without photodiode sync: ','s');
    end
    if strcmp(tmp,'quit')
      return %todo: need to handle this better
    else
      lowThreshold = str2double(tmp);
      subplot(1,2,1);
      h = get(gca,'xlim');
      plot([h(1) h(2)], [lowThreshold lowThreshold]);
      forLegend121 = horzcat(forLegend121,{'low threshold manual'});
      legend(forLegend121);
      subplot(1,2,2);
      h = get(gca,'xlim');
      plot([h(1) h(2)], [lowThreshold lowThreshold]);
      forLegend122 = horzcat(forLegend121,{'low threshold manual'});
      legend(forLegend122,'Location','NorthOutside');
    end
    if numLevels == 3
      tmp = input('Enter the mid threshold, or ''quit'' to continue without photodiode sync: ','s');
      while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
        tmp = input('Invalid input: enter the mid threshold, or ''quit'' to continue without photodiode sync: ','s');
      end
      if strcmp(tmp,'quit')
        return %todo: need to handle this better
      else
        midThreshold = str2double(tmp);
        subplot(1,2,1);
        h = get(gca,'xlim');
        plot([h(1) h(2)], [midThreshold midThreshold]);
        forLegend121 = horzcat(forLegend121,{'mid threshold manual'});
        legend(forLegend121);
        subplot(1,2,2);
        h = get(gca,'xlim');
        plot([h(1) h(2)], [midThreshold midThreshold]);
        forLegend122 = horzcat(forLegend121,{'mid threshold manual'});
        legend(forLegend122);
      end
    end
    if numLevels > 1
      tmp = input('Enter the high threshold, or ''quit'' to continue without photodiode sync: ','s');
      while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
        tmp = input('Invalid input: enter the high threshold, or ''quit'' to continue without photodiode sync: ','s');
      end
      if strcmp(tmp,'quit')
        return %todo: need to handle this better
      else
        highThreshold = str2double(tmp);
        subplot(1,2,1);
        h = get(gca,'xlim');
        plot([h(1) h(2)], [highThreshold highThreshold]);
        forLegend121 = horzcat(forLegend121,{'high threshold manual'});
        legend(forLegend121);
        subplot(1,2,2);
        h = get(gca,'xlim');
        plot([h(1) h(2)], [highThreshold highThreshold]);
        forLegend122 = horzcat(forLegend121,{'high threshold manual'});
        legend(forLegend122);
      end
    end
    % if two levels, check that frames alternate appropriately
    if numLevels == 2 && checkHighLowAlternation
      peakSampleInds = rawPeaks.loc(dataTrace(rawPeaks.loc) > lowThreshold);
      highPeakSampleInds = rawPeaks.loc(dataTrace(rawPeaks.loc) > highThreshold);
      lowPeakSampleInds = rawPeaks.loc(dataTrace(rawPeaks.loc) > lowThreshold & dataTrace(rawPeaks.loc) < highThreshold);
      frameCountDiff = length(highPeakSampleInds) - length(lowPeakSampleInds);
      % with correct alternation, smallest interval between subsequent high peaks must tbe greater than longest frame; same for low peaks 
      frameAlternationError = (min(min(diff(highPeakSampleInds)),min(diff(lowPeakSampleInds))) > max(diff(peakSampleInds)));
      notValid = (abs(frameCountDiff) > 1 || frameAlternationError);
      if notValid
        if abs(frameCountDiff) > 1
          tmp = input(sprintf('Frame count mismatch: Found %d bright-square frames and %d dark-square frames. Press ENTER to re-start manual check, or q to quit: ',...
            length(highPeakSampleInds),length(lowPeakSampleInds)));
        end
        if frameAlternationError
          tmp = input(sprintf('Found light-dark frame alternation error. Press ENTER to re-start manual check, c to continue, or q to quit: '));
        end
        while ~strcmp('tmp',{'','q','c'})
          tmp = input('Invalid input. Press ENTER to re-start manual check, c to continue, or q to quit: ');
        end
        if isempty(tmp)
          continue
        elseif strcmp(tmp,'q')
          return
        else
          break
        end
      end
    else
      notValid = 0;
    end
  end
end
if saveFigures
  saveFigure(outDir,sprintf('%s_Run%s.mat',calibFigFname,runNum), [], 1, 0, sprintf('%s, Run %s',dateSubject,runNum));
end


% rebuild the peaks with final threshold values
peakSampleInds = rawPeaks.loc(dataTrace(rawPeaks.loc) > lowThreshold);
if cleanPeaks
  peakSampleIndsClean = zeros(size(peakSampleInds));
  cleanPeaks_i = 1;
  lookAhead = 1;
  for peak_i = 1:length(peakSampleInds)
    if lookAhead > 1
      lookAhead = lookAhead - 1;
      continue
    end
    peakToKeep = peak_i;
    peakToKeepLevel = dataTrace(peakSampleInds(peak_i));
    while peak_i+lookAhead < length(peakSampleInds) && min(dataTrace(peakSampleInds(peak_i):peakSampleInds(peak_i+lookAhead))) > lowThreshold
      newPeakToKeepLevel = max(peakToKeepLevel,dataTrace(peakSampleInds(peak_i+lookAhead)));
      peakToKeep = peakToKeep*(newPeakToKeepLevel == peakToKeepLevel) + (peak_i+lookAhead)*(newPeakToKeepLevel ~= peakToKeepLevel);
      lookAhead = lookAhead + 1;
    end
    peakSampleIndsClean(cleanPeaks_i) = peakSampleInds(peakToKeep);
    cleanPeaks_i = cleanPeaks_i + 1;
  end
  peakSampleInds = peakSampleIndsClean(peakSampleIndsClean > 0);
end

if useRisingEdge
  peakSampleInds = find(diff(dataTrace > lowThreshold) == 1);
end

if numLevels == 2
  highPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > highThreshold);
  lowPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > lowThreshold & dataTrace(peakSampleInds) < highThreshold);
end

if numLevels == 3
  highPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > highThreshold);
  midPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > midThreshold & dataTrace(peakSampleInds) < highThreshold); 
  lowPeakSampleInds = peakSampleInds(dataTrace(peakSampleInds) > lowThreshold & dataTrace(peakSampleInds) < midThreshold); 
end

% find peak level transitions, and convert from sample inds to time in ms 
strobeTriggers.all = (1000*peakSampleInds/samplingFreq) - peakTimeOffset;
if numLevels == 2
  lowToHighInds = peakSampleInds(logical(horzcat(false, diff((dataTrace(peakSampleInds) > highThreshold)) == 1)));
  highToLowInds = peakSampleInds(logical(horzcat(false, diff((dataTrace(peakSampleInds) < highThreshold)) == 1)));
  
  strobeTriggers.low = (1000*lowPeakSampleInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.high = (1000*highPeakSampleInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.lowToHigh = (1000*lowToHighInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.highToLow = (1000*highToLowInds/samplingFreq) - peakTimeOffset;
end
if numLevels == 3
  
  toHighInds = peakSampleInds(logical(horzcat(false, diff((dataTrace(peakSampleInds) > highThreshold)) == 1)));
  toLowInds = peakSampleInds(logical(horzcat(false, diff((dataTrace(peakSampleInds) < midThreshold)) == 1)));
  toMidInds = peakSampleInds(logical(horzcat(false, diff((dataTrace(peakSampleInds) < highThreshold)) == 1) &...
    (dataTrace(peakSampleInds) > midThreshold) | horzcat(false, diff((dataTrace(peakSampleInds) > lowThreshold)) == 1) &...
    (dataTrace(peakSampleInds) < highThreshold)));
  
  fromHighInds = peakSampleInds(logical(horzcat(false, diff((dataTrace(peakSampleInds) < highThreshold)) == 1)));
  fromLowInds = peakSampleInds(logical(horzcat(false, diff((dataTrace(peakSampleInds) > midThreshold)) == 1)));
  
  fromMidInds = setdiff(union(toHighInds, toLowInds), union(fromHighInds, fromLowInds));
  
  lowToHighInds = intersect(fromLowInds, toHighInds);
  highToLowInds = intersect(fromHighInds, toLowInds);
  midToHighInds = intersect(fromMidInds, toHighInds);
  highToMidInds = intersect(fromHighInds, toMidInds);
  lowToMidInds = intersect(fromLowInds, toMidInds);
  midToLowInds = intersect(fromMidInds, toLowInds);
  
  strobeTriggers.low = (1000*lowPeakSampleInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.mid = (1000*midPeakSampleInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.high = (1000*highPeakSampleInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.lowToHigh = (1000*lowToHighInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.highToLow = (1000*highToLowInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.lowToMid = (1000*lowToMidInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.midToLow = (1000*midToLowInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.midToHigh = (1000*midToHighInds/samplingFreq) - peakTimeOffset;
  strobeTriggers.highToMid = (1000*highToMidInds/samplingFreq) - peakTimeOffset;
end  

% plot the peak levels and computed switch times
figure('Name','Strobe Level Outputs','NumberTitle','off');;
h = plot(strobeTriggers.all,ones(size(strobeTriggers.all)),'linestyle','none','marker','o','color','b');
hold on
% note: build handle and label arrays for legend to defend against Matlab bug that leads to incorrect legends when some data series are empty
handlesForLegend = h;
labelsForLegend = {'all'};
if numLevels == 2
  h = plot(strobeTriggers.high,4*ones(size(strobeTriggers.high)),'linestyle','none','marker','*','color','b');
  if ~isempty(strobeTriggers.high)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'high');
  end
  h = plot(strobeTriggers.low,3*ones(size(strobeTriggers.low)),'linestyle','none','marker','x','color','b');
  if ~isempty(strobeTriggers.low)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'low');
  end
  h = plot(strobeTriggers.highToLow,3.5*ones(size(strobeTriggers.highToLow)),'linestyle','none','marker','o','color','r');
  if ~isempty(strobeTriggers.highToLow)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'highToLow');
  end
  h = plot(strobeTriggers.lowToHigh,3.5*ones(size(strobeTriggers.lowToHigh)),'linestyle','none','marker','o','color','g');
  if ~isempty(strobeTriggers.lowToHigh)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'lowToHigh');
  end
  legend(handlesForLegend, labelsForLegend);
end
if numLevels == 3
  h = plot(strobeTriggers.high,5*ones(size(strobeTriggers.high)),'linestyle','none','marker','*','color','b');
  if ~isempty(strobeTriggers.high)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'high');
  end
  h = plot(strobeTriggers.mid,4*ones(size(strobeTriggers.mid)),'linestyle','none','marker','s','color','b');
  if ~isempty(strobeTriggers.mid)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'mid');
  end
  h = plot(strobeTriggers.low,3*ones(size(strobeTriggers.low)),'linestyle','none','marker','x','color','b');
  if ~isempty(strobeTriggers.low)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'low');
  end
  h = plot(strobeTriggers.highToLow,4*ones(size(strobeTriggers.highToLow)),'linestyle','none','marker','o','color','r');
  if ~isempty(strobeTriggers.highToLow)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'highToLow');
  end
  h = plot(strobeTriggers.lowToHigh,4*ones(size(strobeTriggers.lowToHigh)),'linestyle','none','marker','o','color','g');
  if ~isempty(strobeTriggers.lowToHigh)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'lowToHigh');
  end
  h = plot(strobeTriggers.highToMid,4.5*ones(size(strobeTriggers.highToMid)),'linestyle','none','marker','o','color','y');
  if ~isempty(strobeTriggers.highToMid)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'highToMid');
  end
  h = plot(strobeTriggers.midToHigh,4.5*ones(size(strobeTriggers.midToHigh)),'linestyle','none','marker','o','color','m');
  if ~isempty(strobeTriggers.midToHigh)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'midToHigh');
  end
  h = plot(strobeTriggers.midToLow,3.5*ones(size(strobeTriggers.midToLow)),'linestyle','none','marker','o','color','k');
  if ~isempty(strobeTriggers.midToLow)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'midToLow');
  end
  h = plot(strobeTriggers.lowToMid,3.5*ones(size(strobeTriggers.lowToMid)),'linestyle','none','marker','o','color','c');
  if ~isempty(strobeTriggers.lowToMid)
    handlesForLegend = horzcat(handlesForLegend,h);
    labelsForLegend = horzcat(labelsForLegend,'lowToMid');
  end
  legend(handlesForLegend, labelsForLegend);
end
if displayStats
  fprintf('\n\n\nStrobe timing statistics\n\n')
  fprintf('Mean inter-strobe interval: %.3f ms (%.3f Hz)\n',mean(diff(strobeTriggers.all)), 1000/(mean(diff(strobeTriggers.all))));
  fprintf('Max inter-strobe interval: %.3f ms (%.3f Hz)\n',max(diff(strobeTriggers.all)), 1000/(max(diff(strobeTriggers.all))));
  fprintf('Min inter-strobe interval: %.3f ms (%.3f Hz)\n',min(diff(strobeTriggers.all)), 1000/(min(diff(strobeTriggers.all))));
  fprintf('Mode inter-strobe interval: %.3f ms (%.3f Hz)\n',mode(diff(strobeTriggers.all)), 1000/(mode(diff(strobeTriggers.all))));
  fprintf('Inter-strobe interval std: %.3f ms\n',std(diff(strobeTriggers.all)));
end
if saveFigures
  saveFigure(outDir,sprintf('%s_Run%s.mat',triggersFigFname,runNum), [], 1, 0, sprintf('%s, Run %s',dateSubject,runNum));
end

if saveCalibFile 
  save(outputCalibrationFile, 'strobeTriggers');
end  
end

