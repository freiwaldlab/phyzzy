function [ diodeTriggers ] = preprocessPhotodiodeStrobe( photodiodeFilename, params )
%preprocessPhotodiodeStrobe extracts frame times from a photodiode signal 
%   Inputs:
%   - analogInData: data array returned by preprocessAnalogIn (channels x samples)
%   - params: struct with the fields
%     - needPhotodiode: return immediately if 0 (type: logical)
%     - frameTriggerChannel: this is the index of the frame trigger photodiode channel in the list of analog in
%                           channels specified in the paramFile (type: int)
%     - stimulusTriggerChannel: this is the index of the stimulus trigger photodiode channel in the list of analog in
%                           channels specified in the paramFile (type: int)
%     - centerCornerOffset: this is the offset, in ms, of the photodiode
%                           trigger from the screen-center frame time
%                           (note: will be subtracted, so should be >0 for
%                           diode at bottom of screen, and <0 for diode at top of screen) 
%                           (type: float)
%     - levelCalibType: (type: string) options are:
%       - hardcode
%       -
%       -
%
%
%     - numLevels: (type: int)
%     - hardcodeFromFile: (type: logical)
%     - saveCalibFile: (type: logical)
%     - inputCalibrationFile (required only if hardcodeFromFile): (type: string)
%     - outputCalibrationFile (required only if saveCalibFile): (type:string)
if ~params.needPhotodiode
  diodeTriggers = [];
  return
end

% load photodiode data
assert(logical(exist(photodiodeFilename,'file')),'The photodiode data file you requested does not exist.');
tmp = openNSx(photodiodeFilename,'report','read');
photodiodeHeader = tmp.MetaTags;
sampleRate = photodiodeHeader.SamplingFreq;
params.frameTriggerChannel = find(photodiodeHeader.ChannelID == params.frameTriggerChannel);
frameTriggerDiodeTrace = tmp.Data(params.frameTriggerChannel,:);
stimulusTriggerDiodeTrace = [];
if ~isempty(params.stimulusTriggerChannel)
  params.stimulusTriggerChannel = find(photodiodeHeader.ChannelID == params.stimulusTriggerChannel);
  stimulusTriggerDiodeTrace = tmp.Data(params.stimulusTriggerChannel,:);
end
clear tmp

centerCornerOffset = params.centerCornerOffset; 
frameNumCeiling = ceil(length(frameTriggerDiodeTrace)*params.frameRate/sampleRate);


for diode_i = 1:2
  if diode_i == 1
    diodeTrace = frameTriggerDiodeTrace;
  else
    if isempty(stimulusTriggerDiodeTrace)
      break
    end
    diodeTrace = stimulusTriggerDiodeTrace;
  end
        
  % Preprocessing logic:
  % 1) Find peaks in the inverted diode trace. A peak occurs when the cathode ray passes the diode,
  %    even on black pixels. However, additional small-prominence peaks occur due to noise
  % 2) Sort peaks by their height: since both the the signal and its
  %    second derivative are high at true peaks, they should be the largest
  diodeTrace = -1*(diodeTrace - mean(diodeTrace)); %this makes frames peaks
  rawPeaks = findpeaks(diodeTrace);
  rawPeaksSorted = sort(diodeTrace(rawPeaks.loc),'descend');
  rawPeaksSorted = rawPeaksSorted(1:frameNumCeiling);

  % additional todo: add frame count and alternation checks for manual levels

  if strcmp(params.levelCalibType,'hardcode')  
    if params.hardcodeFromFile
      switch params.numLevels
        case 1
          load(params.inputCalibrationFile,'highThreshold');
        case 2
          load(params.inputCalibrationFile,'highThreshold','lowThreshold');
        case 3
          load(params.inputCalibrationFile,'highThreshold','midThreshold','lowThreshold');
      end
    else
      switch params.numLevels
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

  if any(strcmp(params.levelCalibType, {'auto','autoAndPlot','autoAndCheck'}))
    % Auto-calibration logic: 
    %  The largest jumps in the plot of sorted peak heights should occur between true peak types,  
    %  or between the smallest peak type and the largest noise peak. So find and sort those jumps, 
    %  then place level cut-offs at the mid-point of the largest ones. 

    peakDiffs = diff(rawPeaksSorted);
    [~,peakDiffSortInds] = sort(peakDiffs);
    if ~(strcmp(params.levelCalibType,'hardcodeAndPlot') || strcmp(params.levelCalibType,'hardcodeAndCheck'))
      switch params.numLevels
        case 1
          highThresholdDiffInd = peakDiffSortInds(1);
          highThreshold = (rawPeaksSorted(highThresholdDiffInd + 1) - rawPeaksSorted(highThresholdDiffInd))/2;
        case 2
          highThresholdDiffInd = min(peakDiffSortInds(1:2)); %note: these end up on the high side of the threshold,i.e. the low index side
          lowThresholdDiffInd = max(peakDiffSortInds(1:2));
          highThreshold = (rawPeaksSorted(highThresholdDiffInd + 1) - rawPeaksSorted(highThresholdDiffInd))/2;
          lowThreshold = (rawPeaksSorted(lowThresholdDiffInd + 1) - rawPeaksSorted(lowThresholdDiffInd))/2;
        case 3
          highThresholdDiffInd = min(peakDiffSortInds(1:3)); %note: these end up on the high side of the threshold,i.e. the low index side
          lowThresholdDiffInd = max(peakDiffSortInds(1:3));
          midThresholdDiffInd = setdiff(peakDiffSortInds(1:3),[highThresholdDiffInd, lowThresholdDiffInd]);
          highThreshold = (rawPeaksSorted(highThresholdDiffInd + 1) - rawPeaksSorted(highThresholdDiffInd))/2;
          midThreshold = (rawPeaksSorted(midThresholdDiffInd + 1) - rawPeaksSorted(midThresholdDiffInd))/2;
          lowThreshold = (rawPeaksSorted(lowThresholdDiffInd + 1) - rawPeaksSorted(lowThresholdDiffInd))/2;
      end
    end
  end
  
  if params.numLevels == 1
    frameSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > highThreshold);
    highPeakSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > highThreshold);
  end
  % if two levels, check that frames alternate appropriately
  if params.numLevels == 2
    frameSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > lowThreshold);
    highPeakSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > highThreshold);
    lowPeakSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > lowThreshold & diodeTrace(rawPeaks.loc) < highThreshold);
    frameCountDiff = length(highPeakSampleInds) - length(lowPeakSampleInds);
    % with correct alternation, smallest interval between subsequent high peaks must tbe greater than longest frame; same for low peaks 
    frameAlternationError = (min(min(diff(highPeakSampleInds)),min(diff(lowPeakSampleInds))) > max(diff(frameSampleInds)));
  end

  %make plots, if requested, or if auto calib failed
  if any(strcmp(params.levelCalibType,{'autoAndPlot','autoAndCheck','hardcodeAndPlot','hardcodeAndCheck'})) || (params.numLevels == 2 && (abs(frameCountDiff) > 1 || frameAlternationError))
    if params.numLevels == 2 && abs(frameCountDiff) > 1
      fprintf('Frame count mismatch: Found %d bright-square frames and %d dark-square frames. Starting manual check.\n',...
        length(highPeakSampleInds),length(lowPeakSampleInds));
    end
    if params.numLevels == 2 && frameAlternationError
      fprintf('Found light-dark frame alternation error. Starting manual check.\n');
    end
    if any(strcmp(params.levelCalibType, {'hardcodeAndPlot','hardcodeAndCheck','hardcode'}))
      calibSuffix = 'hard code';
    else
      calibSuffix = 'auto';
    end

    figure();
    subplot(1,2,1);
    framesToPlot = 10;  %magic number; should move to top
    sampleRange = ceil(frameSampleInds(1) - 5*sampleRate/params.frameRate):ceil(frameSampleInds(1) + framesToPlot*sampleRate/params.frameRate);
    sampleRange = sampleRange(sampleRange > 0);
    plot(sampleRange,diodeTrace(sampleRange));
    hold on
    h = get(gca,'xlim');
    plot([h(1) h(2)], [highThreshold highThreshold]);
    forLegend121 = {'raw trace',sprintf('high threshold %s',calibSuffix)}; %note: 121 refers to the subplot
    if params.numLevels == 3
      plot([h(1) h(2)], [midThreshold midThreshold]);
      forLegend121 = horzcat(forLegend121, {sprintf('mid threshold %s',calibSuffix)});
    end
    if params.numLevels > 1
      plot([h(1) h(2)], [lowThreshold lowThreshold]);
      forLegend121 = horzcat(forLegend121, {sprintf('low threshold %s',calibSuffix)});
    end
    plot(frameSampleInds(1:framesToPlot), diodeTrace(frameSampleInds(1:framesToPlot)));
    forLegend121 = horzcat(forLegend121, {'frameTimes'});
    legend(forLegend121);  
    subplot(1,2,2);
    plot(rawPeaksSorted);
    hold on
    h = get(gca,'xlim');
    plot([h(1) h(2)], [highThreshold highThreshold]);
    forLegend122 = {'peak levels',sprintf('high threshold %s',calibSuffix)}; %note: 122 refers to the subplot
    if params.numLevels == 3
      plot([h(1) h(2)], [midThreshold midThreshold]);
      forLegend122 = horzcat(forLegend122, {sprintf('mid threshold %s',calibSuffix)});
    end
    if params.numLevels > 1
      plot([h(1) h(2)], [lowThreshold lowThreshold]);
      forLegend122 = horzcat(forLegend122, {sprintf('low threshold %s',calibSuffix)});
    end
    legend(forLegend122);
    drawnow();

    if any(strcmp(params.levelCalibType,{'autoAndCheck','hardcodeAndCheck'})) || params.numLevels == 2
      tmp = input('Accept high threshold? Enter y or n: ','s');
      if strcmp(tmp,'n')
        tmp = input('Enter new high threshold, or ''quit'' to continue without photodiode sync:','s');
        while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
          tmp = input('Invalid input: enter new high threshold, or ''quit'' to continue without photodiode sync:','s');
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
      if params.numLevels == 3
        tmp = input('Accept mid threshold? Enter y or n: ','s');
        if strcmp(tmp,'n')
          tmp = input('Enter new mid threshold, or ''quit'' to continue without photodiode sync:','s');
          while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
            tmp = input('Invalid input: enter new mid threshold, or ''quit'' to continue without photodiode sync:','s');
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
      end
      if params.numLevels > 1
        tmp = input('Accept low threshold? Enter y or n: ','s');
        if strcmp(tmp,'n')
          tmp = input('Enter new low threshold, or ''quit'' to continue without photodiode sync:','s');
          while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
            tmp = input('Invalid input: enter new low threshold, or ''quit'' to continue without photodiode sync:','s');
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
            legend(forLegend122);
          end
        end
      end
    end
  end

  % if two levels, check that frames alternate appropriately
  if params.numLevels == 2
    frameSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > lowThreshold);
    highPeakSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > highThreshold);
    lowPeakSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > lowThreshold & diodeTrace(rawPeaks.loc) < highThreshold);
    frameCountDiff = length(highPeakSampleInds) - length(lowPeakSampleInds);
    % with correct alternation, smallest interval between subsequent high peaks must tbe greater than longest frame; same for low peaks 
    frameAlternationError = (min(min(diff(highPeakSampleInds)),min(diff(lowPeakSampleInds))) > max(diff(frameSampleInds)));
    notValid = (abs(frameCountDiff) > 1 || frameAlternationError);
    if notValid
      if abs(frameCountDiff) > 1
        tmp = input(sprintf('Frame count mismatch: Found %d bright-square frames and %d dark-square frames. Press ENTER to start manual check, or q to quit.',...
          length(highPeakSampleInds),length(lowPeakSampleInds)));
      end
      if frameAlternationError
        tmp = input(sprintf('Found light-dark frame alternation error. Press ENTER to start manual check, or q to quit.'));
      end
      while ~strcmp('tmp',{'','q'})
        tmp = input('Invalid input. Press ENTER to re-start manual check, or q to quit.');
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

  if strcmp(params.levelCalibType,'manual') || notValid
    notValid = 1;
    while notValid
      figure();
      subplot(1,2,1);
      framesToPlot = 10;  %magic number; should move to top
      sampleRange = ceil(length(diodeTrace)/2 - (framesToPlot/2)*sampleRate/params.frameRate):ceil(length(diodeTrace)/2 - (framesToPlot/2)*sampleRate/params.frameRate);
      plot(sampleRange,diodeTrace(sampleRange));
      forLegend121 = {'raw trace'};
      legend(forLegend121);
      subplot(1,2,2);
      plot(rawPeaksSorted);
      forLegend122 = {'peak levels'};
      legend(forLegend122);
      drawnow();
      tmp = input('Enter the high threshold, or ''quit'' to continue without photodiode sync:','s');
      while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
        tmp = input('Invalid input: enter the high threshold, or ''quit'' to continue without photodiode sync:','s');
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
      if params.numLevels == 3
        tmp = input('Enter the mid threshold, or ''quit'' to continue without photodiode sync:','s');
        while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
          tmp = input('Invalid input: enter the mid threshold, or ''quit'' to continue without photodiode sync:','s');
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
      if params.numLevels > 1
        tmp = input('Enter the low threshold, or ''quit'' to continue without photodiode sync:','s');
        while isempty(str2double(tmp)) && ~strcmp(tmp,'quit')
          tmp = input('Invalid input: enter the low threshold, or ''quit'' to continue without photodiode sync:','s');
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
          legend(forLegend122);
        end
      end
      % if two levels, check that frames alternate appropriately
      if params.numLevels == 2
        frameSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > lowThreshold);
        highPeakSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > highThreshold);
        lowPeakSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > lowThreshold & diodeTrace(rawPeaks.loc) < highThreshold);
        frameCountDiff = length(highPeakSampleInds) - length(lowPeakSampleInds);
        % with correct alternation, smallest interval between subsequent high peaks must tbe greater than longest frame; same for low peaks 
        frameAlternationError = (min(min(diff(highPeakSampleInds)),min(diff(lowPeakSampleInds))) > max(diff(frameSampleInds)));
        notValid = (abs(frameCountDiff) > 1 || frameAlternationError);
        if notValid
          if abs(frameCountDiff) > 1
            tmp = input(sprintf('Frame count mismatch: Found %d bright-square frames and %d dark-square frames. Press ENTER to re-start manual check, or q to quit.',...
              length(highPeakSampleInds),length(lowPeakSampleInds)));
          end
          if frameAlternationError
            tmp = input(sprintf('Found light-dark frame alternation error. Press ENTER to re-start manual check, or q to quit.'));
          end
          while ~strcmp('tmp',{'','q'})
            tmp = input('Invalid input. Press ENTER to re-start manual check, or q to quit.');
          end
          if isempty(tmp)
            continue
          else
            break
          end
        end
      else
        notValid = 0;
      end
    end
  end
  if params.saveCalibFig
    %saveFigure(outDir,sprintf(,runNum), figData, saveFig, exportFig, saveFigData, sprintf('%s, Run %s',dateSubject,runNum));
  end
  
  switch numLevels
    case 1
      error('numLevels = 1 not yet implemented');
    case 2
      frameSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > lowThreshold);
    case 3
      frameSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > lowThreshold);
      highPeakSampleInds = rawPeaks.loc(diodeTrace(rawPeaks.loc) > highThreshold);
  end 
  if diode_i == 1
    photodiodeFrameTimes = (1000*frameSampleInds/samplingFreq) - centerCornerOffset;
    if numLevels == 3
      photodiodeHighFrameTimes = (1000*highPeakSampleInds/samplingFreq) - centerCornerOffset;
    else
      photodiodeHighFrameTimes = [];
    end
  end
  if diode_i == 2
    if numLevels == 2
      highPeakIntervals = diff(highPeakSampleInds)/samplingFreq;
      highStimStartMask = (highPeakIntervals > (1.5/frameRate));
      highStimStartMask = horzcat(1,highStimStartMask);
      highStimStartTimes = 1000*highPeakSampleInds(highStimStartMask)/samplingFreq;
      %
      lowPeakSampleInds = setdiff(frameSampleInds,highPeakSampleInds);
      lowPeakIntervals = diff(lowPeakSampleInds)/samplingFreq;
      lowStimStartMask = (lowPeakIntervals > (1.5/frameRate));
      lowStimStartMask = horzcat(1,lowStimStartMask);
      lowStimStartTimes = 1000*lowPeakSampleInds(lowStimStartMask)/samplingFreq;
    end
    if numLevels == 3
      error('three-level stimulus trigger diode processing not yet implemented');
    end
  end
end
if isempty(stimulusTriggerDiodeTrace)
  highStimStartTimes = [];
  lowStimStartTimes = [];
end
diodeTriggers.frameTimes = photodiodeFrameTimes; 
diodeTriggers.highFrameTimes = photodiodeHighFrameTimes;
diodeTriggers.highStimStartTimes = highStimStartTimes;
diodeTriggers.lowStimStartTimes = lowStimStartTimes;
if params.saveCalibFile %todo: add save for second diode's calib
  switch params.numLevels
    case 1
      save(params.outputCalibrationFilename, 'highThreshold');
    case 2
      save(params.outputCalibrationFilename, 'highThreshold', 'lowThreshold');
    case 3
      save(params.outputCalibrationFilename, 'highThreshold', 'midThreshold', 'lowThreshold');
  end
end  
end

