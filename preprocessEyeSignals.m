function [ analogInData ] = preprocessEyeSignals( analogInData,taskData,params )
%preprocessEyeSignals applies gain, flip, and offset calibration to raw eye data
%   also has methods for computing offsets, and visualizing results
%   Inputs:
%     - analogInData: an nChannels x nTimepoints numeric array
%     - taskData: struct returned by preprocessLogFile
%     - params: struct with fields
%       - method: 
%       - eyeXChannelInd: int, row index into analogInData
%       - eyeYChannelInd: int, row index into analogInData
%       - eyeDChannelInd: int, row index into analogInData, (optional?)
%       - gainX: float, optional if calibMethod is fromFile
%       - gainY: float, optional if calibMethod is fromFile
%       - flipX: bool, optional if calibMethod is fromFile
%       - flipY: bool, optional if calibMethod is fromFile
%       - offsetX: float, applied after gain, optional unless calibMethod is hardcodeZero
%       - offsetY: float, applied after gain, optional unless calibMethod is hardcodeZero
%       - fixOutLag: 
%       - minFixZeroTime:
if ~params.needEyeCal
  return
end

calibMethod = params.method;
switch calibMethod
  case 'fromFile'
    load(params.calFile);
  case  {'hardcodeZero','monkeyLogic'}
    offsetX = taskData.eyeCal.origin(1);
    offsetY = taskData.eyeCal.origin(2);
    gainX = taskData.eyeCal.gain(1);
    gainY = taskData.eyeCal.gain(2);
    flipY = params.flipY;
    flipX = params.flipX;
  case 'zeroEachFixation'
    minFixZeroTime = params.minFixZeroTime; %Don't know what this means
    gainX = taskData.eyeCal.gain(1);
    gainY = taskData.eyeCal.gain(2);
    flipY = params.flipY;
    flipX = params.flipX;
end

if isfield(params, 'fixOutLag')
  fixOutLag = params.fixOutLag;
else
  fixOutLag = 10; %time in eye data samples; typically ms
end

eyeX = analogInData(params.eyeXChannelInd,:);
eyeY = analogInData(params.eyeYChannelInd,:);
eyeD = analogInData(params.eyeDChannelInd,:);

fixInInds = zeros(size(eyeX));
transitionInInds = zeros(size(eyeX));
transitionOutInds = zeros(size(eyeX));
transitionInInds(ceil(taskData.fixationInTimes)) = 1;
fixDurations = max(floor(taskData.fixationOutTimes) - taskData.fixationInTimes(1:length(taskData.fixationOutTimes)) - fixOutLag,0);
transitionOutInds(floor((taskData.fixationInTimes(1:length(taskData.fixationOutTimes))) + fixDurations)) = 1;

assert(all(taskData.fixationOutTimes > taskData.fixationInTimes(1:length(taskData.fixationOutTimes))),...
  'Invalid fixation time alignment: fixation must come in before it goes out');
assert(ismember(length(taskData.fixationInTimes) - length(taskData.fixationOutTimes), [0 1]),'Invalid fixation intervals: must have 0 or 1 more fix out times than fix in times'); 

for fix_i = 1:length(taskData.fixationOutTimes)
  fixInInds(ceil(taskData.fixationInTimes(fix_i)):(floor(taskData.fixationOutTimes(fix_i))-fixOutLag)) = 1;
end
if length(taskData.fixationOutTimes) > length(taskData.fixationInTimes)
  fixInInds(ceil(taskData.fixationInTimes(fix_i)):end) = 1;
end

switch calibMethod
  case 'monkeyLogic' %applies origin, then gain
    eyeX = eyeX - offsetX;
    eyeY = eyeY - offsetY;
    eyeX = gainX*eyeX*(-1*flipX + ~flipX);
    eyeY = gainY*eyeY*(-1*flipY + ~flipY);    
  case 'hardcodeZero' %applies gain, then origin
    eyeX = gainX*eyeX*(-1*flipX + ~flipX);
    eyeY = gainY*eyeY*(-1*flipY + ~flipY);
    eyeX = eyeX - offsetX;
    eyeY = eyeY - offsetY;
  case 'autoZeroSingle'
    eyeX = gainX*eyeX*(-1*flipX + ~flipX);
    eyeY = gainY*eyeY*(-1*flipY + ~flipY);
    eyeX = eyeX - mean(eyeX(fixInInds == 1));
    eyeY = eyeY - mean(eyeY(fixInInds == 1));
  case 'zeroEachFixation'
    eyeX = gainX*eyeX*(-1*flipX + ~flipX);
    eyeY = gainY*eyeY*(-1*flipY + ~flipY);
    for fix_i = 1:length(taskData.fixationInTimes)
      % find the new offset, if we don't have one, or if fixation is long enough
      % note: an alternative method here would be to find the first
      % long-enough interval, then back-apply its offset.
      if fix_i == 1 || (fix_i <= length(taskData.fixationOutTimes) && taskData.fixationOutTimes(fix_i) - taskData.fixationInTimes(fix_i) > minFixZeroTime)
        offsetX = mean(eyeX(ceil(taskData.fixationInTimes(fix_i)):(floor(taskData.fixationOutTimes(fix_i))-fixOutLag)));
        offsetY = mean(eyeY(ceil(taskData.fixationInTimes(fix_i)):(floor(taskData.fixationOutTimes(fix_i))-fixOutLag)));
      end
      % cover the case where final fix interval extends to end of run
      if fix_i > length(taskData.fixationOutTimes) && (fix_i == 1 || length(eyeX) - taskData.fixationInTimes(fix_i) > minFixZeroTime)
        offsetX = mean(eyeX(ceil(taskData.fixationInTimes(fix_i)):end));
        offsetY = mean(eyeY(ceil(taskData.fixationInTimes(fix_i)):end));
      end
      % apply the offset up to the next fixation window, or end, if last window
      if fix_i < length(taskData.fixationInTimes)
        eyeX(ceil(taskData.fixationInTimes(fix_i)):floor(taskData.fixationInTimes(fix_i+1))) =  ...
          eyeX(ceil(taskData.fixationInTimes(fix_i)):floor(taskData.fixationInTimes(fix_i+1))) - offsetX;
        eyeY(ceil(taskData.fixationInTimes(fix_i)):floor(taskData.fixationInTimes(fix_i+1))) =  ...
          eyeY(ceil(taskData.fixationInTimes(fix_i)):floor(taskData.fixationInTimes(fix_i+1))) - offsetY;
      else
        eyeX(ceil(taskData.fixationInTimes(fix_i)):end) = eyeX(ceil(taskData.fixationInTimes(fix_i)):end) - offsetX;
        eyeY(ceil(taskData.fixationInTimes(fix_i)):end) = eyeY(ceil(taskData.fixationInTimes(fix_i)):end) - offsetY;
      end
    end
end

fixInX = eyeX(fixInInds == 1);
fixInY = eyeY(fixInInds == 1);
fixInD = eyeD(fixInInds == 1);

%% monkeyLogic Shift
% If for whatever reason eye signal on Blackrock is unreliable, create
% artificial vectors from monkeyLogic's eye data, placing it at the event
% start times on Blackrock.
if std(eyeX) < 1 || std(eyeY) < 1
  if isfield(params, 'monkeyLogicShift') && params.monkeyLogicShift
    plotShift = 1; %If you want to plot what this part does, switch to 1.
    disp('Low variance in eye signal suggests problem with signal collection - using monkeyLogic eye signal instead')
    %Move MonkeyLogic trial start times using model based on disparities in
    %markers.
    trialStartTimesBlk = predict(taskData.logVsBlkModel, taskData.trialStartTimesMkl);
    
    %initialize arrays
    [eyeXBLK, eyeYBLK] = deal(nan(size(eyeY)));
    samPerMS = 1;
    for ii = 1:length(trialStartTimesBlk) %for every trial start
      %Find the spot in the vector to go (and end) - Mkl times already in ms.
      startInd = round(samPerMS*trialStartTimesBlk(ii)); %POSSIBLE SOURCE OF ISSUE
      if startInd == 0
        startInd = 1;
      end
      endInd = startInd + length(taskData.eyeData(ii).Eye(:,1)) - 1;
      
      %Place the data in the vector
      eyeXBLK(startInd:endInd) = taskData.eyeData(ii).Eye(:,1);
      eyeYBLK(startInd:endInd) = taskData.eyeData(ii).Eye(:,2);
    end
    
    if std(eyeX) < 1
      eyeCorrect = eyeY;
      eyeCorrectBLK = eyeYBLK;
      eyeWrong = eyeXBLK;
      plotTitle = 'eyeY-based model creation, used to make eyeX';
    else
      eyeCorrect = eyeX;
      eyeCorrectBLK = eyeXBLK;
      eyeWrong = eyeYBLK;
      plotTitle = 'eyeX-based model creation, used to make eyeY';
    end
    
    if plotShift
      plotStartInd = round(length(eyeY)/10);
      plotEndInd = plotStartInd + min(1e5, plotStartInd/2);
      figure()
      plot(eyeCorrect(plotStartInd:plotEndInd) , 'k')
      hold on
      plot(eyeCorrectBLK(plotStartInd:plotEndInd), 'r')
    end
    
    %adjust with finddelay.
    eyeCorrectBLK(isnan(eyeCorrectBLK)) = 0; %finddelay doesn't like NaN.
    lag = finddelay(eyeCorrectBLK, eyeY);
    eyeCorrectBLK = padarray(eyeCorrectBLK, [0,lag], 0, 'pre');
    eyeCorrectBLK = eyeCorrectBLK(1:end-lag);
    
    %Create a better fit by using an lm
    eyeCorrectBLK(eyeCorrectBLK == 0) = nan;
    eyeCorrectBLKModel = eyeCorrectBLK;
    eyeCorrectModel = eyeCorrect;
    
    logVsBlkModel = fitlm(eyeCorrectBLKModel, eyeCorrectModel);
    eyeCorrectBLKModel = predict(logVsBlkModel, eyeCorrectBLKModel');
    
    if plotShift
      %Add results to existing plot
      plot(eyeCorrectBLKModel(plotStartInd:plotEndInd), 'b')
      title(plotTitle);
      legend('Blackrock Signal','Reconstructed signal','Reconstructed signal, shifted and fit')
    end
    
    %Assuming the fit is good, now re-create the missing array.
    eyeWrong = padarray(eyeWrong, [0,lag], 0, 'pre');
    eyeWrong = eyeWrong(1:end-lag);
    eyeWrongModel = predict(logVsBlkModel, eyeWrong');
    
    %and overwrite the Blackrock values.
    if std(eyeX) < 1
      eyeY = eyeCorrectBLKModel';
      eyeX = eyeWrongModel';
    else
      eyeY = eyeWrongModel';
      eyeX = eyeCorrectBLKModel';
    end
    
    %Refilter the data based on some of the params of analogInParams. (What
    %follows is copied from preprocessAnalogIn, modified to be appropriate
    %here.
    reFilterEyes = [eyeX; eyeY];
    reFilterEyes(isnan(reFilterEyes)) = 0;
    filterPad = params.analogInParams.filterPad;
    reFilterEyesDecPadded = zeros(size(reFilterEyes,1),ceil(size(reFilterEyes,2)/(params.analogInParams.decimateFactorPass1*params.analogInParams.decimateFactorPass2))+2*filterPad);
    
    for i = 1:size(reFilterEyes,1)
      reFilterEyesDecPadded(i,filterPad+1:end-filterPad) = decimate(decimate(reFilterEyes(i,:),params.analogInParams.decimateFactorPass1),params.analogInParams.decimateFactorPass2);
      %analogInDataDecPadded(i,1:filterPad) = analogInDataDecPadded(i,filterPad+1)*analogInDataDecPadded(i,1:filterPad);
      %analogInDataDecPadded(i,end-(filterPad-1):end) = analogInDataDecPadded(i,end-filterPad)*analogInDataDecPadded(i,end-(filterPad-1):end);
      analogInFilter = params.analogInParams.filters{i};
      if isa(analogInFilter,'digitalFilter')
        disp('using digital filter object');
        reFilterEyesDecPadded(i,:) = filtfilt(analogInFilter, reFilterEyesDecPadded(i,:));
      elseif length(analogInFilter) == 2
        reFilterEyesDecPadded(i,:) = filtfilt(analogInFilter(1),analogInFilter(2),reFilterEyesDecPadded(i,:));
      end
      if plotShift && (isa(analogInFilter,'digitalFilter')) && i == 2
        plot(reFilterEyesDecPadded(i,filterPad+plotStartInd:filterPad+plotEndInd),'color','g');
        legend('Blackrock Signal','Reconstructed signal','Recon signal, shifted and fit','Recon signal, shifted, fit, and filt')
        drawnow;
      end
    end
    reFilterEyes = reFilterEyesDecPadded(:,filterPad+1:end-filterPad);
    eyeX = reFilterEyes(1,:);
    eyeY = reFilterEyes(2,:);
    
  else
    error('Low variance in eye signal suggests problem with signal collection - suggest using monkeyLogic eye signals (set params.monkeyLogicShift to 1)')
  end
elseif std(eyeX) < .1 && std(eyeY) < .1
  error('Low variance in both eye signals suggests problem with signal collection. Blackrock inputs not plugged in')
end

%% Plots and Output

if params.makePlots
  figure('Name','Fix In Scatter Plot','NumberTitle','off')
  scatter(fixInX,fixInY,1,'b');
  hold on
  scatter(eyeX(~fixInInds),eyeY(~fixInInds),1,'r')
  figure('Name','Fixation transitions','NumberTitle','off')
  scatter(eyeX(transitionInInds == 1),eyeY(transitionInInds == 1),20,'b');
  hold on
  scatter(eyeX(transitionOutInds == 1),eyeY(transitionOutInds == 1),20,'r');
  figure('Name','Eye dilation vs Fixation Histogram','NumberTitle','off')
  histogram(fixInD,100,'Normalization','pdf');
  hold on
  histogram(eyeD(~fixInInds),100,'Normalization','pdf');
  title('eyeD')
  legend({'In','Out'});
  figure('Name','Eye Signal Probability Distribution','NumberTitle','off')
  a1 = subplot(1,2,1);
  histogram2(fixInX',fixInY',[50 50],'Normalization','probability','DisplayStyle','tile','BinMethod','integer');
  a2 = subplot(1,2,2);
  histogram2(eyeX(~fixInInds)',eyeY(~fixInInds)',[50 50],'Normalization','probability','DisplayStyle','tile','BinMethod','integer');
end

analogInData(params.eyeXChannelInd,:) = eyeX;
analogInData(params.eyeYChannelInd,:) = eyeY;
analogInData(params.eyeDChannelInd,:) = eyeD;

end

