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
if isfield(params, 'monkeyLogicShift') && params.monkeyLogicShift
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
    endInd = startInd + length(taskData.eyeData(ii).Eye(:,1)') - 1;
    
    %Place the data in the vector
    eyeXBLK(startInd:endInd) = taskData.eyeData(ii).Eye(:,1)';
    eyeYBLK(startInd:endInd) = taskData.eyeData(ii).Eye(:,2)';
  end
  
  figure()
  plot(eyeY, 'k')
  hold on
  plot(eyeYBLK, 'r')
  
  %adjust with finddelay.
  eyeYBLK(isnan(eyeYBLK)) = 0;
  lag = finddelay(eyeYBLK, eyeY);
  eyeYBLK = padarray(eyeYBLK, [0,lag], 0, 'pre');
  eyeYBLK = eyeYBLK(1:end-lag);
  
  %Create a better fit by using an lm
  eyeYModel = eyeY';
  eyeYBLK(eyeYBLK == 0) = nan;
  eyeYBLKModel = eyeYBLK';
  
  logVsBlkModel = fitlm(eyeYBLKModel, eyeYModel);
  eyeYBLKModel = predict(logVsBlkModel, eyeYBLKModel);
  plot(eyeYBLKModel, 'b')
  
  legend('Blackrock Y Signal','Reconstructed Y signal','Reconstructed and Fit Y signal')
  
  %Assuming the fit is good, now re-create the missing x array.
  eyeXBLK = padarray(eyeXBLK, [0,lag], 0, 'pre');
  eyeXBLK = eyeXBLK(1:end-lag);
  eyeXBLKModel = predict(logVsBlkModel, eyeXBLK');

  %and overwrite the Blackrock values.
  eyeY = eyeYBLKModel;
  eyeX = eyeXBLKModel;
  
end

%% Plots and Output

if params.makePlots
  figure()
  scatter(fixInX,fixInY,1,'b');
  hold on
  scatter(eyeX(~fixInInds),eyeY(~fixInInds),1,'r')
  figure()
  scatter(eyeX(transitionInInds == 1),eyeY(transitionInInds == 1),20,'b');
  hold on
  scatter(eyeX(transitionOutInds == 1),eyeY(transitionOutInds == 1),20,'r');
  figure()
  histogram(fixInD,100,'Normalization','pdf');
  hold on
  histogram(eyeD(~fixInInds),100,'Normalization','pdf');
  title('eyeD')
  legend({'In','Out'});
  figure()
  a1 = subplot(1,2,1);
  histogram2(fixInX',fixInY',[50 50],'Normalization','probability','DisplayStyle','tile','BinMethod','integer');
  a2 = subplot(1,2,2);
  histogram2(eyeX(~fixInInds)',eyeY(~fixInInds)',[50 50],'Normalization','probability','DisplayStyle','tile','BinMethod','integer');
end

analogInData(params.eyeXChannelInd,:) = eyeX;
analogInData(params.eyeYChannelInd,:) = eyeY;
analogInData(params.eyeDChannelInd,:) = eyeD;

end

