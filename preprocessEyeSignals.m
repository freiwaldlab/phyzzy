function [ analogInData ] = preprocessEyeSignals( analogInData,taskData,params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~params.needEyeCal
  return
end
calibMethod = params.method;
gainX = params.gainX;
gainY = params.gainY;
flipY = params.flipY;
flipX = params.flipX;
if strcmp(calibMethod,'hardcodeZero')
  offsetX = params.offsetX;
  offsetY = params.offsetY;
end
if strcmp(calibMethod, 'zeroEachFixation')
  minFixZeroTime = params.minFixZeroTime;
end

eyeX = analogInData(params.eyeXChannelInd,:);
eyeY = analogInData(params.eyeYChannelInd,:);
eyeD = analogInData(params.eyeDChannelInd,:);

fixInInds = zeros(size(eyeX));
transitionInInds = zeros(size(eyeX));
transitionOutInds = zeros(size(eyeX));
for fix_i = 1:length(taskData.fixationInTimes)  
  assert(all(taskData.fixationOutTimes > taskData.fixationInTimes(1:length(taskData.fixationOutTimes)))); 
    taskData.fixationInTimes=taskData.fixationInTimes(1:length(taskData.fixationOutTimes));
  fixInInds(ceil(taskData.fixationInTimes(fix_i)):(floor(taskData.fixationOutTimes(fix_i))-10)) = 1;
  transitionInInds(ceil(taskData.fixationInTimes(fix_i))) = 1;
  transitionOutInds(floor(taskData.fixationOutTimes(fix_i))-10) = 1;
end
eyeX = gainX*eyeX*(-1*flipX + ~flipX);
eyeY = gainY*eyeY*(-1*flipY + ~flipY);
if strcmp(calibMethod,'autoZeroSingle')
  eyeX = eyeX - mean(eyeX(fixInInds == 1));
  eyeY = eyeY - mean(eyeY(fixInInds == 1));
end
if strcmp(calibMethod,'zeroEachFixation')
  for fix_i = 1:length(taskData.fixationInTimes)
    if fix_i == 1 || (taskData.fixationOutTimes(fix_i) - taskData.fixationInTimes(fix_i)) > minFixZeroTime
      offsetX = mean(eyeX(ceil(taskData.fixationInTimes(fix_i)):(floor(taskData.fixationOutTimes(fix_i))-10)));
      offsetY = mean(eyeY(ceil(taskData.fixationInTimes(fix_i)):(floor(taskData.fixationOutTimes(fix_i))-10)));
    end
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
if strcmp(calibMethod,'hardcodeZero')
  eyeX = eyeX - offsetX;
  eyeY = eyeY - offsetY;
end
if strcmp(calibMethod, 'ninePoint')
  error('ninePoint eye calibration not yet implemented');
end
eyeD = eyeD - mean(eyeD(fixInInds == 1));
fixInX = eyeX(fixInInds == 1);
fixInY = eyeY(fixInInds == 1);
fixInD = eyeD(fixInInds == 1);
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
end

