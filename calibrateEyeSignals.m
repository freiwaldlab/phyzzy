function [ eyeX,eyeY,eyeD ] = calibrateEyeSignals( eyeX,eyeY,eyeD,taskData,params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
calibMethod = params.method;
gainX = params.gainX;
gainY = params.gainY;
flipY = params.flipY;
flipX = params.flipX;
if strcmp(calibMethod,'auto')
  fixInInds = zeros(size(eyeX));
  transitionInInds = zeros(size(eyeX));
  transitionOutInds = zeros(size(eyeX));
  for fix_i = 1:length(taskData.fixationInTimes-1)
    assert(all(taskData.fixationOutTimes > taskData.fixationInTimes));
    fixInInds(ceil(taskData.fixationInTimes(fix_i)):floor(taskData.fixationOutTimes(fix_i))) = 1;
    transitionInInds(ceil(taskData.fixationInTimes(fix_i))) = 1;
    transitionOutInds(floor(taskData.fixationOutTimes(fix_i))) = 1;
  end
  eyeX = gainX*((eyeX - mean(eyeX(fixInInds == 1)))*(-1*flipX + ~flipX));
  eyeY = gainY*((eyeY - mean(eyeY(fixInInds == 1)))*(-1*flipY + ~flipY));
  eyeD = eyeD - mean(eyeD(fixInInds == 1));
  fixInX = eyeX(fixInInds == 1);
  fixInY = eyeY(fixInInds == 1);
  fixInD = eyeD(fixInInds == 1);
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
  %histogram2(horzcat(fixInX',fixInY'),[50 50],'Normalization','probability');
  [n,c] = hist3(horzcat(fixInX',fixInY'),[50 50]);
  imagesc(c{1},c{2},n);
  a2 = subplot(1,2,2);
  %histogram2(horzcat(eyeX(~fixInInds)',eyeY(~fixInInds)'),[50 50],'Normalization','probability');
  [n,c] = hist3(horzcat(eyeX(~fixInInds)',eyeY(~fixInInds)'),[50 50]);
  imagesc(c{1},c{2},n);
  %linkaxes([a1 a2])
end

end

