function [ analogInByItem ] = alignAnalogIn( analogInData, alignPointsByItem, analogInChannels, alignParams)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if isempty(analogInData)
  analogInByItem = [];
  Output.VERBOSE('not aligning analog in data');
  return
end

samPerMS = alignParams.samPerMS;
samplesPreAlign = samPerMS*alignParams.msPreAlign;
samplesPostAlign = samPerMS*alignParams.msPostAlign;


analogInByItem = cell(length(alignPointsByItem),1); 
for item_i = 1:length(alignPointsByItem)  
  onsets = int32(alignPointsByItem{item_i});
  analogInArray = zeros(1,length(analogInChannels),length(onsets),samplesPreAlign+samplesPostAlign+1); %(1,channel,trial,sample)
  for trial_i = 1:length(onsets)
    analogInArray(1,:,trial_i,:) = analogInData(:,samPerMS*(onsets(trial_i)-samplesPreAlign):samPerMS*(onsets(trial_i)+samplesPostAlign));
  end
  analogInByItem{item_i} = analogInArray;
end
end

