function [ lfpByItem ] = alignLFP( lfpData, alignPointsByItem, lfpChannels, alignParams  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

samPerMS = alignParams.samPerMS;
DCSUB_SAM = alignParams.DCSUB_SAM;
samplesPreAlign = samPerMS*alignParams.msPreAlign;
samplesPostAlign = samPerMS*alignParams.msPostAlign;


lfpByItem = cell(length(alignPointsByItem),1);
for item_i = 1:length(alignPointsByItem)  
  onsets = int32(alignPointsByItem{item_i});
  lfpArray = zeros(1,length(lfpChannels),length(onsets),samplesPreAlign+samplesPostAlign+1); %(1,channel,trial,sample)
  for trial_i = 1:length(onsets)
    lfpArray(1,:,trial_i,:) = lfpData(:,samPerMS*(onsets(trial_i)-samplesPreAlign):samPerMS*(onsets(trial_i)+samplesPostAlign));
    if any(DCSUB_SAM(1,:))
      for ch_i = 1:size(lfpArray,2)
        lfpArray(1,ch_i,trial_i,:) = lfpArray(1,ch_i,trial_i,:) - mean(lfpArray(1,ch_i,trial_i,DCSUB_SAM(1,1):DCSUB_SAM(1,2)),4); %set the first 20ms mean equal across trials and stimuli
      end
    end
    if any(DCSUB_SAM(2,:))
      for ch_i = 1:size(lfpArray,2)
        endVal = mean(lfpArray(1,ch_i,trial_i,DCSUB_SAM(2,1):DCSUB_SAM(2,2)),4);
        lfpArray(1,ch_i,trial_i,:) = squeeze(lfpArray(1,ch_i,trial_i,:)) - linspace(0,endVal,size(lfpArray,4))'; %subtract the linear term to set final values to zero
      end
    end
  end
  lfpByItem{item_i} = lfpArray;
end
end

