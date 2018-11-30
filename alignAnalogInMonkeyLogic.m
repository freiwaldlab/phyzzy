function [ analogInByItem ] = alignAnalogInMonkeyLogic( analogInData, alignPointsByItem, analogInChannels, alignParams  )
% Take MonkeyLogic AnalogIn time series, as broken up by Item, and passes
% this forward. 
%   Detailed explanation goes here
if isempty(analogInData)
  analogInByItem = [];
  Output.VERBOSE('not aligning analog in data');
  return
end

samPerMS = alignParams.samPerMS;
samplesPreAlign = samPerMS*alignParams.msPreAlign;
samplesPostAlign = samPerMS*alignParams.msPostAlign;


if isfield(alignParams,'useMonkeyLogicAnalog') && alignParams.useMonkeyLogicAnalog == 1
  %Load the Logfile again
  [data,MLConfig,TrialRecord] = mlread(alignParams.taskFilename);
  
  %Create the corret array of eye signals
  trueTrialArray = (TrialRecord.TrialErrors == 0);
  tmpEye = [data(trueTrialArray).AnalogData];
  tmpEye = rmfield(tmpEye, {'SampleInterval', 'EyeExtra','Joystick','Mouse','PhotoDiode','General','Button'});
  
  %Issue - Time stamps in monkeylogic data and those in the
  %alignPointsByItem are not in the same clock space. Best way to align
  %them might to find the lowest value in the unaligned data, add them 
  
  
end


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

