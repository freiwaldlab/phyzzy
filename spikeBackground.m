function [ pre, post ] = spikeBackground( spikesByChannel, alignPointsByItem, spikeChannels, params )
%should return full spike time cells; runAnalyses can deal wtih stats later
%   Detailed explanation goes here


spikesByItem = cell(length(alignPointsByItem),1);
psthEmptyByItem = cell(length(alignPointsByItem),1);
for channel_i = 1:length(spikeChannels)
  channelSpikes = spikesByChannel(channel_i);
  for unit_i = 1:length(units)
    tstamps = channelSpikes.times(channelSpikes.units == units(unit_i));
    unitSpikes(trial_i) = struct('times', tstamps(tstamps>=onsets(trial_i) - preAlign & tstamps <= onsets(trial_i)+postAlign) - onsets(trial_i) - refOffset);
    channelSpikes{unit_i} = itemUnitSpikes;
  end
end

end

