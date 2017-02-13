function [ spikeCounts, spikeCountTrialAve, spikeCountTrialErr ] = spikeCounter( spikesByItem, pre, post )
%spikeCounter counts spikes within (pre, post) of the trial alignment time
%   input:
%   - spikesByItem, as returned by alignSpikes; pre and post in ms
%   returns:
%   - spikeCounts: cell array of per-trial rates, indexed {channel}{unit}{item}.counts(trial)
%   - spikeCountTrialAve: trial ave spike rates, indexed {channel}(unit,item)
%   - spikeCountTrialErr: trial sterr spike rates, indexed {channel}(unit,item)
%   todo: 
%   - add mean spike time and spike time stdev, or gaussian fit?
%   - todo: rename 'counts' variable to 'rates'


numChannels = length(spikesByItem{1});
spikeCounts = cell(numChannels);
spikeCountTrialAve = cell(numChannels);
spikeCountTrialErr = cell(numChannels);

for channel_i = 1:numChannels
  channelSpikeCounts = cell(length(spikesByItem{1}{channel_i}));
  channelSpikeCountsTrialAve = zeros(length(spikesByItem{1}{channel_i}),length(spikesByItem));
  channelSpikeCountsTrialErr = zeros(length(spikesByItem{1}{channel_i}),length(spikesByItem));
  for unit_i = 1:length(spikesByItem{1}{channel_i})
    unitSpikeCounts = cell(length(spikesByItem));
    for item_i = 1:length(spikesByItem)
      s.counts = zeros(length(spikesByItem{item_i}{channel_i}{unit_i}),1); 
      for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
        s.counts(trial_i) = (1000/(post-pre))*sum(spikesByItem{item_i}{channel_i}{unit_i}(trial_i).times > pre & spikesByItem{item_i}{channel_i}{unit_i}(trial_i).times < post);
      end
      unitSpikeCounts{channel_i}{unit_i}{item_i} = s;
      channelSpikeCountsTrialAve(unit_i, item_i) = mean(s.counts);
      channelSpikeCountsTrialErr(unit_i, item_i) = std(s.counts)/sqrt(length(s.counts));
    end
    channelSpikeCounts{channel_i}{unit_i} = unitSpikeCounts;
  end
  spikeCounts{channel_i} = channelSpikeCounts;
  spikeCountTrialAve{channel_i} = channelSpikeCountsTrialAve;
  spikeCountTrialErr{channel_i} = channelSpikeCountsTrialErr;
end
end

