function [ spikeCounts, spikeRatesTrialAve, spikeRatesTrialErr] = spikeCounter( spikesByItem, pre, post )
%spikeCounter counts spikes within (pre, post) of the trial alignment time
%   input:
%   - spikesByItem, as returned by alignSpikes:
%   spikesByItem{stim}{channel}{unit}(trial).times
%   - pre and post in ms
%   returns:
%   - spikeCounts: cell array of per-trial counts and rates, indexed {channel}{unit}{item}.counts(trial) and .rates(trial).
%   - spikeCountTrialAve: trial ave spike rates, indexed {channel}(unit,item)
%   - spikeCountTrialErr: trial sterr spike rates, indexed {channel}(unit,item)
%   todo: 
%   - add mean spike time and spike time stdev, or gaussian fit?

numChannels = length(spikesByItem{1});
[spikeCounts, spikeRatesTrialAve, spikeRatesTrialErr] = deal(cell(numChannels,1));

for channel_i = 1:numChannels
  channelSpikeCounts = cell(length(spikesByItem{1}{channel_i}),1);
  [channelSpikeRatesTrialAve, channelSpikeRatesTrialErr] = deal(zeros(length(spikesByItem{1}{channel_i}),length(spikesByItem)));
  for unit_i = 1:length(spikesByItem{1}{channel_i})
    unitSpikeCounts = cell(length(spikesByItem),1);
    for item_i = 1:length(spikesByItem)
      [s.counts, s.rates] = deal(zeros(length(spikesByItem{item_i}{channel_i}{unit_i}),1)); 
      for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
        s.counts(trial_i) = sum(spikesByItem{item_i}{channel_i}{unit_i}(trial_i).times > pre & spikesByItem{item_i}{channel_i}{unit_i}(trial_i).times < post);
        s.rates(trial_i) = (1000/(post-pre))*s.counts(trial_i);
      end
      unitSpikeCounts{item_i} = s;
      channelSpikeRatesTrialAve(unit_i, item_i) = mean(s.rates);
      n = length(s.rates);
      channelSpikeRatesTrialErr(unit_i, item_i) = (std(s.rates)/sqrt(n))*sqrt(2/(n-1))*(gamma(n/2)/gamma((n-1)/2));
    end
    channelSpikeCounts{unit_i} = unitSpikeCounts;
  end
  spikeCounts{channel_i} = channelSpikeCounts;
  spikeRatesTrialAve{channel_i} = channelSpikeRatesTrialAve;
  spikeRatesTrialErr{channel_i} = channelSpikeRatesTrialErr;
end
end

