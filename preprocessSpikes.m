function [ spikesByChannel, pictureTriggers ] = preprocessSpikes( spikeFilename, params )
%UNTITLED5 Summary of this function goes here
%   params is struct with fields
%   - spikeChannels: same length as LFP channels, in the same order, if analyzing both
%   - cPtCal: conversion from spike sample indices to timestep of decimated LFP
%   returns:
%   - TODO add description
Output.VERBOSE('loading blackrock event file');
NEV = openNEV(spikeFilename,'read','nosave','noparse'); %note: add param 'report' for verbose
Output.VERBOSE('parsing serial IO packets');
pictureTriggers = NEV.Data.SerialDigitalIO;

%%%%% remove spike data from non-spike channels (e.g. reference electrodes)
spikesByChannel = repmat(struct('times',[],'units',[],'waveforms',[]),length(params.spikeChannels),1);
for i = 1:length(params.spikeChannels)
  %change units from sample index to ms; type from int32 to double
  tmp.times = params.cPtCal*double(NEV.Data.Spikes.Timestamps(NEV.Data.Spikes.Electrode == params.spikeChannels(i)));
  tmp.units = NEV.Data.Spikes.Unit(NEV.Data.Spikes.Electrode == params.spikeChannels(i));
  tmp.waveforms = NEV.Data.Spikes.Waveform(NEV.Data.Spikes.Electrode == params.spikeChannels(i));
  spikesByChannel(i) = tmp;
end
clear NEV
end

