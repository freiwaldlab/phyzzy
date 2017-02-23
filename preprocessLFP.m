function [ lfpData ] = preprocessLFP( lfpFilename, params )
%UNTITLED Summary of this function goes here
%   

% unpack params fields
needLFP = params.needLFP;  %todo: implement
needSpikes = params.needSpikes; %todo: needed?
lfpChannels = params.lfpChannels;
lfpChannelScaleBy = params.lfpChannelScaleBy; %converts raw values to microvolts
common_ref = params.common_ref; %not yet implemented; will allow software re-refrence across headstages
cPtCal = params.cPtCal; % conversion from spike sample indices to timestep of decimated LFP
decimateFactorPass1 = params.decimateFactorPass1; %note: product of the two decimate factors should be 30, if 1 khz samples desired
decimateFactorPass2 = params.decimateFactorPass2;
samPerMS = params.samPerMS; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)  
lfpFilter = params.filter; %if filtering desired, lfpFilter is a digitalFilter object

% load lfp data
tmp = openNSx(lfpFilename,'report','read');
lfpHeader = tmp.MetaTags;
assert(lfpHeader.SamplingFreq/(decimateFactorPass1*decimateFactorPass2) == 1000, 'error: expected lfp data to decimate to 1ks/sec');
lfpDataRaw = tmp.Data; %note: returns each channel as a row

% sort ns5 data so channel indexing matches indexing in ns5channels array
lfpChannelMap = zeros(length(lfpChannels),1);
for i = 1:length(lfpChannels)
  assert(any(lfpHeader.ChannelID == lfpChannels(i)), strcat('error: requested analysis for unrecorded LFP channel: ',num2str(lfpChannels(i))));
  lfpChannelMap(i) = find(lfpHeader.ChannelID == lfpChannels(i));
end
lfpChannelMap = lfpChannelMap(lfpChannelMap > 0);
lfpData = lfpDataRaw(lfpChannelMap,:);
clear lfpDataRaw
lfpDataDec = zeros(size(lfpData,1),ceil(size(lfpData,2)/(decimateFactorPass1*decimateFactorPass2)));
% convert scaled units to microvolts, and decimate (note: decimation broken
% into two steps, per matlab doc recommendation
disp('decimating, scaling, and filtering LFP');
for i = 1:size(lfpData,1)
  lfpDataDec(i,:) = lfpChannelScaleBy(i)*decimate(decimate(lfpChannelScaleBy(i)*lfpData(i,:),decimateFactorPass1),decimateFactorPass2);
  if isa(lfpFilter,'digitalFilter')
    lfpDataDec(i,:) = filter(lfpFilter, lfpDataDec(i,:));
  end
end
lfpData = lfpDataDec;
for i = 1:size(lfpData,1)
  lfpData(i,:) = lfpData(i,:) - mean(lfpData(i,:));
end
Output.VERBOSE('done decimating, scaling, and filtering LFP');
Output.DEBUG('size LFP data:'); Output.DEBUG(size(lfpData));
Output.DEBUG('size LFP data after decimation:'); Output.DEBUG(size(lfpData));
Output.DEBUG('LFP channel map'); Output.DEBUG(lfpChannelMap);
Output.DEBUG('channel id info'); Output.DEBUG(lfpHeader.ChannelID);
Output.DEBUG('size LFP data'); Output.DEBUG(size(lfpData));
end

