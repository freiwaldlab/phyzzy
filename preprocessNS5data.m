function [ lfpData, diodeTrace, audioTrace ] = preprocessNS5( lfpFilename, ephysParams, photodiodeParams, audioParams )
% Load, decimate, and filter LFP, and index it by order in params. lfpChannels
%    decimation and filtering are optional; specified in params

% unpack params fields
if ~ephysParams.needLFP
  lfpData = [];
  Output.VERBOSE(sprintf('not loading LFP data from %s',lfpFilename));
  return
end
lfpChannels = ephysParams.lfpChannels;
lfpChannelScaleBy = ephysParams.lfpChannelScaleBy; %converts raw values to microvolts
common_ref = ephysParams.commonRef; %not yet implemented; will allow software re-refrence across headstages
cPtCal = ephysParams.cPtCal; % conversion from spike sample indices to timestep of decimated LFP
decimateFactorPass1 = ephysParams.decimateFactorPass1; %note: product of the two decimate factors should be 30, if 1 khz samples desired
decimateFactorPass2 = ephysParams.decimateFactorPass2;
samPerMS = ephysParams.samPerMS; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)  
lfpFilter = ephysParams.filter; %if filtering desired, lfpFilter is a digitalFilter object

% load lfp data
assert(logical(exist(lfpFilename,'file')),'The LFP file you requested does not exist.');
tmp = openNSx(lfpFilename,'report','read');
lfpHeader = tmp.MetaTags;
assert(lfpHeader.SamplingFreq/(decimateFactorPass1*decimateFactorPass2) == 1000, 'error: expected lfp data to decimate to 1ks/sec');
lfpDataRaw = tmp.Data; %note: returns each channel as a row

% sort lfp data so channel indexing matches indexing in lfpChannels array
lfpChannelMap = zeros(length(lfpChannels),1);
for i = 1:length(lfpChannels)
  assert(any(lfpHeader.ChannelID == lfpChannels(i)), strcat('error: requested analysis for unrecorded LFP channel: ',num2str(lfpChannels(i))));
  lfpChannelMap(i) = find(lfpHeader.ChannelID == lfpChannels(i));
end
lfpChannelMap = lfpChannelMap(lfpChannelMap > 0);
lfpData = lfpDataRaw(lfpChannelMap,:);
%clear lfpDataRaw
filterPad = 50000;
lfpDataDecPadded = zeros(size(lfpData,1),ceil(size(lfpData,2)/(decimateFactorPass1*decimateFactorPass2))+2*filterPad);
% convert scaled units to microvolts, and decimate (note: decimation broken
% into two steps, per matlab doc recommendation
disp('decimating, scaling, and filtering LFP');
for i = 1:size(lfpData,1)
  lfpData(i,:) = lfpData(i,:) - mean(lfpData(i,:));
end
for i = 1:size(lfpData,1)
  lfpDataDecPadded(i,filterPad+1:end-filterPad) = lfpChannelScaleBy(i)*decimate(decimate(lfpData(i,:),decimateFactorPass1),decimateFactorPass2);
  lfpDataDecPadded(i,1:filterPad) = lfpDataDecPadded(i,filterPad+1)*lfpDataDecPadded(i,1:filterPad);
  lfpDataDecPadded(i,end-(filterPad-1):end) = lfpDataDecPadded(i,end-filterPad)*lfpDataDecPadded(i,end-(filterPad-1):end);
  if i == 1 && ephysParams.plotFilterResult
    figure();
    plot(lfpDataDecPadded(1,filterPad+100000:filterPad+105000),'color','r');
    hold on
  end
  if isa(lfpFilter,'digitalFilter')
    disp('using digital filter object');
    lfpDataDecPadded(i,:) = filtfilt(lfpFilter, lfpDataDecPadded(i,:));
  else
    if length(lfpFilter) == 2
      lfpDataDecPadded(i,:) = filtfilt(lfpFilter(1),lfpFilter(2),lfpDataDecPadded(i,:));
    end
  end
end
lfpData = lfpDataDecPadded(:,filterPad+1:end-filterPad);
if ephysParams.plotFilterResult
  plot(lfpData(1,100000:105000),'color','b');
  legend({'raw','filtered'});
  drawnow;
end
disp(mean(lfpData(1,:)));

Output.VERBOSE('done decimating, scaling, and filtering LFP');
Output.DEBUG('size LFP data:'); Output.DEBUG(size(lfpData));
Output.DEBUG('size LFP data after decimation:'); Output.DEBUG(size(lfpData));
Output.DEBUG('LFP channel map'); Output.DEBUG(lfpChannelMap);
Output.DEBUG('channel id info'); Output.DEBUG(lfpHeader.ChannelID);
Output.DEBUG('size LFP data'); Output.DEBUG(size(lfpData));

if ~photodiodeParams.needPhotodiode
  photodiodeFrameTimes = [];
  
else
    diodeChannel = find(lfpHeader.ChannelID == photodiodeParams.channels);
    diodeTrace = lfpDataRaw(diodeChannel,:);

end

if ~audioParams.needAudio
  audio = [];
else
    audioChannel = find(lfpHeader.ChannelID == audioParams.channels);
    audioTrace = lfpDataRaw(audioChannel,:);

end

end