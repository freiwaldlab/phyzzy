function [ analogInData ] = preprocessAnalogIn( analogInFilename, params )
% Load, decimate, and filter analogIn data, and index by order in params.analogInChannels
%   Note: decimation and filtering optional; specified in params

% unpack params fields
if ~params.needAnalogIn
  analogInData = [];
  Output.VERBOSE('not loading analog in data');
  return
end
analogInChannels = params.analogInChannels;
analogInChannelScaleBy = params.analogInChannelScaleBy; %converts raw values to microvolts
decimateFactorPass1 = params.decimateFactorPass1; %note: product of the two decimate factors should be (1/1000)*sam/sec, if 1 khz samples desired
decimateFactorPass2 = params.decimateFactorPass2;
samPerMS = params.samPerMS; %THIS IS AFTER DECIMATION,  
analogInFilters = params.filters; %if filtering desired, analogInFilter is a digitalFilter object

% load analogIn data
assert(logical(exist(analogInFilename,'file')),'The analog input file you requested does not exist.');
tmp = openNSx(analogInFilename,'report','read');
analogInHeader = tmp.MetaTags;
assert(analogInHeader.SamplingFreq/(decimateFactorPass1*decimateFactorPass2) == 1000, 'error: expected analogIn data to decimate to 1ks/sec');
analogInDataRaw = tmp.Data; %note: returns each channel as a row

% sort data so channel indexing matches indexing in analogInChannels array
analogInChannelMap = zeros(length(analogInChannels),1);
for i = 1:length(analogInChannels)
  assert(any(analogInHeader.ChannelID == analogInChannels(i)), strcat('error: requested analysis for unrecorded analogIn channel: ',num2str(analogInChannels(i))));
  analogInChannelMap(i) = find(analogInHeader.ChannelID == analogInChannels(i));
end
analogInChannelMap = analogInChannelMap(analogInChannelMap > 0);
analogInData = analogInDataRaw(analogInChannelMap,:);
clear analogInDataRaw

% filter the analogIn data using 2 step process.
if isfield(params, 'filterPad')
  filterPad = params.filterPad;
else
  filterPad = 0;
end
analogInDataDecPadded = zeros(size(analogInData,1),ceil(size(analogInData,2)/(decimateFactorPass1*decimateFactorPass2))+2*filterPad);
% convert scaled units to microvolts, and decimate (note: decimation broken
% into two steps, per matlab doc recommendation
disp('decimating, scaling, and filtering analog inputs');
for i = 1:size(analogInData,1)
  analogInDataDecPadded(i,filterPad+1:end-filterPad) = analogInChannelScaleBy(i)*decimate(decimate(analogInData(i,:),decimateFactorPass1),decimateFactorPass2);
  %analogInDataDecPadded(i,1:filterPad) = analogInDataDecPadded(i,filterPad+1)*analogInDataDecPadded(i,1:filterPad);
  %analogInDataDecPadded(i,end-(filterPad-1):end) = analogInDataDecPadded(i,end-filterPad)*analogInDataDecPadded(i,end-(filterPad-1):end);
  analogInFilter = analogInFilters{i};
  if params.plotFilterResult && (isa(analogInFilter,'digitalFilter') || length(analogInFilter) == 2)
    figure();
    plot(analogInDataDecPadded(1,filterPad+100000:filterPad+105000),'color','r');
    hold on
  end
  if isa(analogInFilter,'digitalFilter')
    disp('using digital filter object');
    analogInDataDecPadded(i,:) = filtfilt(analogInFilter, analogInDataDecPadded(i,:));
  else
    if length(analogInFilter) == 2
      analogInDataDecPadded(i,:) = filtfilt(analogInFilter(1),analogInFilter(2),analogInDataDecPadded(i,:));
    end
    if params.plotFilterResult && (isa(analogInFilter,'digitalFilter') || length(analogInFilter) == 2)
      plot(analogInData(i,filterPad+100000:filterPad+105000),'color','b');
      legend({'raw','filtered'});
      drawnow;
    end
  end
end
analogInData = analogInDataDecPadded(:,filterPad+1:end-filterPad);

Output.VERBOSE('done decimating, scaling, and filtering analogIn');
Output.DEBUG('size analogIn data:'); Output.DEBUG(size(analogInData));
Output.DEBUG('size analogIn data after decimation:'); Output.DEBUG(size(analogInData));
Output.DEBUG('analogIn channel map'); Output.DEBUG(analogInChannelMap);
Output.DEBUG('channel id info'); Output.DEBUG(analogInHeader.ChannelID);
Output.DEBUG('size analogIn data'); Output.DEBUG(size(analogInData));
end

