function [ output_args ] = makePlots(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load(analysisParamFilename);
% I should be able to eliminate most of this
pictureLabelsTmp = pictureLabels; %note: hack to avoid overwriting list of not presented stimuli
load(picParamsFilename);
pictureLabels = pictureLabelsTmp; % conclusion of hack
channelNames = ephysParams.channelNames;
spikeChannels = ephysParams.spikeChannels;
lfpChannels = ephysParams.lfpChannels;
psthPre = psthParams.psthPre;
psthPost = psthParams.psthPost;
smoothingWidth = psthParams.smoothingWidth;

lfpPreAlign = lfpAlignParams.msPreAlign; 
lfpPostAlign = lfpAlignParams.msPostAlign;

lfpPaddedBy = tfParams.movingWin(1)/2;

movingWin = tfParams.movingWin;
specgramRowAve = tfParams.specgramRowAve;
samPerMS = ephysParams.samPerMS;
if frCalcOff < frCalcOn
  frCalcOff = psthImDur+frCalcOn;
end







end

