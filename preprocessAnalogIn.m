function [ analogInData ] = preprocessAnalogIn( analogInFilename, params )
%preprocessAnalogIn is not yet implemented! TODO
%   
analogInData = [];
if ~params.needAccel && ~params.needEyes
  return;
end
tmp = openNSx(analogInFilename,'report','read');
header = tmp.MetaTags;
data = ns2obj.Data;
end

