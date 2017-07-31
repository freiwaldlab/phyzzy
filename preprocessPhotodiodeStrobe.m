function [ analogInData ] = preprocessPhotodiodeStrobe( analogInData, params )
%preprocessPhotodiodeStrobe cleans a photodiode signal
%   Detailed explanation goes here
if ~params.needPhotodiode
  return
end
error('photodiode strobe preprocessing not yet implemented');

end

