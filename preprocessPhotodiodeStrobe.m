function [ diodeTriggers ] = preprocessPhotodiodeStrobe( inputData, params )
%preprocessPhotodiodeStrobe gives backward compatibility for the old name
%and run switch syntax of preprocessStrobe
if params.needPhotodiode
  diodeTriggers = preprocessStrobe(inputData, params);
else
  diodeTriggers = [];
end
end

