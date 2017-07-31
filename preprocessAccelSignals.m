function [ analogInData ] = preprocessAccelSignals( analogInData,params ) 
%Applies volts to m/s^2 conversion for accelerometer signals.
%   - Optionally: applies coordinate system conversion based on calibration runs (not yet implemented)   
%   Note: supports multiple accelerometers
if ~params.needAccelCal
  return
end
assert(~any(strcmp(params.calMethods,'calFile')),'Requested file-based accelerometer calibration; not yet implemented');
  

for accel_i = 1:length(params.accelChannels)
  analogInData(params.accelChannels{accel_i}) = analogInData(params.accelChannels{accel_i}) - mean(analogInData(params.accelChannels{accel_i}),2);
  analogInData(params.accelChannels{accel_i}) = params.channelGains{accel_i} .* analogInData(params.accelChannels{accel_i});
end

end

