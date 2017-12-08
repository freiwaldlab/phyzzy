function [  ] = ninePtEyeCalVisiko( outputFilename )
%ninePointCal determines x and y center offsets and gains from a 9-point
%calibration task
%   assumes the stimulation system is visiko
%   assumes the ephys system is Blackrock
assert(0, 'ninePtEyeCalVisiko is not yet implemented');

% note for serene: example eye data are in TEST_DATA/171205TEST, run 002
% also note: brief matlab session seemed to suggest that triggers are fix
% in, fix out, and some other stuff in the log file, but that only some of
% those get sent to Blackrock. Eg there were 325 each of fix in and fix
% out, and 486 blackrock SDIO packets. Note that the 9th and 10th bits of
% the SDIO packets were nonuniform across packets, but that there was no
% obvious e.g. toggle. Bit 9 had 351 1s and thus 135 zeros, while bit 10
% had 205 1s and thus 281 zeros. Plan is to try the tast with Alan, and see
% whether with real eye data it's easier to see what's going on in the
% logs. 
machine = 'laptop';

switch machine
  case 'rig'
    stimulusLogFile = '/Volumes/Users/FreiwaldLab/Desktop/';
    eyeDataFile = '/Volumes/Users-1/User/Desktop';
    triggerFile = '/Volumes/Users-1/User/Desktop';
  case 'laptop'
    eyeDataFile = '/Users/stephenserene/Desktop/Freiwald/TEST_DATA/Blackrock/171205TEST'; 
    stimulusLogFile = '/Users/stephenserene/Desktop/Freiwald/TEST_DATA/Visiko/171205TEST';
    triggerFile = '/Users/stephenserene/Desktop/Freiwald/TEST_DATA/Blackrock/171205TEST00';
end

eyeChX = 138;
eyeChY = 139;

gainX = 1;
gainY = 1;
offsetX = 1;
offsetY = 1;
flipX = 0;
flipY = 0;

save(outputFilename, 'gainX','gainY','offsetX','offsetY','flipX','flipY');
end

