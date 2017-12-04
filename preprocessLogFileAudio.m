function [ taskData, stimTiming ] = preprocessLogFileAudio(logFilename,audioTrace,params )
% Load, decimate, and filter LFP, and index it by order in params. lfpChannels
%    decimation and filtering are optional; specified in params


% parse log file
disp('Loading visiko log file and converting to matlab structure');
assert(logical(exist(logFilename,'file')),'The stimulus log file you requested does not exist.');
logStruct = xml2struct(logFilename);
assert(isfield(logStruct.VISIKOLOG,'EndStimulation'),'Error: stimulation end not included in log file');
assert(strcmp(logStruct.VISIKOLOG.Attributes.tasktype,'Bitmap Continuous'),'Error: unknown Visiko task type. Must be Bitmap Continuous');


if isa(logStruct.VISIKOLOG.DOCDATA,'cell')
    s = input(sprintf('Visiko parameters changed during task; %d parameter sets found. Enter the number to analyze, or n to quit: ',length(logStruct.VISIKOLOG.DOCDATA)),'s');
    if strcmp(s,'n')
        return;
    end
    logStruct.VISIKOLOG.DOCDATA = logStruct.VISIKOLOG.DOCDATA{str2double(s)};
end

if isfield(logStruct.VISIKOLOG.Object{1,1},'Video')
    
    stimTiming.shortest = 1000*str2double(logStruct.VISIKOLOG.DOCDATA.BLOCKPARAMS_BCONT{1,1}.blockTime.Text);
    stimTiming.longest = 1000*str2double(logStruct.VISIKOLOG.DOCDATA.BLOCKPARAMS_BCONT{1,1}.blockTime.Text);
    stimTiming.ISI = 1000*str2double(logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT{1, 1}.blockStartPause  );
    stimParams = logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT;
    stimFilenames = {};
    % note: in case of videos, 'Objects' in the log file are either 'Video' fields
    stimFramesLost = -1*ones(length(logStruct.VISIKOLOG.Object),1);
    stimJumps = zeros(length(logStruct.VISIKOLOG.Object),2); %note: jumps can be negative, so we will use startTime for logical indexing
    stimStartTimesLog = -1*ones(length(logStruct.VISIKOLOG.Object),1); % Log affix signifies stimulation computer reference frame
    stimEndTimesLog = -1*ones(length(logStruct.VISIKOLOG.Object),1);
    fixationInTimesLog = -1*ones(length(logStruct.VISIKOLOG.Trigger),1);
    fixationOutTimesLog = -1*ones(length(logStruct.VISIKOLOG.Trigger),1);
    juiceOnTimesLog = -1*ones(length(logStruct.VISIKOLOG.Trigger),1);
    juiceOffTimesLog = -1*ones(length(logStruct.VISIKOLOG.Trigger),1);
    if isfield(logStruct.VISIKOLOG,'FixspotFlash')
        fixSpotFlashStartTimesLog = zeros(length(logStruct.VISIKOLOG.FixspotFlash),1);
        fixSpotFlashEndTimesLog = zeros(length(logStruct.VISIKOLOG.FixspotFlash),1);
    else
        fixSpotFlashStartTimesLog = 0;
        fixSpotFlashEndTimesLog = 0;
    end
    
    for i = 1:length(logStruct.VISIKOLOG.Object)
        stimulusStruct = logStruct.VISIKOLOG.Object{i};
        if isfield(stimulusStruct,'PictureCompleted')
            continue
        end
        tmp = regexp(stimulusStruct.Video.Filename.Text,'\','split');
        stimFilenames = vertcat(stimFilenames, strtrim(tmp{end})); %note: for some reason, filenames have trailing whitespace; trim it off
        stimFramesLost(i) = str2double(stimulusStruct.Frameslost.Text);
        stimStartTimesLog(i) = str2double(stimulusStruct.Start.Text);
        stimEndTimesLog(i) = str2double(stimulusStruct.End.Text);
    end
    Output.VERBOSE(sprintf('number of stimulus trials: %s',length(stimFilenames)));
    stimJumps = stimJumps(stimStartTimesLog >= 0,:); %note: jumps can be negative, so use startTime for logical indexing
    stimFramesLost = stimFramesLost(stimFramesLost >= 0);
    stimStartTimesLog = stimStartTimesLog(stimStartTimesLog >= 0);
    stimEndTimesLog = stimEndTimesLog(stimEndTimesLog >= 0);
    
    for i = 1:length(logStruct.VISIKOLOG.Trigger)
        triggerStruct = logStruct.VISIKOLOG.Trigger{i};
        switch triggerStruct.Name.Text
            case 'Fixation_in'
                fixationInTimesLog(i) = str2double(triggerStruct.Time.Text);
            case 'Fixation_out'
                fixationOutTimesLog(i) = str2double(triggerStruct.Time.Text);
            case 'RewardOn'
                juiceOnTimesLog(i)= str2double(triggerStruct.Time.Text);
            case 'RewardOff'
                juiceOffTimesLog(i)= str2double(triggerStruct.Time.Text);
            otherwise
                disp(strcat('unknown trigger name, ',TriggerStruct.Name.Text));
        end
    end
    fixationInTimesLog = fixationInTimesLog(fixationInTimesLog >= 0);
    fixationOutTimesLog = fixationOutTimesLog(fixationOutTimesLog >= 0);
    juiceOnTimesLog = juiceOnTimesLog(juiceOnTimesLog >= 0);
    juiceOffTimesLog = juiceOffTimesLog(juiceOffTimesLog >= 0);
    
    if isfield(logStruct.VISIKOLOG,'FixspotFlash')
        if iscell(logStruct.VISIKOLOG.FixspotFlash)
            for i = 1:length(logStruct.VISIKOLOG.FixspotFlash)
                flashStruct = logStruct.VISIKOLOG.FixspotFlash{i};
                fixSpotFlashStartTimesLog(i) = str2double(flashStruct.Time.Text);
                % note: flash duration is measured in sec in the log file; convert to msec
                fixSpotFlashEndTimesLog(i) = fixSpotFlashStartTimesLog(i) + 1000*str2double(flashStruct.Duration.Text);
            end
        else
            fixSpotFlashStartTimesLog(1) = str2double(logStruct.VISIKOLOG.FixspotFlash.Time.Text);
            fixSpotFlashEndTimesLog(1) = fixSpotFlashStartTimesLog(1) + 1000*str2double(logStruct.VISIKOLOG.FixspotFlash.Duration.Text);
        end
    end
end


disp('scaling audio');

for i = 1:size(audioTrace,1)
    audio(i,:) = audioTrace(i,:) - mean(audioTrace(i,:));
end

disp(['Session duration: ' num2str(length(audio)/30000/60) ' minutes'])
figure; hold on; plot(audio)
disp(['Considering first 2 seconds of BR signal for baseline'])
baseMean=abs(mean(audio(1:60000)));
baseStd=std(audio(1:60000));
localmax=abs(audio)>abs(baseMean+baseStd*10);
plot(localmax*500)


% first automatically detect 1st onset and add duration of video to find
% the rest of the onsets. Check in figure and if happy, continue. If not
% good, detect onsets from waveforms.
% onset = find(localmax(1,startAudioSearch:end),1);
duration_sec =  stimTiming.shortest/1000-0.4;
%
% for i=1:length(stimFilenames)
%     onseti = onset + duration_sec*30000*(i-1);
%     onsets(onseti)=1;
% end

startAudioSearch = uint64(1);
onsets = zeros(size(audio)); offsets = onsets;

% now, calculate visiko-to-blackrock conversion
for i=1:length(stimFilenames)
    if ~strcmp(stimFilenames{i}(end-6:end),'_SV.avi')
        remainingAudio = localmax(1,startAudioSearch:end);
        onset = find(remainingAudio,1);
        
        %[wave, fs] = audioread(fullfile(wavFiles,aFiles{i}));
        %duration_sec = length(wave)/fs;
        offset=round(onset+duration_sec*30000);
        
        onsets(1,int64(onset+startAudioSearch-1)) = 1;
        offsets(1,int64(offset+startAudioSearch-1)) = 1;
        
        startAudioSearch = uint64(offset+startAudioSearch);
    else
        remainingAudio = localmax(1,startAudioSearch:end);
        onset=startAudioSearch;
        onsets(1,int64(onset)) = 1;
        offset=startAudioSearch+(duration_sec*30000-1);
        offsets(1,int64(offset)) = 1;
        
        startAudioSearch = uint64(offset);
        
    end
end

figure;
plot(audio,'g'); hold on; plot(onsets*1000,'r');
plot(-offsets*1000,'b')

stimStartTimesBlk=find(onsets==1)/30000*1000; % in msec
stimEndTimesBlk=find(offsets==1)/30000*1000;


% now, calculate visiko-to-blackrock conversion
% first, if nev has one more start trigger than log, throw out final nev
% trigger (this is a known visiko bug, according to Michael Borisov's code)
assert(length(stimStartTimesBlk) - length(stimStartTimesLog) <= 1, 'Error: Start triggers missing from log file');
stimStartTimesBlk = stimStartTimesBlk(1:length(stimStartTimesLog));
%note: don't use first trigger in fit; sometimes off (known visiko bug)
logVsBlkModel = fitlm(stimStartTimesBlk(2:end), stimStartTimesLog(2:end));
disp(logVsBlkModel);
m = logVsBlkModel.Coefficients.Estimate(2);
y0 = logVsBlkModel.Coefficients.Estimate(1);
% for debugging
stimStartTimesFit = (1/m)*(stimStartTimesLog - y0);
disp(strcat('Max magnitude fit residual, msec: ',num2str(max(abs(stimStartTimesBlk-stimStartTimesFit)))));
% end for debugging
%stimEndTimesBlk = (1/m)*(stimEndTimesLog - y0);
fixationInTimesBlk = (1/m)*(fixationInTimesLog - y0);
fixationOutTimesBlk = (1/m)*(fixationOutTimesLog - y0);
juiceOnTimesBlk = (1/m)*(juiceOnTimesLog - y0);
juiceOffTimesBlk = (1/m)*(juiceOffTimesLog - y0);
fixSpotFlashStartTimesBlk = (1/m)*(fixSpotFlashStartTimesLog - y0);
fixSpotFlashEndTimesBlk = (1/m)*(fixSpotFlashEndTimesLog - y0);

% finally, build the output structure
taskData.stimFilenames = stimFilenames;
taskData.stimJumps = stimJumps;
taskData.stimFramesLost = stimFramesLost;
taskData.stimStartTimes = stimStartTimesBlk;
taskData.stimEndTimes = stimEndTimesBlk;
taskData.fixationInTimes = fixationInTimesBlk;
taskData.fixationOutTimes = fixationOutTimesBlk;
taskData.juiceOnTimes = juiceOnTimesBlk;
taskData.juiceOffTimes = juiceOffTimesBlk;
taskData.fixSpotFlashStartTimes = fixSpotFlashStartTimesBlk;
taskData.fixSpotFlashEndTimes = fixSpotFlashEndTimesBlk;
taskData.stimParams = stimParams;
taskData.RFmap = 0;
