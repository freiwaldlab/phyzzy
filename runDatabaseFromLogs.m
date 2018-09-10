function [  ] = runDatabaseFromLogs( dateSubjList, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

stimulusLogVolume = '/Volumes/Users/FreiwaldLab/Desktop';
if ~isempty(varargin)
  db = varargin{1};
else
  db = struct();
end
for dateSubj_i = 1:length(dateSubjList)
  dateSubject = dateSubjList{dateSubj_i};
  noRunsFor = 0;
  runNum = 0;
  while noRunsFor < 3
    runNum = runNum + 1;
    if runNum < 10
      runNumChar = sprintf('00%d',runNum);
    elseif runNum < 100
      runNumChar = sprintf('0%d',runNum);
    elseif runNum < 1000
      runNumChar = sprintf('%d',runNum);
    else
      break
    end
    logFilename = sprintf('%s/%s/%s0%s.log',stimulusLogVolume,dateSubject,dateSubject,runNumChar);
    try
      logStruct = xml2struct(logFilename);
    catch
      noRunsFor = noRunsFor + 1;
      continue
    end
    stimulusParams = logStruct.VISIKOLOG.DOCDATA.OBJECTPARAMS_BCONT;
    taskParams = logStruct.VISIKOLOG.DOCDATA.BCONTPARAMS;
    runParams = logStruct.VISIKOLOG.SESSIONDATA;
    fixationParams = logStruct.VISIKOLOG.DOCDATA.FIXSPOTPARAMS_MARIA;
    
    db.dateSubj.(strcat('run',runNumChar)).stimList = stimulusParams.imageListName.Text;
    db.dateSubj.(strcat('run',runNumChar)).stimPosRho = str2num(stimulusParams.stimPosRho.Text);
    db.dateSubj.(strcat('run',runNumChar)).stimPosPhi = str2num(stimulusParams.stimPosPhi.Text);
    db.dateSubj.(strcat('run',runNumChar)).stimFixedSize = str2num(stimulusParams.stimFixedSize.Text);
    if strcmp(stimulusParams.activateJumping.Text,'1')
      db.dateSubj.(strcat('run',runNumChar)).paradigm = 'RF map';
      db.dateSubj.(strcat('run',runNumChar)).rfMapWidth = str2num(stimulusParams.jumpExtentHorizontal.Text);
      db.dateSubj.(strcat('run',runNumChar)).rfMapHeight = str2num(stimulusParams.jumpExtentVertical.Text);
      db.dateSubj.(strcat('run',runNumChar)).rfMapGridSpacing = str2num(stimulusParams.jumpGridSize.Text);
    end  
    if str2num(stimulusParams.pictureTimeFrom.Text) == str2num(stimulusParams.pictureTimeTo.Text)
      db.dateSubj.(strcat('run',runNumChar)).stimDuration = str2num(stimulusParams.pictureTimeTo.Text);
    else
      db.dateSubj.(strcat('run',runNumChar)).stimDurationMin = str2num(stimulusParams.pictureTimeFrom.Text);
      db.dateSubj.(strcat('run',runNumChar)).stimDurationMax = str2num(stimulusParams.pictureTimeTo.Text);
    end
    db.dateSubj.(strcat('run',runNumChar)).ISI = str2num(stimulusParams.pauseTime.Text);
    db.dateSubj.(strcat('run',runNumChar)).isiJitter = max(str2num(stimulusParams.randomPauseVariation.Text),str2num(stimulusParams.maxCumulativeVariation.Text));
    
    %db.dateSubj.(strcat('run',runNumChar)).backgroundLevel = str2num(taskParams.backgroundColor.Text);
    db.dateSubj.(strcat('run',runNumChar)).rewardPeriod = str2num(taskParams.rewardPeriod.Text);
    db.dateSubj.(strcat('run',runNumChar)).frameRate = str2num(runParams.Framerate.Text);
    
    db.dateSubj.(strcat('run',runNumChar)).fixWinWidth = str2num(fixationParams.ewWidth.Text);
    db.dateSubj.(strcat('run',runNumChar)).fixWinHeight = str2num(fixationParams.ewHeight.Text);
    db.dateSubj.(strcat('run',runNumChar)).fixSpotSize = str2num(fixationParams.sizePoint.Text);
    
    db.dateSubj.(strcat('run',runNumChar)).numBlocks = length(logStruct.VISIKOLOG.BLOCKSTART);
    
    db.dateSubj.(strcat('run',runNumChar)).numBlocks = length(logStruct.VISIKOLOG.BLOCKSTART);
    db.dateSubj.(strcat('run',runNumChar)).numBlocks = length(logStruct.VISIKOLOG.Object);
    
    if ~strcmp(stimulusParams.distortionType.Text,'0')
      %todo: bubbles handling
    end
    if (~strcmp(stimulusParams.jumpSimBitmaps.Text,'1')) && strcmp(stimulusParams.activateJumping.Text,'1')
    end
    if ~strcmp(stimulusParams.stimScalingType.Text,'1')
    else
    end
    if ~strcmp(stimulusParams.stimAlpha.Text,'100')
    end
     
  end
end
end














