function [] = processRunBatch(runList, analysisParamFile)
%Processes a batch of runs, using the parameters in analysisParamFile
%   runList syntax is {{'180314MO; {'002';'003'}};{'180315MO; {'002';'003';'004'}}}
%
%   to do: add support for parallel execution (need to maintain graceful
%   compaitbility with matlab licenses that lack parallel computing toolbox
load(analysisParamFile);
for dateSubj_i = 1:lenght(runList)
  dateSubj = runList{dateSubj_i}{1};
  for run_i = 1:length(runList{dateSubj_i}{2})
  runNum = runList{dateSubj_i}{2}{run_i};
  analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);   %#ok
  lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);        
  spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
  taskFilename = sprintf('%s/%s/%s%s.mat',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
  outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
  analysisParamFilename = strcat(outDir,analysisParamFilenameStem);
  preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok
  if ~exist(outDir,'dir')
    mkdir(outDir);
  end
  save(analysisParamFilename); 
  processRun('preprocessed', analysisParamFilename);
end
end

