function spikeDataBankPath = processBatchAnalysis( varargin )
% Function uses batchAnalysisParamFile to coordinate the compilation of the
% desired 'spikeDataBank'. 
% Input
% - string of name of buildBatchAnalysisParamFile function
% - usePreprocessedSpikeData: a 1 or 0 indicating whether to regenerate the
% spikeDataBank (0), or use what is found (1). Default = 1.
% Output
% - paths to relevant .mat files (sliced up spikeDataBank + other
% variables). spikeDataBank files end with '_N', while the other variables
% end with 'Vars'. Vars file is 1st in cell array of file paths. 

%% Set up
addpath(genpath('dependencies'));
rmpath(genpath('dependencies/mvgc_v1.0')); %note: Not sure if this is appropriate replacement for genpath_exclude. previous line caused issues in parallel runs.
addpath('buildAnalysisParamFileLib');

switch length(varargin)
  case 1
    batchAnalysisParamFile  = varargin{1};
    usePreprocessedSpikeData = 1;
  case 2
    batchAnalysisParamFile  = varargin{1};
    usePreprocessedSpikeData  = varargin{2};
  otherwise
    error('incorrect number of arguments')
end

addpath('buildAnalysisParamFileLib');
batchAnalysisParamFilename = feval(batchAnalysisParamFile);
load(batchAnalysisParamFilename);

%% Get the relevant data
spikeDataBaseFile = [outputDir '/' preprocessParams.spikeDataFileName];
files = dir([spikeDataBaseFile '*.mat']);

if ~isempty(files)
  disp('spikeDataBank found, Returning paths')
  spikeDataBankPath = cell(length(files),1);
  for file_ind = 1:length(files)
    spikeDataBankPath{file_ind} = fullfile(files(file_ind).folder,files(file_ind).name);
  end
else
  disp('spikeDataBank Not found, generating now...')
  %Cycle through analyzedData.mat files, store and organize the relevant structures.
  preprocessedList = dir([analysisDirectory filesep '**' filesep 'preprocessedData.mat']);
  analyzedList = dir([analysisDirectory filesep '**' filesep 'analyzedData.mat']);
  assert(length(preprocessedList) == length(analyzedList), 'Lists arent same length, confirm every preprocessed file has an analyzed file')
  
  allStimuliVec = {};
  sessionList = cell(length(preprocessedList),1);
  spikeDataBank = struct();
  
  for ii = 1:length(preprocessedList)
    tmp = load([preprocessedList(ii).folder filesep preprocessedList(ii).name],'spikesByEvent','eventIDs','eventCategories','preAlign','postAlign');
    tmp2 = load([analyzedList(ii).folder filesep analyzedList(ii).name], 'analysisParamFilename','dateSubject', 'runNum', 'groupLabelsByImage','psthByImage','attendedObjData');
    tmp3 = load([analyzedList(ii).folder filesep 'AnalysisParams.mat'], 'psthParams');
    
    sessField = sprintf('S%s%s', tmp2.dateSubject, tmp2.runNum);
    allStimuliVec = [allStimuliVec; tmp.eventIDs];
    sessionList{ii} = [tmp2.dateSubject tmp2.runNum];
    spikeDataBank.(sessField).dateSubject = tmp2.dateSubject;
    spikeDataBank.(sessField).runNum = tmp2.runNum;
    spikeDataBank.(sessField).spikesByEvent = tmp.spikesByEvent;
    spikeDataBank.(sessField).psthByImage = tmp2.psthByImage;
    spikeDataBank.(sessField).attendedObjData = tmp2.attendedObjData;
    spikeDataBank.(sessField).eventIDs = tmp.eventIDs;
    spikeDataBank.(sessField).eventCategories = tmp.eventCategories;
    spikeDataBank.(sessField).groupLabelsByImage = tmp2.groupLabelsByImage;
    spikeDataBank.(sessField).start = -tmp3.psthParams.psthPre;
    spikeDataBank.(sessField).stimDur = tmp3.psthParams.psthImDur;
    spikeDataBank.(sessField).end = tmp3.psthParams.psthImDur + tmp3.psthParams.psthPost;
    %spikeDataBank.(sessField).groupLabel = target; %Now in
    %slidingANOVAparams.
    spikeDataBank.(sessField).figDir = preprocessedList(ii).folder;
  end
  allStimuliVec = unique(allStimuliVec);
  spikeDataBankPath = saveSpikeDataBank(spikeDataBank, 5, 'save',outputDir);
  clearvars spikeDataBank
  nonSpikeSaveName = fullfile(outputDir, [preprocessParams.spikeDataFileName 'Vars']);
  spikeDataBankPath = [nonSpikeSaveName; spikeDataBankPath];
  save(nonSpikeSaveName)
end

%Return relevant paths to batchAnalysis.
runBatchAnalysis(spikeDataBankPath)

end
