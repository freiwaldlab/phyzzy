% Merge runs in which exactly the same task was run, with same timing
% conditions and same preprocessing was done, with same units defined

clear all; close all
inputDir='/Volumes/FreiwaldShares/slandi/ephys/tests_2017/Analyzed/171115Buster/Basic/';
runs = {'036';'037';'038'}; 
name='';
for i=1:numel(runs)
    dataRun{i}=load([inputDir runs{i} '/preprocessedData.mat']);    
    name=[name runs{i}];
end

outputDir=[inputDir name '/'];
if ~exist(outputDir)
    mkdir(outputDir)
end

% take first run info for all constants between merged runs
mergeRuns=dataRun{1};
load(mergeRuns.analysisParamFilename);
runNum=name;
mergeRuns.analysisParamFilename=[outputDir '/analysisParams.mat'];
analysisParamFilename=[outputDir '/analysisParams.mat'];
analysisParamFilenameParts{1}=[outputDir '/analysisParams'];
analysisParamFilenameParts{2}='mat';
analysisParamFilenameStem=analysisParamFilenameParts{1};
outDir=outputDir;
save(mergeRuns.analysisParamFilename)

for irun=2:numel(runs)
    
    nImg = numel(dataRun{irun}.spikesByImage);
    
    for im=1:nImg
        for iCh=1:numel(mergeRuns.channelUnitNames{1})
            mergeRuns.spikesByImage{im}{1}{iCh} = [mergeRuns.spikesByImage{im}{1}{iCh}; dataRun{irun}.spikesByImage{im}{1}{iCh}];            
        end
        
        mergeRuns.lfpByImage{im} = cat(3,mergeRuns.lfpByImage{im}, dataRun{irun}.lfpByImage{im});  
        
    end
    
    for icat=1:numel(mergeRuns.categoryList)
        
        for iCh=1:numel(mergeRuns.channelUnitNames{1})
            mergeRuns.spikesByCategory{icat}{1}{iCh} = [mergeRuns.spikesByCategory{icat}{1}{iCh}; dataRun{irun}.spikesByCategory{icat}{1}{iCh}];            
        end
        
        mergeRuns.lfpByCategory{icat} = cat(3,mergeRuns.lfpByCategory{icat}, dataRun{irun}.lfpByCategory{icat});  
        
    end
    
end
                
save([outDir 'merge.mat'])
%%
runAnalyses(mergeRuns)