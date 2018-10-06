function [outputArg1,outputArg2] = testProcessRunsAsOne(runAnalysisInputsCmb,runAnalysisInputs)
%unit test suite for processRunsAsOne
%   
assert(length(runAnalysisInputsCmb.lfpData) == 2*length(runAnalysisInputs.lfpData),'failed test: LFP data length')
assert(length(runAnalysisInputsCmb.analogInData) == 2*length(runAnalysisInputs.analogInData),'failed test: analogIn data length')
assert(length(runAnalysisInputsCmb.taskData.taskEventStartTimes) == 2*length(runAnalysisInputs.taskData.taskEventStartTimes)+1,'failed test: taskData.taskEventStartTimes length')
assert(all(runAnalysisInputsCmb.trialIDsByEvent{1} == vertcat(runAnalysisInputs.trialIDsByEvent{1},runAnalysisInputs.trialIDsByEvent{1}+length(runAnalysisInputs.taskData.taskEventStartTimes))),'failed test: trialIDsByEvent');

end

