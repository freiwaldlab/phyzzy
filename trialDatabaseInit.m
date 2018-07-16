function [ db ] = trialDatabaseInit( dateSubj, runNum, numTrials )
%Initialize trial database struct
% dateSubj and runNum are strings
db.(sprintf('%srun%s',dateSubj,runNum)).fieldValues = struct();
db.(sprintf('%srun%s',dateSubj,runNum)).numTrials = numTrials;
end

