function [ db ] = trialDatabaseInit( dateSubj, runNum )
%Initialize trial database struct
% dateSubj and runNum are strings
db.(sprintf('%srun%s',dateSubj,runNum)).trialCounter = 0;
db.(sprintf('%srun%s',dateSubj,runNum)).fieldValuesToTrialIDs = struct();
end

