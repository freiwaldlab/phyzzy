function [ db, trialID ] = trialDatabaseAddTrial( db, dateSubj, runNum )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
runID = sprintf('%srun%s',dateSubj,runNum);
db.(runID).trialCounter = db.(runID).trialCounter + 1;
trialID = num2str(db.(runID).trialCounter);
db.(runID).trials.(trialID) = struct;
end

