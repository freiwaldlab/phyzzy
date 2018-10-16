function [ db ] = trialDatabaseInit( dateSubj, runNum, numTrials, varargin )
%Initialize trial database struct
% dateSubj and runNum are strings
if isempty(varargin)
  db = struct();
else
  db = varargin{1};
end
db.(sprintf('sess%srun%s',dateSubj,runNum)).fields = struct();
db.(sprintf('sess%srun%s',dateSubj,runNum)).numTrials = numTrials;
end

