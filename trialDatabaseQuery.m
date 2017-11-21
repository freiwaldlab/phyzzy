function [ trialIDs ] = trialDatabaseQuery( fieldName, fieldValue, db, varargin )
%returns trialIDs of the trials matching a given name-value pair
%   - fieldName: a string
%   - fieldValue: a string, a double, or a 2 element array of doubles
%   - db: a trial database struct
%   - varargin: 
%     - to specify dateSubj and runNum for possibly multi-run databases,
%       varargin should be e.g.'dateSub','171030ALAN','runNum','012',
%     - to automatically plot the trial summary of the identified trial(s),
%       include 'plotTrialSummary' in varargin
%     - to supress command prompt output, include 'silent' in varargin
%
%   todo: implement range lookup for numerical values

if length(varargin) > 4 && strcmp(varargin{1},'dateSubj') && strcmp(varargin{3},'runNum')
  dateSubj = varargin{1};
  runNum = varargin{4};
  runID = (sprintf('%srun%s',dateSubj,runNum));
else
  fields = fieldnames(db);
  if length(fields) > 1
    disp(fields);
    dateSubj = input('The trial database given contains multiple runs (see above). Please specify a run as date and subject (typically e.g. 171030ALAN)');
    runNum = input('Now specify a run number, typically a three-character string like 012');
    runID = (sprintf('%srun%s',dateSubj,runNum));
  else
    runID = fields{1};
  end
end
assert(isfield(db,runID),'database given does not contain data from the run requested');
assert(isfield(db.(runID).fieldValuesToTrialIDs,fieldName),'invalid field');
switch class(fieldValue)
  case 'char'
    fieldValueStr = fieldValue;
  case 'double'  %important: this line sets the query precision to 4 digits
    if numel(fieldValue) == 1
      fieldValueStr = num2str(fieldValue,4); %important: this line sets the query precision to 4 digits
    else
      if numel == 2
        fieldValueStr = sprintf('%s_%s',num2str(fieldValue(1),4),num2str(fieldValue(2),4));
      else
        error('invalid field value type');
      end
    end
end
if iskey(db.(runID).fieldValuesToTrialIDs.(fieldName),fieldValueStr)
  trialIDs = db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr);
  if length(trialIDs) == 1 && ~any(strcmp('silent',varargin))
    sprintf('found a trial with the requested field-value pair: ID is %s\n',trialIDs{1});
  end
else  
  disp('No trial found with the requested field-value pair');
end

if any(strcmp('plotTrialSummary',varargin))
  for trial_i = 1:length(trialIDs)
    plotTrialSummary(db,trialIDs{trial_i},runID)
  end
end

if any(strcmp('plotGroupSummary',varargin))
  % this will plot a summary of all trials that share a given property wtih
  % the querried trial
end

if any(strcmp('highlightTrialOnFigure',varargin))
end

end

