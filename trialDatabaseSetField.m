function [ db ] = trialDatabaseSetField( fieldName, fieldValue, db, trialID, varargin )
%trialDatabaseSetField adds a field-value pair to the trial database
%   Places the key value pair into a dictionary for O(1) look-up
%   - db: the trial database struct; initialize wtih trialDatabaseInit
%   - fieldValue: acceptable types are string, double, 1x2 or 2x1 array of
%     doubles, or Nx2 array of doubles. 
%   - fieldName: a string
%   - dateSubj: a string
%   - runNum: a string
%   - trialID: a string, returned by trialDatabaseAddTrial
%   
%   Todo: add sorted list option for double value types, to allow range lookup

if ~isempty(varargin)
  dateSubj = varargin{1};
  runNum = varargin{2};
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

if ~isfield(db.(runID).fieldValuesToTrialIDs,fieldName)
  db.(runID).fieldValuesToTrialIDs.(fieldName) = containers.Map();
end
switch class(fieldValue)
  case 'char'
    fieldValueStr = fieldValue;
    if iskey(db.(runID).fieldValuesToTrialIDs.(fieldName),fieldValueStr)
      db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr) = [db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr);trialID];
    else  
      db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr) = {trialID};
    end
  case 'double'  %important: this line sets the query precision to 4 digits
    if numel(fieldValue) == 1
      fieldValueStr = num2str(fieldValue,4); %important: this line sets the query precision to 4 digits
      if iskey(db.(runID).fieldValuesToTrialIDs.(fieldName),fieldValueStr)
        db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr) = [db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr);trialID];
      else
        db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr) = {trialID};
      end
    else
      if numel == 2
        fieldValueStr = sprintf('%s_%s',num2str(fieldValue(1),4),num2str(fieldValue(2),4));
        if iskey(db.(runID).fieldValuesToTrialIDs.(fieldName),fieldValueStr)
          db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr) = [db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr);trialID];
        else
          db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr) = {trialID};
        end
      else
        if size(fieldValue,2) == 2
          for pair_i = 1:length(fieldValue)
            fieldValueStr = sprintf('%s_%s',num2str(fieldValue(1),4),num2str(fieldValue(2),4));
            if iskey(db.(runID).fieldValuesToTrialIDs.(fieldName), fieldValueStr)
              db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr) = [db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr);trialID];
            else
              db.(runID).fieldValuesToTrialIDs.(fieldName)(fieldValueStr) = {trialID};
            end
          end
        else
          error('Database insertion not implemented for the data field shape given');
        end
      end
    end
  otherwise
    error('Database insertion not implemented for values of the type given');
end
db.(runID).trials.(trialID).(fieldName) = fieldValue;
end

