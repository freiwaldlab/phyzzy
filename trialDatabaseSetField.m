function [ db ] = trialDatabaseSetField( fieldName, fieldValue, db, trialID, varargin )
%trialDatabaseSetField adds a field-value pair to the trial database
%   Places the key value pair into a trial database
%   - db: the trial database struct; initialize wtih trialDatabaseInit
%   - fieldValue: acceptable types are string, numeric array, cell array
%   - fieldName: a string; must be valid struct field name
%   - dateSubj: a string
%   - runNum: a string
%   - trialID: a positive integer
%   


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

switch class(fieldValue)
  case 'char'
    if ~isfield(db.(runID).fieldValues,fieldName)
      db.(runID).fieldValues.(fieldName) = cell(numTrials,1);
    end
    db.(runID).fieldValues.(fieldName){trialID} = fieldValue;
  case 'double'
    if ~isfield(db.(runID).fieldValues,fieldName)
      if numel(fieldValue) == 1
        db.(runID).fieldValues.(fieldName) = NaN(numTrials,1);
      elseif length(size(fieldValue)) == 2 && min(size(fieldValue)) == 1
        db.(runID).fieldValues.(fieldName) = NaN(numTrials,length(fieldValue));
      else
        db.(runID).fieldValues.(fieldName) = NaN([numTrials,size(fieldValue)]);
      end
    end
    
    if numel(fieldValue) == 1
      db.(runID).fieldValues.(fieldName)(trialID) = fieldValue;
    elseif length(size(fieldValue)) == 2 && min(size(fieldValue)) == 1
      db.(runID).fieldValues.(fieldName)(trialID,:) = fieldValue;
    else
      db.(runID).fieldValues.(fieldName)(trialID,:) = fieldValue(:);
    end
  case 'cell'
    if ~isfield(db.(runID).fieldValues,fieldName)
      db.(runID).fieldValues.(fieldName) = cell(numTrials,1);
    end
    db.(runID).fieldValues.(fieldName){trialID} = fieldValue;
  case 'struct'
    assert(length(fieldnames(fieldValue)) == 1, 'Error: database insertion not implemented for multi-field structs');
    if ~isfield(db.(runID).fieldValues,fieldName)
      db.(runID).fieldValues.(fieldName) = repmat(struct(fieldnames(fieldValue),[]),numTrials,1);
    end
    db.(runID).fieldValues.(fieldName){trialID} = fieldValue;
  otherwise
    error('Database insertion not implemented for values of the type given');
end
end

