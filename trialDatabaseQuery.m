function [ trialIDs ] = trialDatabaseQuery( fieldName, fieldValue, db, varargin )
%returns trialIDs of the trials matching a given name-value pair
%   - fieldName: a string
%   - fieldValue: string, number, numeric array, or function handle
%                 - for string and number
%                   - if field type is string or number, will check equality
%                   - if field type is cell or struct, will check for any equality
%                 - for other numeric array, will use fieldValue(1:end-1)
%                   as index, and check equality against fieldValue(end)
%                 - for function handle, handle is evaluated for each
%                   trial's cell, struct, or subarray
%   - db: a trial database struct
%   - varargin: 
%     - to specify dateSubj and runNum for possibly multi-run databases,
%       varargin should be e.g.'dateSub','171030ALAN','runNum','012',


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
assert(isfield(db.(runID).fieldValues,fieldName),'invalid field');
switch class(fieldValue)
  case 'char'
    if iscell(db.(runID).fieldValues.(fieldName))
      trialIDs = find(cellfun(@ (x) any(strcmp(x,fieldValue),db.(runID).fieldValues.(fieldName))));
    elseif isstruct(db.(runID).fieldValues.(fieldName))
      subfieldNames = fieldnames(db.(runID).fieldValues.(fieldName));
      subfieldName = subfieldNames{1};
      trialIDs = find(structfun(@ (x) any(strcmp(x,fieldValue),db.(runID).fieldValues.(fieldName).(subfieldName))));
    else
      trialIDs = find(strcmp(db.(runID).fieldValues.(fieldName),fieldValue));
    end
  case 'double'  
    if iscell(db.(runID).fieldValues.(fieldName))
      trialIDs = find(cellfun(@ (x) any(isequal(x,fieldValue),db.(runID).fieldValues.(fieldName))));
    elseif isstruct(db.(runID).fieldValues.(fieldName))
      subfieldNames = fieldnames(db.(runID).fieldValues.(fieldName));
      subfieldName = subfieldNames{1};
      trialIDs = find(structfun(@ (x) any(isequal(x,fieldValue),db.(runID).fieldValues.(fieldName).(subfieldName))));
    else
      if numel(fieldValue) == 1
        assert(size(db.(runID).fieldValues.(fieldName),2) == 1, 'error: invalid query shape for requested field shape');
        trialIDs = find(db.(runID).fieldValues.(fieldName) == fieldValue);
      else
        assert(numel(fieldValue) == size(db.(runID).fieldValues.(fieldName)), 'error: invalid query shape for requested field shape');
        trialIDs = find(squeeze(db.(runID).fieldValues.(fieldName)(:,fieldValue(1:end-1))) == fieldValue(end));
      end
    end
  case 'function_handle'
    if iscell(db.(runID).fieldValues.(fieldName))
      trialIDs = find(cellfun(@ (x) any(fieldValue(x),db.(runID).fieldValues.(fieldName))));
    elseif isstruct(db.(runID).fieldValues.(fieldName))
      subfieldNames = fieldnames(db.(runID).fieldValues.(fieldName));
      subfieldName = subfieldNames{1};
      trialIDs = find(structfun(@ (x) any(fieldValue(x),db.(runID).fieldValues.(fieldName).(subfieldName))));
    else
      trialIDs = find(squeeze(fieldValue(db.(runID).fieldValues.(fieldName))));
    end
  otherwise
    trialIDs = [];
end
end

