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
%   - varargin: name-value pairs
%     - dateSubj: a string, e.g. '180426ALAN'
%     - runNum: a string, e.g. '002'
%     - rawIndices: logical array, Nx1, N is number of indices into
%                   fieldName.data; if 1, treat indices in query as
%                   indices, even if index-value conversion given in metadata
%     - (todo) closest: int, if defined and > 0, return the closest N matches 


assert(mod(length(varargin),2) == 0, 'trialDatabaseQuery takes optional arguments as name-value pairs; odd number provided');
for argPair_i=1:length(varargin)/2
  argName = varargin{1+2*(argPair_i-1)};
  argVal = varargin{2+2*(argPair_i-1)};
  if strcmp(argName,'dateSubj')
    dateSubj = argval;
  elseif strcmp(argName,'runNum')
    runNum = argVal;
  elseif strcmp(argName,'rawIndices')
    assert(length(argVal) == length(fieldValue)-1,'rawIndices argument must have length length(fieldValue)-1');
    useRawIndices = argVal;
  end
end
if exist('dateSubj','var') && exist('runNum','var')
  runID = (sprintf('sess%srun%s',dateSubj,runNum));
else
  dateSubjList = fieldnames(db);
  if length(dateSubjList) > 1
    disp(dateSubjList);
    dateSubj = input('The trial database given contains multiple runs (see above). Please specify a run as date and subject (typically e.g. 171030ALAN)');
    runNum = input('Now specify a run number, typically a three-character string like 012');
    runID = (sprintf('%srun%s',dateSubj,runNum));
  else
    runID = dateSubjList{1};
  end
end
assert(isfield(db,runID),'database given does not contain data from the run requested');

splitField = strsplit(fieldName,'.');
assert(length(splitField) <=2, 'invalid fieldname %s: only one level of subfield currently implemented',fieldName);
fieldName = splitField{1};
if length(splitField) == 2
  subfieldName = splitField{2};
  fieldDepth = 2;
else
  fieldDepth = 1;
end

assert(isfield(db.(runID).fields,fieldName),'invalid field: %s',fieldName);
if fieldDepth == 2
  assert(isfield(db.(runID).fields.(fieldName),subfieldName),'invalid subfield: %s',subfieldName);
end

for ind_i = 1:length(fieldValue)-1
  if fieldDepth == 2 && isfield(db.(runID).fields.(fieldName).(subfieldName).meta.(sprintf('axis%d',ind_i))) && exist('useRawIndices','var') && useRawIndices(ind_i)
    newInd = find(db.(runID).fields.(fieldName).(subfieldName).meta.(sprintf('axis%d',ind_i)) == fieldValue(ind_i));
    assert(~isempty(newInd),'queried value at point not measured'); 
    fieldValue(ind_i) = newInd;
  elseif isfield(db.(runID).fields.(fieldName).meta.(sprintf('axis%d',ind_i))) && exist('useRawIndices','var') && useRawIndices(ind_i)
    newInd = find(db.(runID).fields.(fieldName).meta.(sprintf('axis%d',ind_i)) == fieldValue(ind_i));
    assert(~isempty(newInd),'queried value at point not measured'); 
    fieldValue(ind_i) = newInd;
  end
end


if fieldDepth == 1
  fieldData = db.(runID).fields.(fieldName).data;
else
  fieldData = db.(runID).fields.(fieldName).(subfieldName).data;
end
switch class(db.(runID).fields.(fieldName).data)
  case 'cell'
    if ischar(fieldValue)
      assert(ischar(fieldData{1}),'querried with char, but values are %s',class(db.(runID).fields.(fieldName){1}));
      trialIDs = find(cellfun(@(s) strcmp(s,fieldValue),fieldData));
    elseif iscell(fieldValue)
      assert(false,'query by value not currently implemented for values of class cell array');
    elseif isa(fieldValue,'function_handle')
      assert(false,'query by value not currently implemented for values of class function_handle');
    else
      assert(false,'query by value not currently implemented for values of type %s in fields whose data is of type cell',class(fieldValue));
    end
  case 'double'  
    if numel(fieldValue) == 1
      assert(size(fieldData,2) == 1 && ismatrix(fieldData), 'error: invalid query shape for requested field shape');
      trialIDs = find(fieldData == fieldValue);
    else 
      %assert(size(db.(runID).fields.(fieldName).data,2) == 2 && ismatrix(db.(runID).fields.(fieldName).data), 'error: invalid query shape for requested field shape');
      trialIDs = find(squeeze(fieldData(:,fieldValue(1:end-1))) == fieldValue(end));
    end
  case 'struct'
%     subfieldNames = fieldnames(db.(runID).fields.(fieldName));
%     subfieldName = subfieldNames{1};
%     trialIDs = find(structfun(@ (x) any(isequal(x,fieldValue),db.(runID).fields.(fieldName).(subfieldName))));
    trialIDs = [];
  otherwise
    trialIDs = [];
end
end

