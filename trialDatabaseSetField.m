function [ db ] = trialDatabaseSetField( fieldName, fieldValue, db, trialID, varargin )
%trialDatabaseSetField adds a field-value pair to the trial database
%   Places the key value pair into a trial database
%   - if the database contains data from multiple days, subjects, or runs, must specify dateSubj and runNum to which you are adding a field/value
%   - maintains setBy.function and setBy.line fields for each field-trial pair
%   - saves the commit hash and git status by field, at the time of the most recent update
%   Params:
%   - fieldName: a string; must be either a valid struct field name, or of
%                the form s1.s2, where s1 and s2 are valid struct field
%                names.
%   - fieldValue: acceptable types are string, numeric array, cell array, and one-field struct
%   - db: the trial database struct; initialize wtih trialDatabaseInit
%   - trialID: a positive integer; the unique trialID within the run, assigned in processRun 
%   - varargin: name-value pairs
%     - dateSubj: a string
%     - runNum: a string
%     - metaFieldname: name of meta data field; in this case, trialID must
%           be 'all'; special handling in trialDatabaseQuery when
%           metaFieldname is axis1, axis2, etc.
%

assert(mod(length(varargin),2) == 0, 'trialDatabaseSetField takes optional arguments as name-value pairs; odd number provided');
for argPair_i=1:length(varargin)/2
  argName = varargin{1+2*(argPair_i-1)};
  argVal = varargin{2+2*(argPair_i-1)};
  if strcmp(argName,'dateSubj')
    dateSubj = argVal;
  elseif strcmp(argName,'runNum')
    runNum = argVal;
  elseif strcmp(argName,'metaFieldname')
    metaFieldname = argVal;
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

if exist('metaFieldname','var')
  assert(ischar(trialID) && strcmp(trialID,'all'),'trial datebase insertion of trial-wise metadata not implemented; for metadata insertion, trialID must be all');
  db.(runID).fields.(metaFieldname) = fieldValue;
  return
end

splitField = strsplit(fieldName,'.');
assert(length(splitField) <=2, 'invalid fieldname %s: only one level of subfield currently implemented',fieldName);
fieldName = splitField{1};
if length(splitField) == 2
  subfieldName = splitField{2};
  fieldDepth = 2;
else
  fieldDepth = 1;
end

if fieldDepth == 1
  newField = ~isfield(db.(runID).fields,fieldName);
else
  newField = ~isfield(db.(runID).fields,fieldName) || ~isfield(db.(runID).fields.(fieldName),subfieldName);
end

if newField
  system('git log -1 > tmp.txt');
  tmp = textscan(fopen('tmp.txt','r'),'%q');
  commitHash = tmp{1}{2};
  system('git status > tmp.txt');
  tmp = textscan(fopen('tmp.txt','r'),'%q');
  fclose('all');
  delete('tmp.txt');
  gitStatus = tmp{1};
  if fieldDepth == 1
    db.(runID).fields.(fieldName).commit = commitHash;
    db.(runID).fields.(fieldName).gitStatus = gitStatus; %todo: format better
  else
    db.(runID).fields.(fieldName).(subfieldName).commit = commitHash;
    db.(runID).fields.(fieldName).(subfieldName).gitStatus = gitStatus; %todo: format better
  end
end

if fieldDepth == 1
  newDataField = ~isfield(db.(runID).fields.(fieldName),'data');
else
  newDataField = ~isfield(db.(runID).fields.(fieldName).(subfieldName),'data');
end
if newDataField
  if fieldDepth == 1
    db.(runID).fields.(fieldName).setBy.function = cell(db.(runID).numTrials,1);
    db.(runID).fields.(fieldName).setBy.line = NaN(db.(runID).numTrials,1);
  else
    db.(runID).fields.(fieldName).(subfieldName).setBy.function = cell(db.(runID).numTrials,1);
    db.(runID).fields.(fieldName).(subfieldName).setBy.line = NaN(db.(runID).numTrials,1);
  end
end
callStack = dbstack();
if fieldDepth == 1
  db.(runID).fields.(fieldName).setBy.function{trialID} = callStack(2).name;
  db.(runID).fields.(fieldName).setBy.line(trialID) = callStack(2).line;
else
  db.(runID).fields.(fieldName).(subfieldName).setBy.function{trialID} = callStack(2).name;
  db.(runID).fields.(fieldName).(subfieldName).setBy.line(trialID) = callStack(2).line;
end
switch class(fieldValue)
  case 'char'
    if newDataField
      if fieldDepth == 1
        db.(runID).fields.(fieldName).data = cell(db.(runID).numTrials,1);
      else
        db.(runID).fields.(fieldName).(subfieldName).data = cell(db.(runID).numTrials,1);
      end
    end
    if fieldDepth == 1
      db.(runID).fields.(fieldName).data{trialID} = fieldValue;
    else
      db.(runID).fields.(fieldName).(subfieldName).data{trialID} = fieldValue;
    end
  case 'double'
        
    if newDataField
      if numel(fieldValue) == 1  
        if fieldDepth == 1
          db.(runID).fields.(fieldName).data = NaN(db.(runID).numTrials,1);
        else
          db.(runID).fields.(fieldName).(subfieldName).data = NaN(db.(runID).numTrials,1);
        end
      elseif length(size(fieldValue)) == 2 && min(size(fieldValue)) == 1  %for 1xN or Nx1 fieldValue, make db.(runID).numTrials x N array
        if fieldDepth == 1
          db.(runID).fields.(fieldName).data = NaN(db.(runID).numTrials,length(fieldValue));
        else
          db.(runID).fields.(fieldName).(subfieldName).data = NaN(db.(runID).numTrials,length(fieldValue));
        end
      else
        if fieldDepth == 1
          db.(runID).fields.(fieldName).data = NaN([db.(runID).numTrials,size(fieldValue)]);
        else
          db.(runID).fields.(fieldName).(subfieldName).data = NaN([db.(runID).numTrials,size(fieldValue)]);
        end
      end
    end
    
    if numel(fieldValue) == 1
      if fieldDepth == 1
        db.(runID).fields.(fieldName).data(trialID) = fieldValue;
      else
        db.(runID).fields.(fieldName).(subfieldName).data(trialID) = fieldValue;
      end
    elseif length(size(fieldValue)) == 2 && min(size(fieldValue)) == 1
      if fieldDepth == 1
        db.(runID).fields.(fieldName).data(trialID,:) = fieldValue;  %matlab handles the implicit reshape when needed
      else
        db.(runID).fields.(fieldName).(subfieldName).data(trialID,:) = fieldValue;  %matlab handles the implicit reshape when needed
      end
    else
      if fieldDepth == 1
        db.(runID).fields.(fieldName).data(trialID,:) = fieldValue(:);
      else
        db.(runID).fields.(fieldName).(subfieldName).data(trialID,:) = fieldValue(:);
      end
    end
  case 'cell'
    if newDataField
      if fieldDepth == 1
        db.(runID).fields.(fieldName).data = cell(db.(runID).numTrials,1);
      else
        db.(runID).fields.(fieldName).(subfieldName).data = cell(db.(runID).numTrials,1);
      end
    end
    if fieldDepth == 1
      db.(runID).fields.(fieldName).data{trialID} = fieldValue;
    else
      db.(runID).fields.(fieldName).(subfieldName).data{trialID} = fieldValue;
    end
      
  case 'struct'
    assert(length(fieldnames(fieldValue)) == 1, 'Error: database insertion not implemented for multi-field structs');
    if newDataField
      if fieldDepth == 1
        db.(runID).fields.(fieldName).data = repmat(struct(fieldnames(fieldValue),NaN),db.(runID).numTrials,1);
      else
        db.(runID).fields.(fieldName).(subfieldName).data = repmat(struct(fieldnames(fieldValue),NaN),db.(runID).numTrials,1);
      end
    end
    if fieldDepth == 1
      db.(runID).fields.(fieldName).data(trialID) = fieldValue;
    else
      db.(runID).fields.(fieldName).(subfieldName).data(trialID) = fieldValue;
    end
      
  case 'function_handle'
    if newDataField
      if fieldDepth == 1
        db.(runID).fields.(fieldName).data = cell(db.(runID).numTrials,1);  %note: matlab doesn't allow non-cell arrays of function handles
      else
        db.(runID).fields.(fieldName).(subfieldName).data = cell(db.(runID).numTrials,1);  %note: matlab doesn't allow non-cell arrays of function handles
      end
    end
    if fieldDepth == 1
      db.(runID).fields.(fieldName).data{trialID} = fieldValue;
    else
      db.(runID).fields.(fieldName).(subfieldName).data{trialID} = fieldValue;
    end
  otherwise
    error('Database insertion not implemented for values of the type given');
end
end

