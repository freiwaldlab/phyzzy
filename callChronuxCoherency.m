function [ varargout ] = callChronuxCoherency( funcHandle, varargin)
%callChronuxCoherency runs Chronux coherency functions such that nargout is
%independent of the error calculation type, as specified in params.err(1)
%   Chronux demands that the number of output arguments match the error
%   computation type. This means that you need to write out a separate call
%   for every error computation type. This function removes that need. 
%   Inputs: funcHandle - a function handle to one of:  
%           coherencyc, coherencycpt, coherencycpb, coherencypt, coherencypb, cohgramc, 
%           cohgramcpt, cohgramcpb, cohgrampt, cohgrampb
%   Returns: the longest possible set of chronux coherency outputs
%           - for field-field, [ C,phi,S12,S1,S2,f,confC,phistd,Cerr ]
%           - for all others, [ C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr ]
%           - uncomputed errors are set to NaN
%   Example: 
%           - [ C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr ] = callChronuxCoherency(@coherencycpb,lfpData,spikeData,chronuxParams)

funcString = func2str(funcHandle);
assert(strcmp(funcString,'coherencyc') || strcmp(funcString,'coherencycpt') || strcmp(funcString,'coherencycpb')|| ...
  strcmp(funcString,'cohgramc') || strcmp(funcString,'cohgramcpt') || strcmp(funcString,'cohgramcpb') || ...
  strcmp(funcString,'coherencypt') || strcmp(funcString,'coherencypb') || ...
  strcmp(funcString,'cohgrampt') || strcmp(funcString,'cohgrampb'), ...
  sprintf('Invalid function handle @%s. Acceptable values are: coherencyc, coherencycpt, coherencycpb, coherencypt, coherencypb, cohgramc, cohgramcpt, cohgramcpb, cohgrampt, cohgrampb',...
  funcString));
assert(isstruct(varargin{end}) && isfield(varargin{end},'err') && any(varargin{end}.err(1) == [0 1 2]), 'Invalid chronux params');

if ~isempty(regexp(funcString,'gram','ONCE'))
  if ~isempty(regexp(funcString,'pt','ONCE')) || ~isempty(regexp(funcString,'pb','ONCE'))
    if varargin{end}.err(1) == 2
      [ C,phi,S12,S1,S2,t,f,zerosp,confC,phistd,Cerr ] = funcHandle(varargin{:});
      varargout = {C,phi,S12,S1,S2,t,f,zerosp,confC,phistd,Cerr};
      return
    end
    if varargin{end}.err(1) == 1
      [ C,phi,S12,S1,S2,t,f,zerosp,confC,phistd ] = funcHandle(varargin{:});
      Cerr = NaN;
      varargout = {C,phi,S12,S1,S2,t,f,zerosp,confC,phistd,Cerr};
      return 
    else
      [ C,phi,S12,S1,S2,t,f,zerosp ] = funcHandle(varargin{:});
      Cerr = NaN;
      confC = NaN;
      phistd = NaN;
      varargout = {C,phi,S12,S1,S2,t,f,zerosp,confC,phistd,Cerr};
      return 
    end
  else
    if varargin{end}.err(1) == 2
      [ C,phi,S12,S1,S2,t,f,confC,phistd,Cerr ] = funcHandle(varargin{:});
      varargout = {C,phi,S12,S1,S2,t,f,confC,phistd,Cerr};
      return
    end
    if varargin{end}.err(1) == 1
      [ C,phi,S12,S1,S2,t,f,confC,phistd ] = funcHandle(varargin{:});
      Cerr = NaN;
      varargout = {C,phi,S12,S1,S2,t,f,confC,phistd,Cerr};
      return 
    else
      [ C,phi,S12,S1,S2,t,f ] = funcHandle(varargin{:});
      Cerr = NaN;
      confC = NaN;
      phistd = NaN;
      varargout = {C,phi,S12,S1,S2,t,f,confC,phistd,Cerr};
      return 
    end
  end
else
  if ~isempty(regexp(funcString,'pt','ONCE')) || ~isempty(regexp(funcString,'pb','ONCE'))
    if varargin{end}.err(1) == 2
      [ C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr ] = funcHandle(varargin{:});
      varargout = {C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr};
      return
    end
    if varargin{end}.err(1) == 1
      [ C,phi,S12,S1,S2,f,zerosp,confC,phistd ] = funcHandle(varargin{:});
      Cerr = NaN;
      varargout = {C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr};
      return 
    else
      [ C,phi,S12,S1,S2,f,zerosp ] = funcHandle(varargin{:});
      Cerr = NaN;
      confC = NaN;
      phistd = NaN;
      varargout = {C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr};
      return 
    end
  else
    if varargin{end}.err(1) == 2
      [ C,phi,S12,S1,S2,f,confC,phistd,Cerr ] = funcHandle(varargin{:});
      varargout = {C,phi,S12,S1,S2,f,confC,phistd,Cerr};
      return
    end
    if varargin{end}.err(1) == 1
      [ C,phi,S12,S1,S2,f,confC,phistd ] = funcHandle(varargin{:});
      Cerr = NaN;
      varargout = {C,phi,S12,S1,S2,f,confC,phistd,Cerr};
      return 
    else
      [ C,phi,S12,S1,S2,f ] = funcHandle(varargin{:});
      Cerr = NaN;
      confC = NaN;
      phistd = NaN;
      varargout = {C,phi,S12,S1,S2,f,confC,phistd,Cerr};
      return 
    end
  end
end
end

