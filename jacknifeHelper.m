function [m,jsd,conf]=jacknifeHelper(x,f, varargin)
% Compute jackknife estimates of the mean and standard deviation of input
% data x, run through functional handle f, with n draws of m samples
% 
% Inputs:
%   - x : data in the form trials (nonsingleton) x values (possibly singleton)
%   - f : function handle (for simple jacknife mean estimate, use f = @mean)
%         Note: f must return a 1-d (possibly singleton) array
%   - varargin: valid name value pairs are:
%       - 'draws', int: default = size(x,1)
%       - 'trialsPerDraw', int, default = size(x,1)-1
%       - 'confPercentile', float in [0,100]: error bars will be +/- floor;
%          if confPercentile == 0, does not compute percentile cutoffs
%         of the rank converted from this perentile
%       - CURRENTLY NOT IMPLEMENTED: 'parallelWorkers',
% Outputs:
% m : estimate of the mean (across trials)
% jsd: jackknife estimate of the standard deviation (across trials)
% Notes:
%   - if trialsPerDraw == numTrials-1 (default), samples without
%     replacement. Otherwise, samples with replacement across draws.

[trials,vars] = size(x);
draws = trials;
trialsPerDraw = trials-1;
confPercentile = NaN;
parallelWorkers = 0;
separableF = 0;
auxiliaryFunctions = {};
assert(trials > 1,'Need multiple trials');
assert(mod(length(varargin),2) == 0, 'varargin for jacknifeHelper must consist of name-value pairs');
for i = 1:length(varargin)/2
  switch varargin{2*i-1}
    case 'draws'
      draws = min(varargin{2*i},draws);
    case 'trialsPerDraw'
      trialsPerDraw = min(varargin{2*i},trialsPerDraw);
    case 'confPercentile'
      confPercentile = varargin{2*i};
      assert(confPercentile >= 0 && confPercentile <= 100, 'confidence interval percentile must be between 0 and 100');
    case 'parallelWorkers'
      parallelWorkers = varargin{2*i};
      if parallelWorkers > 0
        disp('parallel jacknifeHelper not currently implemented; using serial');
      end
    case 'separableF'
      separableF = varargin{2*i};
    case 'auxiliaryFunction'
      auxiliaryFunction = varargin{2*i};
  end
end
m = f(x);
assert(numel(m) == max(size(m)), 'jacknifeHelper cannot handle greater than 1D outputs');
assert(trialsPerDraw == trials-1 || ~separableF,'undefined behavior when trialsPerDraw ~= trials-1 and separableF == 1');
theta = zeros(draws,length(m));
estimates = zeros(draws,length(m));
if separableF
  if isempty(auxiliaryFunctions)
    for draw = 1:draws
      y = m - f(x(draw,:));
      estimates(draw,:) = y;
      theta(draw,:) = (trialsPerDraw+1)*m-trialsPerDraw*y; % pseudo values
    end
  else
    mAuxVars = auxiliaryFunction(x);
    for draw = 1:draws
      y = f(mAuxVars - auxiliaryFunction(x(draw,:)));
      estimates(draw,:) = y;
      theta(draw,:) = (trialsPerDraw+1)*m-trialsPerDraw*y; % pseudo values
    end
  end
else
  if trialsPerDraw == trials-1
      idx = 2:trials; % drop 1 trial
  end
  %parfor (draw = 1:draws, parallelWorkers)
  for draw = 1:draws
    if trialsPerDraw == trials-1 && draw > 1
      idx(draw-1) = draw-1; % drop 1 trial
    else
      idx = randperm(trials,trialsPerDraw);
    end
    y = f(x(idx,:));
    estimates(draw,:) = y;
    theta(draw,:) = (trialsPerDraw+1)*m-trialsPerDraw*y; % pseudo values\
  end
end
jm = mean(theta,1);
jsd = sqrt(sum((theta-jm).^2,1)/(draws*(draws-1)));
if confPercentile > 0
  conf = zeros(size(m,2),2);
  confWinInd = ceil(trials*(1-.01*confPercentile/2));
  for var_i = 1:length(m)
    sortedValues = sort(estimates(:,var_i),'ascend');
    conf(var_i,:) = [sortedValues(confWinInd,var_i)-m, sortedValues(end-confWinInd+1,var_i)-m];
  end
else
  conf = NaN(size(m,2),2);
end
end
