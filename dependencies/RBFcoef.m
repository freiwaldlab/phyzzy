%   RBFcoef.m
%   =========
%   find radial basis function coefficients for interpolation
%   done by RBFipol.m
%
%   This function takes a set of N coordinates in DIM dimensions,
%   and a scalar function f(N) defined over that set of sample coordinates.
%   It returns the coefficient lambda(N) necessary to do a Radial Basis
%   Function (RBF) interpolation. The width sigm of the RBF should roughly match
%   the typical distance between the data points.
%   The (optional) output arguments return the maximum and minimum distances
%   between each data point. The (optional) output arguments return
%   the regularized data points, that means, the original data points minus
%   data points which have the same sample coordinates.
%
%   Usage: [lambda, mindists, maxdists, x_regular, f_regular] = ...
%                   RBFcoef(f, x, sigm, flag_regular);
%       x(DIM, N)              : sample coordinates (data points)
%       f(N, 1)                : sample values
%       sigm                   : global width of Gaussians
%       flag_regular           : eliminate double data points by averaging
%                                (OPTIONAL, default=0)
%
%       lambda(N, 1)           : coefficients for RBFipol.m
%       mindists(1, N)         : minimum distances to nearest neighbors
%       maxdists(1, N)         : maximum distances to next elements
%

function [lambda, mindists, maxdists, x_regular, f_regular] = ...
	RBFcoef(f, x, sigm, flag_regular)


% process optional arguments
if ~exist('flag_regular') flag_regular=0; end


% get size and length information
s = size(x);
dim = s(1);
n = s(2);
s = size(f);
m = s(1);

if ((m ~= n) | (s(2) ~= 1))
	error('Function must match vectors, f=f(x)');
end


% calculate the squared euklidean distances
d = zeros([n n]);
idx = 1:n;
one = ones([1 n]);
fprintf('Calculating %i distances: ', n);
for i=idx
	fprintf('.');
	d(i, :) = sum((x(:, idx)-x(:, i*one)).^2, 1);
end
fprintf('\n');


% omit the case of zero distances (ambiguous functional values)
if (flag_regular == 1)
	fprintf('Searching for zero distances: ');
	flags = ones([1 n]);
	for i=idx
		fprintf('.');
		if (flags(i) == 1)
			j = find(d(i, [1:i-1 i+1:n]) == 0);
			if (length(j) > 0)
				idcs = j+(j >= i);
				flags(idcs) = 0;
				f(i) = mean(f([i idcs]));
			end
		end
	end
	idcs_extract = find(flags);
	f_regular = f(idcs_extract, 1);
	x_regular = x(:, idcs_extract);
	d = d(idcs_extract, idcs_extract);
	n = length(f_regular);
	idx = 1:n;
	fprintf('\n');
else
	x_regular = x;
	f_regular = f;
end

% if optional output arg is given, also return minimum/maximum distances
if (nargout > 1)
	fprintf('Calculating minimum/maximum distances...\n');
	mindists = zeros([1 n]);
	maxdists = zeros([1 n]);
	for i=idx
		fprintf('.');
		mindists(i) = min(d(i, [1:i-1 i+1:n]));
		maxdists(i) = max(d(i, [1:i-1 i+1:n]));
	end
	fprintf('\n');
end
fprintf('\n');

% finally (and this is the heart of the algorithm), return coefficients
norm = 1/(0.5*sigm.^2);
d = exp(-norm*d);
fprintf('Doing matrix inversion...\n');
lambda = d\f_regular;
